// =======================================================================================
// Parameter definitions
// =======================================================================================

// Basic dimensions (m)
D_e = 0.125e-3;   // Electrode diameter
D_a = 20.0e-3;     // Outer domain diameter

L_e = 0.5e-3;     // Electrode length
L_a = 10.0e-3;    // Total domain length

thetaDgr = 1;     // Wedge angle in degrees
theta    = thetaDgr * Pi / 180;

// Radii
R_e = D_e / 2.0;
R_a = D_a / 2.0;

// Control arc radius — mesh is finer inside R_u, coarser outside.
// Must satisfy: R_e + buffR < R_u < R_a
R_u = 2.0e-3;

// Axial positions
Z_b = -L_e;           // Inlet
Z_e =  0;             // Electrode tip axial position
Z_a =  L_a - L_e;    // Outlet
Z_u =  2.0e-3;        // Axial control point — mesh is finer near electrode tip, coarser beyond

// =======================================================================================
// Structured buffer layer around electrode
// =======================================================================================
nLayers           = 21;
nPoints           = nLayers + 1;
layerRatio        = 1.0;
minLayerThickness = 2.0e-6;

bufferThickness = 0.0;
h_i = minLayerThickness;
For i In {1:nPoints}
  h_i = h_i * layerRatio;
  bufferThickness += h_i;
EndFor

buffR = bufferThickness;  // Radial thickness of the structured buffer

// =======================================================================================
// Mesh settings
// =======================================================================================
Mesh.Algorithm          = 5;  // Delaunay 2D
Mesh.RecombineAll       = 1;
Mesh.RecombinationAlgorithm = 0;

// Characteristic lengths
lc     = 0.2 * R_e;    // Fine  — on the buffer outer boundary
lc_u   = lc;           // Medium — on the control arc at R_u  (adjust to taste)
lc_zu  = lc;           // Fine  — on the axial control point at Z_u  (adjust to taste)
lc_top = 5*lc;         // Coarse — on the outer boundary at R_a

// Transfinite divisions for the structured buffer
nAxialEle    = 200;
nArc         = 70;
progAxialEle = 1.0;

// =======================================================================================
// Points
// =======================================================================================

// Centre of electrode tip arc
c1  = newp; Point(c1)  = {0,          0,         Z_e - R_e,  lc     };

// Electrode (inner wall) and buffer
p1  = newp; Point(p1)  = {0,  R_e,           Z_b,        lc     };
p2  = newp; Point(p2)  = {0,  R_e,           Z_e - R_e,  0.5*lc     };
p3  = newp; Point(p3)  = {0,  0,             Z_e,        lc     };
p4  = newp; Point(p4)  = {0,  0,             Z_e + buffR, 0.5*lc    };
p6  = newp; Point(p6)  = {0,  R_e + buffR,   Z_b,        lc     };
p7  = newp; Point(p7)  = {0,  R_e + buffR,   Z_e - R_e,  lc     };

// Axis at outlet
p5  = newp; Point(p5)  = {0,  0,             Z_a,        0.5*lc   };

// Axial control point on the axis (y=0, z=Z_u) — splits l_axis for denser mesh near tip
pZu = newp; Point(pZu) = {0,  0,             Z_u,        0.5*lc    };

// Control arc corners at r = R_u
pU1 = newp; Point(pU1) = {0,  R_u,           Z_b,        lc_u   };
pU2 = newp; Point(pU2) = {0,  R_u,           Z_a,        lc_u   };

// Outer boundary corners at r = R_a
p20 = newp; Point(p20) = {0,  R_a,           Z_b,        lc_top };
p21 = newp; Point(p21) = {0,  R_a,           Z_a,        lc_top };

// =======================================================================================
// Lines
// =======================================================================================

// --- Surface 1: structured buffer rectangle (electrode shaft) ---
l1 = newl; Line(l1) = {p1, p2};         // electrode shaft wall
l6 = newl; Line(l6) = {p2, p7};         // shared with s2
l7 = newl; Line(l7) = {p6, p7};         // buffer outer (radial direction)
l5 = newl; Line(l5) = {p1, p6};         // inlet face of buffer

// --- Surface 2: structured buffer arc (electrode tip) ---
l2 = newl; Circle(l2) = {p2, c1, p3};   // electrode tip arc wall
l3 = newl; Line(l3)   = {p3, p4};       // axis segment at tip
l8 = newl; Circle(l8) = {p4, c1, p7};   // buffer outer arc

// --- Shared boundary for unstructured regions ---
l_in_i  = newl; Line(l_in_i)  = {p6,  pU1};  // inlet   : buffer outer → R_u
l_ctrl  = newl; Line(l_ctrl)  = {pU1, pU2};  // control arc at r = R_u
l_out_i  = newl; Line(l_out_i)  = {pU2, p5};   // outlet  : R_u → axis
l_axis1  = newl; Line(l_axis1)  = {p5,  pZu};  // axis    : outlet → Z_u control point
l_axis2  = newl; Line(l_axis2)  = {pZu, p4};   // axis    : Z_u control point → buffer end
l_in_o  = newl; Line(l_in_o)  = {pU1, p20};  // inlet   : R_u → R_a
l_outer = newl; Line(l_outer) = {p20, p21};  // outer radial boundary
l_out_o = newl; Line(l_out_o) = {p21, pU2};  // outlet  : R_a → R_u

// =======================================================================================
// Surfaces
// =======================================================================================

// Surface 1 — structured buffer rectangle
loop1 = newll; Line Loop(loop1) = {l1, l6, -l7, -l5};
s1    = news;  Plane Surface(s1) = {loop1};

// Surface 2 — structured buffer arc at electrode tip
loop2 = newll; Line Loop(loop2) = {l2, l3, l8, -l6};
s2    = news;  Plane Surface(s2) = {loop2};

// Surface 3 — inner unstructured region (buffer outer → R_u, fine mesh)
loop3 = newll; Line Loop(loop3) = {l_in_i, l_ctrl, l_out_i, l_axis1, l_axis2, l8, -l7};
s3    = news;  Plane Surface(s3) = {loop3};

// Surface 4 — outer unstructured region (R_u → R_a, coarser mesh)
loop4 = newll; Line Loop(loop4) = {l_in_o, l_outer, l_out_o, -l_ctrl};
s4    = news;  Plane Surface(s4) = {loop4};

// =======================================================================================
// Transfinite (structured) mesh on the buffer
// =======================================================================================

// Surface 1
Transfinite Line {l1, l7} = nAxialEle Using Progression progAxialEle;
Transfinite Line {l5, l6} = nPoints   Using Progression layerRatio;
Transfinite Surface {s1};
Recombine   Surface {s1};

// Surface 2
Transfinite Line {l2, l8} = nArc    Using Progression 1;
Transfinite Line {l3}     = nPoints Using Progression layerRatio;
Transfinite Surface {s2};
Recombine   Surface {s2};

// =======================================================================================
// 2D → 3D wedge extrusion
// =======================================================================================
dz = 1e-4;  // slab thickness in x

v1[] = Extrude {dz, 0, 0}{Surface{s1}; Layers{1}; Recombine;};
v2[] = Extrude {dz, 0, 0}{Surface{s2}; Layers{1}; Recombine;};
v3[] = Extrude {dz, 0, 0}{Surface{s3}; Layers{1}; Recombine;};
v4[] = Extrude {dz, 0, 0}{Surface{s4}; Layers{1}; Recombine;};

// =======================================================================================
// OpenFOAM boundary patches
// Lateral-face index mapping (Extrude returns [0]=back, [1]=vol, [2..]=curves in loop order)
//   s1 loop {l1, l6, -l7, -l5}  → v1[2]=electrode, v1[3]=internal, v1[4]=internal, v1[5]=inlet
//   s2 loop {l2, l3, l8, -l6}   → v2[2]=electrode, v2[3]=axis-tip,  v2[4]=internal, v2[5]=internal
//   s3 loop {l_in_i, l_ctrl, l_out_i, l_axis1, l_axis2, l8, -l7}
//                                → v3[2]=inlet-i,  v3[3]=internal,  v3[4]=outlet-i, v3[5]=axis1, v3[6]=axis2
//   s4 loop {l_in_o, l_outer, l_out_o, -l_ctrl}
//                                → v4[2]=inlet-o,  v4[3]=topBound,  v4[4]=outlet-o, v4[5]=internal
// =======================================================================================
Physical Surface("frontWedge")  = {s1, s2, s3, s4};
Physical Surface("backWedge")   = {v1[0], v2[0], v3[0], v4[0]};
Physical Surface("electrode")   = {v1[2], v2[2]};
Physical Surface("inlet")       = {v1[5], v3[2], v4[2]};
Physical Surface("outlet")      = {v3[4], v4[4]};
Physical Surface("axis")        = {v2[3], v3[5], v3[6]};
Physical Surface("topBoundary") = {v4[3]};

Physical Volume("region0") = {v1[1], v2[1], v3[1], v4[1]};
