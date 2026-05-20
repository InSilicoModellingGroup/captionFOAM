// =======================================================================================
// Parameter definitions
// =======================================================================================
// Basic dimensions (m)
scale = 1.0e-3;

D_e = scale*0.125;
h_e = scale*0.5;
L_e = scale*1.0;

R_e = D_e / 2.0;

X0 = -h_e/2;
Y0 = 0;
Z0 = 0;

lc = 0.1*R_e;
w = 0.1;

// =======================================================================================
// Points, Lines, Loops, Surfaces, Volumes
// =======================================================================================
p1 = newp; Point(p1) = {X0, Y0, Z0-R_e, lc};
p2 = newp; Point(p2) = {X0, Y0, Z0, w*lc};
p3 = newp; Point(p3) = {X0, Y0+R_e, Z0-R_e, w*lc};
p4 = newp; Point(p4) = {X0, Y0+R_e, Z0-L_e, lc};
p5 = newp; Point(p5) = {X0, Y0, Z0-L_e, lc};

l1 = newl; Line(l1) = {p1, p2};
l2 = newl; Circle(l2) = {p2, p1, p3};
l3 = newl; Line(l3) = {p3, p4};
l4 = newl; Line(l4) = {p4, p5};
l5 = newl; Line(l5) = {p5, p1};
loop1 = newll; Line Loop(loop1) = {l1, l2, l3, l4, l5};
s1 = news; Plane Surface(s1) = {loop1};

v1[] = Extrude {h_e, 0, 0} {Surface{s1}; Layers{50}; };
