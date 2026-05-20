#!/bin/bash
cd "${0%/*}" || exit
#------------------------------------------------------------------------------

CAPTIONFOAM_BASHRC="$PWD/etc/bashrc"

# Check that etc/bashrc exists
if [ ! -f "$CAPTIONFOAM_BASHRC" ]; then
    echo "Error: $CAPTIONFOAM_BASHRC not found"
    exit 1
fi

# Check that OpenFOAM is sourced
if [ -z "$WM_PROJECT_DIR" ]; then
    echo "Error: OpenFOAM environment not found. Source OpenFOAM's etc/bashrc first."
    exit 1
fi

# Source captionFOAM environment
source "$CAPTIONFOAM_BASHRC"

# Add to ~/.bashrc if not already there
if ! grep -q "$CAPTIONFOAM_BASHRC" "$HOME/.bashrc"; then
    echo "" >> "$HOME/.bashrc"
    echo "# captionFOAM environment" >> "$HOME/.bashrc"
    echo "source $CAPTIONFOAM_BASHRC" >> "$HOME/.bashrc"
    echo "Added captionFOAM to ~/.bashrc"
else
    echo "captionFOAM already in ~/.bashrc"
fi

# Compile libraries and solvers
wmake all src
wmake all applications/solvers

#------------------------------------------------------------------------------