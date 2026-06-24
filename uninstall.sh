#!/bin/bash
cd "${0%/*}" || exit
#------------------------------------------------------------------------------

CAPTIONFOAM_BASHRC="$PWD/etc/bashrc"

# Remove from ~/.bashrc
if grep -q "$CAPTIONFOAM_BASHRC" "$HOME/.bashrc"; then
    sed -i "\|# captionFOAM environment|d" "$HOME/.bashrc"
    sed -i "\|$CAPTIONFOAM_BASHRC|d" "$HOME/.bashrc"
    echo "Removed captionFOAM from ~/.bashrc"
else
    echo "captionFOAM not found in ~/.bashrc"
fi

# Clean libraries and solvers
wclean all src
wclean all applications/solvers

# Unset environment variables
unset CAPTIONFOAM_DIR
unset CAPTIONFOAM_SRC

#------------------------------------------------------------------------------