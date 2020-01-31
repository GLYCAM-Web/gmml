#!/bin/bash

printf "Testing superimposition_Eigen... "
g++ -std=c++0x -I $GEMSHOME/gmml/includes/ -L$GEMSHOME/gmml/bin/ -Wl,-rpath,$GEMSHOME/gmml/bin/ tests/superimposition_Eigen.cc -lgmml -o superimposition_Eigen
./superimposition_Eigen > /dev/null 2>&1
if [ -f moved.pdb ]; then
    if ! cmp moved.pdb tests/correct_outputs/moved.pdb > /dev/null 2>&1; then
        printf "Test FAILED!.\n"
        return 1
    else
        printf "Test passed.\n"
        rm moved.pdb superimposition_Eigen > /dev/null 2>&1
        return 0
    fi
else
    printf "Test FAILED!.\n"
    return 1
fi

