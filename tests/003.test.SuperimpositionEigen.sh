#!/bin/bash

printf "Testing 003.superimpositionEigen.cc... "
g++ -std=c++17 -I $GEMSHOME/gmml/ -L$GEMSHOME/gmml/bin/ -Wl,-rpath,$GEMSHOME/gmml/bin/ tests/003.superimpositionEigen.cc -lgmml -pthread -o superimposition_Eigen
./superimposition_Eigen > /dev/null 2>&1
if [ -f moved.pdb ]; then
    if ! cmp moved.pdb tests/correct_outputs/moved.pdb > /dev/null 2>&1; then
        printf "Test FAILED!\n"
        echo "Exit Code: 1"
        return 1
    else
        printf "Test passed\n"
        rm moved.pdb superimposition_Eigen > /dev/null 2>&1
        echo "Exit Code: 0"
        return 0
    fi
else
    printf "Test FAILED!\n"
    echo "Exit Code: 1"
    return 1
fi

