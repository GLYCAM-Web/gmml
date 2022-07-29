#!/bin/bash

printf "Testing 008.pdb2glycam.cc and molecule subgraph matching... "
g++ -std=c++17 -I $GEMSHOME/gmml/ -L$GEMSHOME/gmml/bin/ -Wl,-rpath,$GEMSHOME/gmml/bin/ tests/008.pdb2glycam.cc -lgmml -pthread -o pdb2glycam
./pdb2glycam tests/inputs/pdb2glycam_4YG0.pdb
if [ -f pdb2glycam_output.pdb ]; then
    if ! cmp pdb2glycam_output.pdb tests/correct_outputs/pdb2glycam_4YG0_output.pdb > /dev/null 2>&1; then
        printf "\nTest FAILED!\n pdb2glycam_output.pdb does not match tests/correct_outputs/pdb2glycam_4YG0_output.pdb"
        echo "Exit Code: 1"
        return 1
    else
        printf "Test passed.\n"
        rm ring_conformations.txt pdb2glycam_output.pdb pdb2glycam > /dev/null 2>&1
        echo "Exit Code: 0"
        return 0
    fi
else
    printf "\nTest FAILED!\n"
    echo "Exit Code: 1"
    return 1
fi
