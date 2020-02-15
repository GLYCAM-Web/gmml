#!/bin/bash

printf "Testing pdb2glycam and molecule subgraph matching... "
g++ -std=c++0x -I $GEMSHOME/gmml/includes/ -L$GEMSHOME/gmml/bin/ -Wl,-rpath,$GEMSHOME/gmml/bin/ tests/pdb2glycam.cc -lgmml -o pdb2glycam
./pdb2glycam tests/inputs/pdb2glycam_4YG0.pdb
if [ -f pdb2glycam_output.pdb ]; then
    if ! cmp pdb2glycam_output.pdb tests/correct_outputs/pdb2glycam_4YG0_output.pdb > /dev/null 2>&1; then
        printf "Test FAILED!.\n"
        return 1
    else
        printf "Test passed.\n"
        rm ring_conformations.txt pdb2glycam_output.pdb pdb2glycam > /dev/null 2>&1
        return 0
    fi
else
    printf "Test FAILED!.\n"
    return 1
fi
