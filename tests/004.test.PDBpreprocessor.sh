#!/bin/bash

printf "Testing PDBPreprocessor... "
g++ -std=c++0x -I $GEMSHOME/gmml/ -L$GEMSHOME/gmml/bin/ -Wl,-rpath,$GEMSHOME/gmml/bin/ tests/004.PDBPreprocessor.cc -lgmml -pthread -o PDBPreprocessor
./PDBPreprocessor > /dev/null 2>&1
if [ -f Processed.pdb ]; then
    if ! cmp Processed.pdb tests/correct_outputs/Processed.pdb > /dev/null 2>&1; then
        printf "Test FAILED!\n"
        echo "Exit Code: 1"
        return 1
    else
        printf "Test passed.\n"
        rm Processed.pdb PDBPreprocessor > /dev/null 2>&1
        echo "Exit Code: 0"
        return 0
    fi
else
    printf "Test FAILED!\n"
    echo "Exit Code: 1"
    return 1
fi
