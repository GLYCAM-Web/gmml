#!/bin/bash

printf "Testing PDBpreprocessor... "
g++ -std=c++0x -I $GEMSHOME/gmml/includes/ -L$GEMSHOME/gmml/bin/ -Wl,-rpath,$GEMSHOME/gmml/bin/ tests/PDB_preprocessor.cc -lgmml -o PDB_preprocessor
./PDB_preprocessor > /dev/null 2>&1
if [ -f Processed.pdb ]; then
    if ! cmp Processed.pdb tests/correct_outputs/Processed.pdb > /dev/null 2>&1; then
        printf "Test FAILED!.\n"
        return 1
    else
        printf "Test passed.\n"
        rm Processed.pdb PDB_preprocessor > /dev/null 2>&1
        return 0
    fi
else
    printf "Test FAILED!.\n"
    return 1
fi
