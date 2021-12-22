#!/bin/bash

printf "Testing PDBPreprocessor... "
g++ -std=c++17 -I $GEMSHOME/gmml/ -L$GEMSHOME/gmml/bin/ -Wl,-rpath,$GEMSHOME/gmml/bin/ tests/004.PDBPreprocessor.cc -lgmml -pthread -o PDBPreprocessor 
./PDBPreprocessor tests/inputs/preprocessor_input.pdb > /dev/null 2>&1
if [ -f Processed.pdb ]; then
    if ! cmp Processed.pdb tests/correct_outputs/Processed.pdb > /dev/null 2>&1; then
        printf "Test FAILED!.\nINFO: These files are different: Processed.pdb tests/correct_outputs/Processed.pdb\n"
        return 1
    else
        printf "Test passed.\n"
        rm Processed.pdb PDBPreprocessor > /dev/null 2>&1
        return 0
    fi
else
    printf "Test FAILED!.\nProcessed.pdb was not created.\n"
    return 1
fi
