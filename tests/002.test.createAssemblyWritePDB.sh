#!/bin/bash

printf "Testing create_Assembly_WritePDB... "
g++ -std=c++0x -I $GEMSHOME/gmml/includes/ -L$GEMSHOME/gmml/bin/ -Wl,-rpath,$GEMSHOME/gmml/bin/ tests/create_Assembly_WritePDB.cc -lgmml -o create_Assembly_WritePDB
./create_Assembly_WritePDB > /dev/null 2>&1
if [ -f test-NLN.pdb ]; then
    if ! cmp test-NLN.pdb tests/correct_outputs/test-NLN.pdb > /dev/null 2>&1; then
        printf "Test FAILED!.\n"
        return 1
    else
        printf "Test passed.\n"
        rm test-NLN.pdb create_Assembly_WritePDB > /dev/null 2>&1
        return 0
    fi
else
    printf "Test FAILED!.\n"
    return 1
fi

