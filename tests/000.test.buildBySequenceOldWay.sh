#!/bin/bash

printf "Testing buildBySequence... "
g++ -std=c++0x -I $GEMSHOME/gmml/includes/ -L$GEMSHOME/gmml/bin/ -Wl,-rpath,$GEMSHOME/gmml/bin/ tests/buildBySequence.cc -lgmml -o buildBySequence
./buildBySequence > /dev/null 2>&1
#./buildBySequence
if [ -f buildBySequence.pdb ]; then
    if ! cmp buildBySequence.pdb tests/correct_outputs/buildBySequence.pdb > /dev/null 2>&1; then
        printf "Test FAILED!.\n"
        return 1;
    else
        printf "Test passed.\n"
        rm buildBySequence buildBySequence.pdb
        return 0;
    fi
else
    printf "Test FAILED!.\n"
    return 1;
fi
