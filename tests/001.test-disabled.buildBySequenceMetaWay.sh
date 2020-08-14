#!/bin/bash

printf "Testing buildBySequenceMeta... "
g++ -std=c++0x -I $GEMSHOME/gmml/includes/ -L$GEMSHOME/gmml/bin/ -Wl,-rpath,$GEMSHOME/gmml/bin/ tests/buildBySequenceMeta.cc -lgmml -o buildBySequenceMeta
./buildBySequenceMeta > /dev/null 2>&1
if [ -f default.pdb ]; then
    if ! cmp default.pdb tests/correct_outputs/buildBySequenceMeta.pdb > /dev/null 2>&1; then
        printf "Test FAILED!.\n"
        return 1;
    else
        printf "Test passed.\n"
        rm default.pdb buildBySequenceMeta
        return 0;
    fi
else
    printf "Test FAILED!.\n"
    return 1;
fi

