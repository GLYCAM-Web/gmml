#!/bin/bash

printf "Testing buildBySequenceRotamer... "
g++ -std=c++0x -I $GEMSHOME/gmml/includes/ -L$GEMSHOME/gmml/bin/ -Wl,-rpath,$GEMSHOME/gmml/bin/ tests/010.buildBySequenceRotamer.cc -lgmml -o buildBySequenceRotamer
./buildBySequenceRotamer > /dev/null 2>&1
if [ -f output.pdb ] && [ -f output.off ]; then
    if ! cmp output.pdb tests/correct_outputs/010.buildBySequenceRotamer.pdb > /dev/null 2>&1; then
        printf "Test FAILED!. PDB file different\n"
        return 1;
    elif ! cmp output.off tests/correct_outputs/010.buildBySequenceRotamer.off > /dev/null 2>&1; then
        printf "Test FAILED!. Off file different.\n"
        return 1;
    else
        printf "Test passed.\n"
        rm output.pdb output.off buildBySequenceRotamer
        return 0;
    fi
else
    printf "Test FAILED!.\n"
    return 1;
fi

