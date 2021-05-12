#!/bin/bash

printf "Testing buildBySequenceRotamer... "
g++ -std=c++0x -I $GEMSHOME/gmml/ -I $GEMSHOME/gmml/ -L$GEMSHOME/gmml/bin/ -Wl,-rpath,$GEMSHOME/gmml/bin/ tests/010.buildBySequenceRotamer.cc -lgmml -pthread -o buildBySequenceRotamer
./buildBySequenceRotamer > 010.output_buildBySequenceRotamer.txt
if [ -f structure.pdb ] && [ -f structure.off ]; then
    if ! cmp structure.pdb tests/correct_outputs/010.buildBySequenceRotamer.pdb > /dev/null 2>&1; then
        printf "Test FAILED!. PDB file different\n"
        return 1;
    #elif ! cmp structure.off tests/correct_outputs/010.buildBySequenceRotamer.off > /dev/null 2>&1; then
    #    printf "Test FAILED!. Off file different.\n"
    #    return 1;
    elif ! cmp  010.output_buildBySequenceRotamer.txt tests/correct_outputs/010.output_buildBySequenceRotamer.txt > /dev/null 2>&1; then
        printf "Test FAILED!. Output file different\n"
        return 1;
    else
        printf "Test passed.\n"
        rm structure.pdb structure.off buildBySequenceRotamer 010.output_buildBySequenceRotamer.txt
        return 0;
    fi
else
    printf "Test FAILED!.\n"
    return 1;
fi

