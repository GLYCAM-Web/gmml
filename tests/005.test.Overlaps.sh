#!/bin/bash

printf "Testing Overlaps function... "
g++ -std=c++0x -I $GEMSHOME/gmml/includes/ -L$GEMSHOME/gmml/bin/ -Wl,-rpath,$GEMSHOME/gmml/bin/ tests/overlaps.cc -lgmml -o overlaps
./overlaps > overlaps.txt 2> /dev/null
if [ -f overlaps.txt ]; then
    if ! cmp overlaps.txt tests/correct_outputs/overlaps.txt > /dev/null 2>&1; then
        printf "Test FAILED!.\n"
        return 1
    else
        printf "Test passed.\n"
        rm overlaps overlaps.txt > /dev/null 2>&1
        return 0
    fi
else
    printf "Test FAILED!.\n"
    return 1
fi
