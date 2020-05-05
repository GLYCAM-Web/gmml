#!/bin/bash

printf "Testing BFMP Ring Shape Calculation... "
g++ -std=c++0x -I $GEMSHOME/gmml/includes/ -L$GEMSHOME/gmml/bin/ -Wl,-rpath,$GEMSHOME/gmml/bin/ tests/ring_shape_detection.cc -lgmml -o ring_shape_detection
./ring_shape_detection > ring_shape_detection.txt
if [ -f ring_conformations.txt ]; then
    if ! cmp ring_conformations.txt tests/correct_outputs/ring_conformations.txt > /dev/null 2>&1; then
        printf "Test FAILED!.\nPlease compare ring_conformations.txt to tests/correct_outputs/ring_conformations.txt"
        return 1
    else
        printf "Test passed.\n"
        rm ring_shape_detection ring_shape_detection.txt ring_conformations.txt > /dev/null 2>&1
        return 0
    fi
else
    printf "Test FAILED!.\nPlease compare ring_conformations.txt to tests/correct_outputs/ring_conformations.txt"
    return 1
fi
