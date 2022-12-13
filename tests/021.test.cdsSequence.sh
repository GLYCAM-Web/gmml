#!/bin/bash

printf "Testing 021.test.cdsSequence.cpp... "
g++ -std=c++17 -I $GEMSHOME/gmml/ -I $GEMSHOME/gmml/ -L$GEMSHOME/gmml/bin/ -Wl,-rpath,$GEMSHOME/gmml/bin/ tests/021.test.cdsSequence.cpp -lgmml -pthread -o cdsSequence
./cdsSequence > 021.output_cdsSequence.txt 2>&1
if ! cmp  021.output_cdsSequence.txt tests/correct_outputs/021.output_cdsSequence.txt > /dev/null 2>&1; then
    printf "Test FAILED! 021.output_cdsSequence.txt different from tests/correct_outputs/021.output_cdsSequence.txt\n"
    echo "Exit Code: 1"
    return 1
else
    printf "Test passed.\n"
    rm cdsSequence 021.output_cdsSequence.txt 
    echo "Exit Code: 0"
    return 0
fi

