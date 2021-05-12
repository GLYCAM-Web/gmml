#!/bin/bash

printf "Testing 014.test.SequenceParser.cc... "
g++ -std=c++17 -I $GEMSHOME/gmml/ -I $GEMSHOME/gmml/ -L$GEMSHOME/gmml/bin/ -Wl,-rpath,$GEMSHOME/gmml/bin/ tests/014.test.SequenceParser.cc -lgmml -pthread -o sequenceParser
./sequenceParser > 014.output_sequenceParser.txt
if ! cmp  014.output_sequenceParser.txt tests/correct_outputs/014.output_sequenceParser.txt > /dev/null 2>&1; then
    printf "Test FAILED!. Output file different\n"
    return 1;
else
    printf "Test passed.\n"
    rm sequenceParser 014.output_sequenceParser.txt
    return 0;
fi

