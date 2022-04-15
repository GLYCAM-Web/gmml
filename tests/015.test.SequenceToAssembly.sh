#!/bin/bash

printf "Testing 015.test.SequenceAssembly.cpp... "
g++ -std=c++17 -I $GEMSHOME/gmml/ -I $GEMSHOME/gmml/ -L$GEMSHOME/gmml/bin/ -Wl,-rpath,$GEMSHOME/gmml/bin/ tests/015.test.SequenceAssembly.cpp -lgmml -pthread -o sequenceAssembly
./sequenceAssembly > 015.output_sequenceAssembly.txt 2>&1
if ! cmp  015.output_sequenceAssembly.txt tests/correct_outputs/015.output_sequenceAssembly.txt > /dev/null 2>&1; then
    printf "Test FAILED!. Output file different\n"
    return 1;
else
    printf "Test passed.\n"
    rm sequenceAssembly 015.output_sequenceAssembly.txt 
    return 0;
fi

