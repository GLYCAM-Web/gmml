#!/bin/bash

printf "Testing readJSON... "
g++ -std=c++0x -I $GEMSHOME/gmml/includes/ -L$GEMSHOME/gmml/bin/ -Wl,-rpath,$GEMSHOME/gmml/bin/ tests/readJSON.cc -lgmml -o readJSON
./readJSON
#./readJSON > /dev/null 2>&1
#if [ -f readJSON.pdb ]; then
#    if ! cmp readJSON.pdb tests/correct_outputs/readJSON.pdb > /dev/null 2>&1; then
#        printf "Test FAILED!.\n"
#        return 1;
#    else
#        printf "Test passed.\n"
#        ((tests_passed++))
#        return 0;
#    fi
#else
#    printf "Test FAILED!.\n"
#    return 1;
#fi

