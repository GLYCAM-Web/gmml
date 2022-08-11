#!/bin/bash

printf "Testing 020.test.parameterFiles.cpp... "
g++ -std=c++17 -I $GEMSHOME/gmml/ -L$GEMSHOME/gmml/bin/ -Wl,-rpath,$GEMSHOME/gmml/bin/ tests/020.test.parameterFiles.cpp -lgmml -pthread -o loadParameters
./loadParameters

