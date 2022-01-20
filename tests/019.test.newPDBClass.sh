#!/bin/bash
printf "Testing 019.test.newPDBClass.cpp... "
g++ -std=c++17 -I $GEMSHOME/gmml/ -L$GEMSHOME/gmml/bin/ -Wl,-rpath,$GEMSHOME/gmml/bin/ tests/019.test.newPDBClass.cpp -lgmml -pthread -o newPdbClass
./newPdbClass tests/inputs/4mbz.pdb

