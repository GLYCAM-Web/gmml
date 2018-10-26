#!/bin/bash

# Required for compiling
# Can do this, or with root copy bin/* to /usr/lib/
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${GEMSHOME}/gmml/lib
export LD_LIBRARY_PATH


###################### Test 07 ######################
printf "Testing buildBySequence... "
g++ -std=c++11 -I$GEMSHOME/gmml/includes -L$GEMSHOME/gmml/lib -Wall buildBySequence.cc -o buildBySequence -lgmml
./buildBySequence 
#rm buildBySequence.pdb buildBySequence > /dev/null 2>&1

