#!/bin/bash

###################### Test 08 ######################
printf "Testing detectSugars... "
g++ -std=c++0x -I $GEMSHOME/gmml/includes/ -L$GEMSHOME/gmml/bin/ -Wl,-rpath,$GEMSHOME/gmml/bin/ tests/detect_sugars.cc -lgmml -o detect_sugars
if [ -f gmmo.ttl ]; then
   mv gmmo.ttl gmmoBeforeTests.ttl > /dev/null 2>&1
fi
./detect_sugars tests/inputs/4mbz.pdb > /dev/null 2>&1
if [ -f gmmo.ttl ]; then
   if ! cmp gmmo.ttl tests/correct_outputs/gmmo.ttl > /dev/null 2>&1; then
       printf "Test FAILED!.\n"
       return 1
   else
       printf "Test passed.\n"
       rm gmmo.ttl ring_conformations.txt detect_sugars > /dev/null 2>&1
       if [ -f gmmoBeforeTests.ttl ]; then
          mv gmmoBeforeTests.ttl gmmo.ttl > /dev/null 2>&1
       fi
       return 0
   fi
else
   printf "Test FAILED!.\n"
   return 1
fi

