#!/bin/bash

###################### Test 08 ######################
printf "Testing 007.detectSugars.cc... "
g++ -std=c++17 -I $GEMSHOME/gmml/ -L$GEMSHOME/gmml/bin/ -Wl,-rpath,$GEMSHOME/gmml/bin/ tests/007.detectSugars.cc -lgmml -pthread -o detect_sugars
if [ -f gmmo.ttl ]; then
   mv gmmo.ttl gmmoBeforeTests.ttl > /dev/null 2>&1
fi
./detect_sugars tests/inputs/4mbz.pdb > /dev/null 2>&1
if [ -f gmmo.ttl ]; then
   if ! cmp gmmo.ttl tests/correct_outputs/gmmo.ttl > /dev/null 2>&1; then
       printf "\nTest FAILED! gmml.ttl differs from tests/correct_outputs/gmmo.ttl\n"
       echo "Exit Code: 1"
       return 1
   else
       printf "Test passed.\n"
       rm gmmo.ttl ring_conformations.txt detect_sugars > /dev/null 2>&1
       if [ -f gmmoBeforeTests.ttl ]; then
          mv gmmoBeforeTests.ttl gmmo.ttl > /dev/null 2>&1
       fi
       echo "Exit Code: 0"
       return 0
   fi
else
   printf "\nTest FAILED!\n"
   echo "Exit Code: 1"
   return 1
fi

