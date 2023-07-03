#!/bin/bash

GMML_ROOT_DIR=$(git rev-parse --show-toplevel)

if [[ "${GMML_ROOT_DIR}" != *"gmml" ]] ; then
            echo "Test 007 failed, we think our GMML root directory is:\t${GMML_ROOT_DIR}\n"
            exit 1
fi

###################### Test 08 ######################
printf "Testing 007.detectSugars.cc... "
g++ -std=c++17 -I "${GMML_ROOT_DIR}" -L"${GMML_ROOT_DIR}"/bin/ -Wl,-rpath,"${GMML_ROOT_DIR}"/bin/ tests/007.detectSugars.cc -lgmml -lstdc++fs -pthread -o detect_sugars
if [ -f gmmo.ttl ]; then
   mv gmmo.ttl gmmoBeforeTests.ttl > /dev/null 2>&1
fi
./detect_sugars tests/inputs/4mbz.pdb > /dev/null 2>&1
if [ -f gmmo.ttl ]; then
   if ! cmp gmmo.ttl tests/correct_outputs/gmmo.ttl > /dev/null 2>&1; then
       printf "\nTest FAILED! gmmo.ttl differs from tests/correct_outputs/gmmo.ttl\n"
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

