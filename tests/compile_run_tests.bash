#!/bin/bash

#Manually change this number as you add tests:
number_of_tests=4
tests_passed=0

# Required for compiling
# Can do this, or with root copy bin/* to /usr/lib/
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:../bin
export LD_LIBRARY_PATH

###################### Test 01 ###################### 
printf "Testing create_Assembly_WritePDB... "
g++ -I../includes/* -L../bin/ tests/create_Assembly_WritePDB.cc -lgmml -o create_Assembly_WritePDB
./create_Assembly_WritePDB > /dev/null 2>&1
if [ -f test-NLN.pdb ]; then
    if ! cmp test-NLN.pdb tests/correct_outputs/test-NLN.pdb > /dev/null 2>&1; then
        printf "Test FAILED!.\n"
    else 
        printf "Test passed.\n"
        ((tests_passed++))
    fi
else
    printf "Test FAILED!.\n"
fi
rm test-NLN.pdb create_Assembly_WritePDB > /dev/null 2>&1


###################### Test 02 ######################
printf "Testing superimposition_Eigen... "
g++ -I../includes/* -L../bin/ tests/superimposition_Eigen.cc -lgmml -o superimposition_Eigen
./superimposition_Eigen > /dev/null 2>&1
if [ -f moved.pdb ]; then
    if ! cmp moved.pdb tests/correct_outputs/moved.pdb > /dev/null 2>&1; then
        printf "Test FAILED!.\n"
    else
        printf "Test passed.\n"
        ((tests_passed++))
    fi
else
    printf "Test FAILED!.\n"
fi
rm moved.pdb superimposition_Eigen > /dev/null 2>&1

###################### Test 03 ######################
printf "Testing PDBpreprocessor... "
g++ -I../includes/* -L../bin/ tests/PDB_preprocessor.cc -lgmml -o PDB_preprocessor
./PDB_preprocessor > /dev/null 2>&1
if [ -f Processed.pdb ]; then
    if ! cmp Processed.pdb tests/correct_outputs/Processed.pdb > /dev/null 2>&1; then
        printf "Test FAILED!.\n"
    else
        printf "Test passed.\n"
        ((tests_passed++))
    fi
else
    printf "Test FAILED!.\n"
fi
rm Processed.pdb PDB_preprocessor > /dev/null 2>&1

###################### Test 04 ######################
printf "Testing Overlaps function... "
g++ -I../includes/* -L../bin/ tests/overlaps.cc -lgmml -o overlaps
./overlaps > overlaps.txt
if [ -f overlaps.txt ]; then
    if ! cmp overlaps.txt tests/correct_outputs/overlaps.txt > /dev/null 2>&1; then
        printf "Test FAILED!.\n"
    else
        printf "Test passed.\n"
        ((tests_passed++))
    fi
else
    printf "Test FAILED!.\n"
fi
rm overlaps overlaps.txt > /dev/null 2>&1


############# Allow git push ########################
if [[ $tests_passed -eq $number_of_tests ]]; then
   exit 0 #All tests passed
else
   exit 1
fi

