#!/bin/bash

# Required for compiling
# Can do this, or with root copy bin/* to /usr/lib/
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:../bin
export LD_LIBRARY_PATH

###################### Test 01 ###################### 
printf "Testing create_Assembly_WritePDB... "
g++ -I../includes/* -L../bin/ tests/create_Assembly_WritePDB.cc -lgmml -o create_Assembly_WritePDB
./create_Assembly_WritePDB > /dev/null 2>&1
if [ -f test-NLN.pdb ]; then
    if [ `diff test-NLN.pdb tests/correct_outputs/test-NLN.pdb` ]; then
        printf "Test failed.\n"
    else 
        printf "Test passed.\n"
    fi
else
    printf "Test failed.\n"
fi
rm test-NLN.pdb create_Assembly_WritePDB > /dev/null 2>&1


###################### Test 02 ######################
printf "Testing superimposition_Eigen... "
g++ -I../includes/* -L../bin/ tests/superimposition_Eigen.cc -lgmml -o superimposition_Eigen
./superimposition_Eigen > /dev/null 2>&1
if [ -f moved.pdb ]; then
    if [ `diff moved.pdb tests/correct_outputs/moved.pdb` ]; then
        printf "Test failed.\n"
    else
        printf "Test passed.\n"
    fi
else
    printf "Test failed.\n"
fi
rm moved.pdb superimposition_Eigen > /dev/null 2>&1

###################### Test 03 ######################
printf "Testing PDBpreprocessor... "
g++ -I../includes/* -L../bin/ tests/PDB_preprocessor.cc -lgmml -o PDB_preprocessor
./PDB_preprocessor > /dev/null 2>&1
if [ -f Processed.pdb ]; then
    if [ `diff Processed.pdb tests/correct_outputs/Processed.pdb` ]; then
        printf "Test failed.\n"
    else
        printf "Test passed.\n"
    fi
else
    printf "Test failed.\n"
fi
rm Processed.pdb PDB_preprocessor > /dev/null 2>&1
