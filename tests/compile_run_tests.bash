#!/bin/bash

#Manually change this number as you add tests:
number_of_tests=9
tests_passed=0

# Required for compiling
# Can do this, or with root copy bin/* to /usr/lib/
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:../bin
export LD_LIBRARY_PATH

printf "$number_of_tests tests will be run.\n"
###################### Test 01 ###################### 
printf "Testing create_Assembly_WritePDB... "
g++ -std=c++0x -I $GEMSHOME/gmml/includes/ -L$GEMSHOME/gmml/bin/ -Wl,-rpath,$GEMSHOME/gmml/bin/ tests/create_Assembly_WritePDB.cc -lgmml -o create_Assembly_WritePDB
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
g++ -std=c++0x -I $GEMSHOME/gmml/includes/ -L$GEMSHOME/gmml/bin/ -Wl,-rpath,$GEMSHOME/gmml/bin/ tests/superimposition_Eigen.cc -lgmml -o superimposition_Eigen
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
g++ -std=c++0x -I $GEMSHOME/gmml/includes/ -L$GEMSHOME/gmml/bin/ -Wl,-rpath,$GEMSHOME/gmml/bin/ tests/PDB_preprocessor.cc -lgmml -o PDB_preprocessor
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
g++ -std=c++0x -I $GEMSHOME/gmml/includes/ -L$GEMSHOME/gmml/bin/ -Wl,-rpath,$GEMSHOME/gmml/bin/ tests/overlaps.cc -lgmml -o overlaps
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

###################### Test 05 ######################
printf "Testing to make sure there are no using namespace declarations... "
namespace_count=$( grep -r --exclude-dir=External_Libraries "using namespace" ../includes/ ../src/ | wc -l )
if [[ $namespace_count -eq 0 ]]; then
	printf "Test passed.\n"
	((tests_passed++))
else
	printf "Test FAILED!\n"
fi

###################### Test 06 ######################
printf "Testing BFMP Ring Shape Calculation... "
g++ -std=c++0x -I $GEMSHOME/gmml/includes/ -L$GEMSHOME/gmml/bin/ -Wl,-rpath,$GEMSHOME/gmml/bin/ tests/ring_shape_detection.cc -lgmml -o ring_shape_detection
./ring_shape_detection > ring_shape_detection.txt
if [ -f ring_conformations.txt ]; then
    if ! cmp ring_conformations.txt tests/correct_outputs/ring_conformations.txt > /dev/null 2>&1; then
        printf "Test FAILED!.\nPlease compare ring_conformations.txt to tests/correct_outputs/ring_conformations.txt"
    else
        printf "Test passed.\n"
        ((tests_passed++))
	rm ring_shape_detection ring_shape_detection.txt ring_conformations.txt > /dev/null 2>&1
    fi
else
    printf "Test FAILED!.\nPlease compare ring_conformations.txt to tests/correct_outputs/ring_conformations.txt"
fi

###################### Test 07 ######################
printf "Testing buildBySequence... "
g++ -std=c++0x -I $GEMSHOME/gmml/includes/ -L$GEMSHOME/gmml/bin/ -Wl,-rpath,$GEMSHOME/gmml/bin/ tests/buildBySequence.cc -lgmml -o buildBySequence
./buildBySequence > /dev/null 2>&1
if [ -f buildBySequence.pdb ]; then
    if ! cmp buildBySequence.pdb tests/correct_outputs/buildBySequence.pdb > /dev/null 2>&1; then
        printf "Test FAILED!.\n"
    else
        printf "Test passed.\n"
        ((tests_passed++))
	rm buildBySequence.pdb buildBySequence > /dev/null 2>&1
    fi
else
    printf "Test FAILED!.\n"
fi
#rm buildBySequence.pdb buildBySequence > /dev/null 2>&1

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
   else
       printf "Test passed.\n"
       ((tests_passed++))
       rm gmmo.ttl ring_conformations.txt detect_sugars > /dev/null 2>&1
       if [ -f gmmoBeforeTests.ttl ]; then
          mv gmmoBeforeTests.ttl gmmo.ttl > /dev/null 2>&1
       fi
   fi
else
   printf "Test FAILED!.\n"
fi
###################### Test 09 ######################
printf "Testing pdb2glycam and molecule subgraph matching... "
g++ -std=c++0x -I $GEMSHOME/gmml/includes/ -L$GEMSHOME/gmml/bin/ -Wl,-rpath,$GEMSHOME/gmml/bin/ tests/pdb2glycam.cc -lgmml -o pdb2glycam
./pdb2glycam tests/inputs/pdb2glycam_4YG0.pdb
if [ -f pdb2glycam_output.pdb ]; then
    if ! cmp pdb2glycam_output.pdb tests/correct_outputs/pdb2glycam_4YG0_output.pdb > /dev/null 2>&1; then
	printf "Test FAILED!.\n"
    else
	printf "Test passed.\n"
	((tests_passed++))
    fi
else
    printf "Test FAILED!.\n"
fi
rm ring_conformations.txt pdb2glycam_output.pdb pdb2glycam > /dev/null 2>&1

printf "Completed $number_of_tests tests.\n"
#####################################################

# ############# Allow git push ########################
if [[ $tests_passed -eq $number_of_tests ]]; then
   rm log.log > /dev/null 2>&1
   printf "All tests passed.\n"
   exit 0
else
   printf "\nError:  Only $tests_passed tests passed!\n"
   printf   "        Examine file log.log to investigate the failures.\n\n"
   exit 1
fi

