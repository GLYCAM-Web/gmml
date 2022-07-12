#!/bin/bash

#Manually change this number as you add tests:
number_of_tests=1
tests_passed=0

##################### Test 1 ########################
echo "Testing Glycoprotein Builder..."
./bin/gpBuilder input.txt tests/simple/ > test1_output
#DIFF=$(diff test1_output tests/simple/output.txt)
#if [ "$DIFF" != "" ]; then
if grep -q "Program got to end ok" test1_output; then 
    if cmp tests/simple/savedOutput_GlycoProtein_All_Resolved.pdb tests/simple/GlycoProtein_All_Resolved.pdb > /dev/null 2>&1; then
        echo "Test passed."
        ((tests_passed++))
        rm test1_output
	rm ASN_*_glycan.pdb
    else
        echo "Structure files different, ( tests/simple/savedOutput_GlycoProtein_All_Resolved.pdb Vs tests/simple/GlycoProtein_All_Resolved.pdb ) test FAILED!"
    fi
else 
   echo "Test FAILED, program did not get to the end ok"
fi

############# Allow git Pushes ###################
if [[ $tests_passed == $number_of_tests ]]; then
    exit 0
    echo "All tests passed"
else
    exit 1
fi
