#!/bin/bash

required_passing_tests=$(/bin/ls -1 *.test.*.sh | wc -l)
echo """
Number of tests found: ${required_passing_tests}
Beginning testing.
"""

run_test() 
{
    sh $1
    return $?
}

##
##   To disable a test:
##
##        Change the filename so that it 
##        doesn't match the pattern *.test.*.sh
##

tests_attempted=0
tests_passed=0
for i in $(/bin/ls *.test.*.sh) ; do 
	printf "Using test file:  ${i} \n"
	tests_attempted=$((tests_attempted+1))
	if run_test ${i} ; then tests_passed=$((tests_passed+1)); fi
done
echo """
$tests_attempted tests were attempted
$tests_passed tests passed 
$required_passing_tests were required"

if [ "$tests_passed" -ge "$required_passing_tests" ]; then 
    #echo "The required number of tests passed"
    exit 0
else
    #echo "Some required tests did not pass."
    exit 1
fi
