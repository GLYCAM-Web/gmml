#!/bin/bash

#Manually change this number as you add tests:
required_passing_tests=9
tests_passed=0

run_test() 
{
    sh $1
    return $?
}

# Comment out tests you don't want to run. Add additional tests at the bottom. 
# Change required_passing_tests to equal number of tests.

#if run_test 001.test.detectSugar.sh; then tests_passed=$(($tests_passed + 1)); fi
if run_test 000.test.buildBySequenceOldWay.sh; then tests_passed=$(($tests_passed + 1)); fi
if run_test 001.test.buildBySequenceMetaWay.sh; then tests_passed=$(($tests_passed + 1)); fi
if run_test 002.test.createAssemblyWritePDB.sh; then tests_passed=$(($tests_passed + 1)); fi
if run_test 003.test.SuperimpositionEigen.sh; then tests_passed=$(($tests_passed + 1)); fi
if run_test 004.test.PDBpreprocessor.sh; then tests_passed=$(($tests_passed + 1)); fi
if run_test 005.test.Overlaps.sh; then tests_passed=$(($tests_passed + 1)); fi
if run_test 006.test.BFMP-RingShapeCalculation.sh; then tests_passed=$(($tests_passed + 1)); fi
if run_test 007.test.DetectSugars.sh; then tests_passed=$(($tests_passed + 1)); fi
if run_test 008.test.PDB2GlycamAndSubgraphMatching.sh; then tests_passed=$(($tests_passed + 1)); fi

echo "$tests_passed tests passed of $required_passing_tests"

if [[ "$tests_passed" -ge "$required_passing_tests" ]]; then
    exit 0
    echo "All tests passed"
else
    exit 1
fi
