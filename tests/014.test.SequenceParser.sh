#!/bin/bash

GMML_ROOT_DIR=$(git rev-parse --show-toplevel)

if [[ "${GMML_ROOT_DIR}" != *"gmml" ]] ; then
            echo "Test 014 failed, we think our GMML root directory is:\t${GMML_ROOT_DIR}\n"
            exit 1
fi

printf "Testing 014.test.SequenceParser.cc... "
g++ -std=c++17 -I "${GMML_ROOT_DIR}"/ -I "${GMML_ROOT_DIR}"/ -L"${GMML_ROOT_DIR}"/bin/ -Wl,-rpath,"${GMML_ROOT_DIR}"/bin/ tests/014.test.SequenceParser.cc -lgmml -pthread -o sequenceParser
./sequenceParser > 014.output_sequenceParser.txt 2>&1
if ! cmp  014.output_sequenceParser.txt tests/correct_outputs/014.output_sequenceParser.txt > /dev/null 2>&1; then
    printf "Test FAILED! Output file different, do:\ndiff 014.output_sequenceParser.txt tests/correct_outputs/014.output_sequenceParser.txt\n"
    echo "Exit Code: 1"
    return 1
else
    printf "Test passed.\n"
    rm sequenceParser 014.output_sequenceParser.txt
    echo "Exit Code: 0"
    return 0
fi

