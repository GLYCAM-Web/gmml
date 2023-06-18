#!/bin/bash

GMML_ROOT_DIR=$(git rev-parse --show-toplevel)

if [[ "${GMML_ROOT_DIR}" != *"gmml" ]] ; then
            echo "Test 005 failed, we think our GMML root directory is:\t${GMML_ROOT_DIR}\n"
            exit 1
fi

printf "Testing 005.overlaps.cc... "
g++ -std=c++17 -I "${GMML_ROOT_DIR}"/ -L"${GMML_ROOT_DIR}"/bin/ -Wl,-rpath,"${GMML_ROOT_DIR}"/bin/ tests/005.overlaps.cc -lgmml -pthread -o overlaps
./overlaps > overlaps.txt 2> /dev/null
if [ -f overlaps.txt ]; then
    if ! cmp overlaps.txt tests/correct_outputs/overlaps.txt > /dev/null 2>&1; then
        printf "Test FAILED!\n"
        echo "Exit Code: 1"
        return 1
    else
        printf "Test passed\n"
        rm overlaps overlaps.txt > /dev/null 2>&1
        echo "Exit Code: 0"
        return 0
    fi
else
    printf "Test FAILED!\n"
    echo "Exit Code: 1"
    return 1
fi
