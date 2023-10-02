#!/bin/bash

GMML_ROOT_DIR=$(git rev-parse --show-toplevel)

if [[ "${GMML_ROOT_DIR}" != *"gmml" ]]; then
    echo -e "Test 002 failed, we think our GMML root directory is:\t${GMML_ROOT_DIR}\n"
    exit 1
fi

printf "Testing 002.create_Assembly_WritePDB.cc... "
g++ -std=c++17 -I "${GMML_ROOT_DIR}"/ -L"${GMML_ROOT_DIR}"/bin/ -Wl,-rpath,"${GMML_ROOT_DIR}"/bin/ tests/002.create_Assembly_WritePDB.cc -lgmml -pthread -o create_Assembly_WritePDB
./create_Assembly_WritePDB >/dev/null 2>&1
if [ -f test-NLN.pdb ]; then
    if ! cmp test-NLN.pdb tests/correct_outputs/test-NLN.pdb >/dev/null 2>&1; then
        printf "Test FAILED!\n"
        echo "Exit Code: 1"
        return 1
    else
        printf "Test passed\n"
        rm test-NLN.pdb create_Assembly_WritePDB >/dev/null 2>&1
        echo "Exit Code: 0"
        return 0
    fi
else
    printf "Test FAILED!\n"
    echo "Exit Code: 1"
    return 1
fi
