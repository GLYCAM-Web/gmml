#!/bin/bash

GMML_ROOT_DIR=$(git rev-parse --show-toplevel)

if [[ "${GMML_ROOT_DIR}" != *"gmml" ]]; then
    echo -e "Test 006 failed, we think our GMML root directory is:\t${GMML_ROOT_DIR}\n"
    exit 1
fi

printf "Testing 006.ringShapeDetection.cc (BFMP Ring Shape Calculation)... "
g++ -std=c++17 -I "${GMML_ROOT_DIR}"/ -L"${GMML_ROOT_DIR}"/bin/ -Wl,-rpath,"${GMML_ROOT_DIR}"/bin/ tests/006.ringShapeDetection.cc -lgmml -pthread -o ring_shape_detection
./ring_shape_detection >ring_shape_detection.txt
if [ -f ring_conformations.txt ]; then
    if ! cmp ring_conformations.txt tests/correct_outputs/ring_conformations.txt >/dev/null 2>&1; then
        printf "Test FAILED!\nPlease compare ring_conformations.txt to tests/correct_outputs/ring_conformations.txt"
        echo "Exit Code: 1"
        return 1
    else
        printf "Test passed.\n"
        rm ring_shape_detection ring_shape_detection.txt ring_conformations.txt >/dev/null 2>&1
        echo "Exit Code: 0"
        return 0
    fi
else
    printf "Test FAILED!\nPlease compare ring_conformations.txt to tests/correct_outputs/ring_conformations.txt"
    echo "Exit Code: 1"
    return 1
fi
