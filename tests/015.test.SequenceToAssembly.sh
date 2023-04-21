#!/bin/bash

GMML_ROOT_DIR=$(git rev-parse --show-toplevel)

if [ "$(git config --get remote.origin.url)" != "https://github.com/GLYCAM-Web/gmml" ]; then
            exit 1
fi

printf "Testing 015.test.SequenceAssembly.cc... "
g++ -std=c++17 -I "${GMML_ROOT_DIR}" -I "${GMML_ROOT_DIR}" -L"${GMML_ROOT_DIR}"/bin/ -Wl,-rpath,"${GMML_ROOT_DIR}"/bin/ tests/015.test.SequenceAssembly.cc -lgmml -lstdc++fs -pthread -o sequenceAssembly
./sequenceAssembly > 015.output_sequenceAssembly.txt 2>&1
if ! cmp  015.output_sequenceAssembly.txt tests/correct_outputs/015.output_sequenceAssembly.txt > /dev/null 2>&1; then
    printf "Test FAILED! Output file different\n"
    echo "Exit Code: 1"
    return 1
else
    printf "Test passed.\n"
    rm sequenceAssembly 015.output_sequenceAssembly.txt 
    echo "Exit Code: 0"
    return 0
fi
