#!/bin/bash

GMML_ROOT_DIR=$(git rev-parse --show-toplevel)

if [ "$(git config --get remote.origin.url)" != "https://github.com/GLYCAM-Web/gmml.git" ]; then
            exit 1
fi

printf "Testing 008.pdb2glycam.cc and molecule subgraph matching... "
g++ -std=c++17 -I "${GMML_ROOT_DIR}" -L"${GMML_ROOT_DIR}"/bin/ -Wl,-rpath,"${GMML_ROOT_DIR}"/bin/ "${GMML_ROOT_DIR}"/internalPrograms/convertPdbToGlycam.cpp -lgmml -pthread -o convertPdbToGlycam
./convertPdbToGlycam tests/inputs/008.convertPdbToGlycam_4YG0.pdb convertPdbToGlycam_output
if [ -f convertPdbToGlycam_output.pdb ]; then
    if ! cmp convertPdbToGlycam_output.pdb tests/correct_outputs/008.convertPdbToGlycam_4YG0_output.pdb > /dev/null 2>&1; then
        printf "\nTest FAILED!\n convertPdbToGlycam_output.pdb does not match tests/correct_outputs/008.convertPdbToGlycam_4YG0_output.pdb"
        echo "Exit Code: 1"
        return 1
        exit 1
    else
        printf "Test passed.\n"
        rm ring_conformations.txt convertPdbToGlycam_output.pdb convertPdbToGlycam > /dev/null 2>&1
        echo "Exit Code: 0"
        return 0
        exit 0
    fi
fi
printf "\nTest FAILED!\nconvertPdbToGlycam_output.pdb does not exist.\n"
echo "Exit Code: 1"
return 1
exit 1