#!/bin/bash

GMML_ROOT_DIR=$(git rev-parse --show-toplevel)

if [ "$(git config --get remote.origin.url)" != "https://github.com/GLYCAM-Web/gmml.git" ]; then
            exit 1
fi

printf "Testing 010.buildBySequenceRotamer.cc... "
g++ -std=c++17 -I "${GMML_ROOT_DIR}"/ -I "${GMML_ROOT_DIR}"/ -L"${GMML_ROOT_DIR}"/bin/ -Wl,-rpath,"${GMML_ROOT_DIR}"/bin/ tests/010.buildBySequenceRotamer.cc -lgmml -pthread -o buildBySequenceRotamer
./buildBySequenceRotamer > 010.output_buildBySequenceRotamer.txt
if [ -f structure.pdb ] && [ -f structure.off ]; then
    if ! cmp structure.pdb tests/correct_outputs/010.buildBySequenceRotamer.pdb > /dev/null 2>&1; then
        printf "Test FAILED! structure.pdb file different from tests/correct_outputs/010.buildBySequenceRotamer.pdb\n"
        echo "Exit Code: 1"
        return 1
    #elif ! cmp structure.off tests/correct_outputs/010.buildBySequenceRotamer.off > /dev/null 2>&1; then
    #    printf "Test FAILED!. Off file different.\n"
    #    return 1;
    elif ! cmp  010.output_buildBySequenceRotamer.txt tests/correct_outputs/010.output_buildBySequenceRotamer.txt > /dev/null 2>&1; then
        printf "Test FAILED! Output file different\n"
        echo "Exit Code: 1"
        return 1
    else
        printf "Test passed.\n"
        rm structure.pdb structure.off buildBySequenceRotamer 010.output_buildBySequenceRotamer.txt
        echo "Exit Code: 0"
        return 0
    fi
else
    printf "Test FAILED!\nstructure.pdb doesn't exist\n"
    echo "Exit Code: 1"
    return 1
fi

