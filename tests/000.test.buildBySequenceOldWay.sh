#!/bin/bash

GMML_ROOT_DIR=$(git rev-parse --show-toplevel)

if [ "$(git config --get remote.origin.url)" != "https://github.com/GLYCAM-Web/gmml" ]; then
            exit 1
fi

printf "Testing 000.buildBySequence.cc... "
g++ -std=c++17 -I "${GMML_ROOT_DIR}"/ -L"${GMML_ROOT_DIR}"/bin/ -Wl,-rpath,"${GMML_ROOT_DIR}"/bin/ tests/000.buildBySequence.cc -lgmml -pthread -o buildBySequence
./buildBySequence > /dev/null 2>&1
#./buildBySequence
if [ -f buildBySequence.pdb ]; then
    if ! cmp buildBySequence.pdb tests/correct_outputs/buildBySequence.pdb > /dev/null 2>&1; then
        printf "Test FAILED!\nbuildBySequence.pdb is different from tests/correct_outputs/buildBySequence.pdb\n"
        echo "Exit Code: 1"
        return 1
    else
        printf "Test passed\n"
        rm buildBySequence buildBySequence.pdb
        echo "Exit Code: 0"
        return 0
    fi
else
    printf "Test FAILED!\n buildBySequence.pdb does not exist.\n"
    echo "Exit Code: 1"
    return 1
fi
