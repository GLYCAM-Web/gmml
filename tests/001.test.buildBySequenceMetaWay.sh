#!/bin/bash

printf "Testing buildBySequenceMeta... "
g++ -std=c++17 -I $GEMSHOME/gmml/ -L$GEMSHOME/gmml/bin/ -Wl,-rpath,$GEMSHOME/gmml/bin/ tests/001.buildBySequenceMeta.cc -lgmml -pthread -o buildBySequenceMeta
./buildBySequenceMeta > /dev/null 2>&1      
if [ -f structure.pdb ]; then
    if ! cmp structure.pdb tests/correct_outputs/001.buildBySequenceMeta.pdb > /dev/null 2>&1; then
        printf "Test FAILED!\nstructure.pdb is different from tests/correct_outputs/001.buildBySequenceMeta.pdb\n"
        echo "Exit Code: 1"
        return 1
    else
        printf "Test passed.\n"
        rm structure.pdb buildBySequenceMeta
        echo "Exit Code: 0"
        return 0
    fi
else
    printf "Test FAILED!\n"
    echo "Exit Code: 1"
    return 1
fi

