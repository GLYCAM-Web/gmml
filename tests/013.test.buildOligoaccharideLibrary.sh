#!/bin/bash

printf "Testing buildOligosaccharide library... "
g++ -std=c++17 -I $GEMSHOME/gmml/includes/ -L$GEMSHOME/gmml/bin/ -Wl,-rpath,$GEMSHOME/gmml/bin/ tests/013.buildOligosaccharideLibrary.cc -lgmml -o buildOligosaccharideLibrary
./buildOligosaccharideLibrary tests/inputs/013.smallLibrary.txt _ 013.outputs ../dat/prep/GLYCAM_06j-1.prep > 013.output_buildOligosaccharideLibrary.txt
if [ -f structure.pdb ] && [ -f structure.off ]; then
    if ! cmp  013.output_buildOligosaccharideLibrary.txt tests/correct_outputs/013.output_buildOligosaccharideLibrary.txt > /dev/null 2>&1; then
        printf "Test FAILED!. Output file different\n"
        return 1;
    else
        printf "Test passed.\n"
        rm -r 013.outputs/ buildOligosaccharideLibrary 013.output_buildOligosaccharideLibrary.txt
        return 0;
    fi
else
    printf "Test FAILED!.\n"
    return 1;
fi

