#!/bin/bash
printf "Testing buildOligosaccharide library... "
g++ -std=c++17 -I $GEMSHOME/gmml/ -L$GEMSHOME/gmml/bin/ -Wl,-rpath,$GEMSHOME/gmml/bin/ tests/013.buildOligosaccharideLibrary.cpp -lgmml -pthread -o buildOligosaccharideLibrary
rm -r 013.outputs/ >/dev/null 2>&1
mkdir 013.outputs/
./buildOligosaccharideLibrary tests/inputs/013.smallLibrary.txt _ 013.outputs ../dat/prep/GLYCAM_06j-1.prep > 013.output_buildOligosaccharideLibrary.txt 2>&1

for i in `cut -d _ -f1 tests/inputs/013.smallLibrary.txt`;
do
    if [ -f 013.outputs/$i.pdb ]; then
        echo "$i.pdb succesfully created." >> 013.output_buildOligosaccharideLibrary.txt
        if ! cmp 013.outputs/$i.pdb tests/correct_outputs/013.outputs/$i.pdb > /dev/null 2>&1; then
        	echo "Test FAILED!. Created pdb file 013.outputs/$i.pdb is different from tests/correct_outputs/013.outputs/$i.pdb"
        	return 1;
        fi
    else
        echo "$i.pdb not created." >> 013.output_buildOligosaccharideLibrary.txt
        if [ -f tests/correct_outputs/013.outputs/$i.pdb ]; then
        	echo "Test FAILED!. Did not create $i.pdb, yet it exists in tests/correct_outputs/013.outputs/$i.pdb"
        	return 1;
        fi
    fi
done
if ! cmp  013.output_buildOligosaccharideLibrary.txt tests/correct_outputs/013.output_buildOligosaccharideLibrary.txt > /dev/null 2>&1; then
    printf "Test FAILED!. Output file different\n"
    return 1;
else
    printf "Test passed.\n"
    rm -r 013.outputs/ buildOligosaccharideLibrary 013.output_buildOligosaccharideLibrary.txt
    return 0;
fi


