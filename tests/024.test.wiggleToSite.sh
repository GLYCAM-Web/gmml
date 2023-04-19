#!/bin/bash
printf "Testing 024.wiggleToSite..."
g++ -std=c++17 -I $GEMSHOME/gmml/ -L$GEMSHOME/gmml/bin/ -Wl,-rpath,$GEMSHOME/gmml/bin/ ../internalPrograms/WiggleToSite/wiggleToSiteDriver.cpp -lgmml -pthread -o wiggleToSite

./wiggleToSite tests/inputs/024.input.txt > 024.output_wiggleToSite.txt 2>&1

if ! cmp 024.output_wiggleToSite.txt tests/correct_outputs/024.output_wiggleToSite.txt > /dev/null 2>&1; then
    printf "Test FAILED! Output file %s different from %s \n" 024.output_wiggleToSite.txt tests/correct_outputs/024.output_wiggleToSite.txt
    echo "Exit Code: 1"
    return 1
    exit 1
fi
if ! [ -f finished.pdb ]; then
  	echo "Test FAILED! Did not create finished.pdb"
   	echo "Exit Code: 1"
   	return 1
   	exit 1
fi
#if ! cmp finished.pdb tests/correct_outputs/024.finished.pdb > /dev/null 2>&1; then
#  	echo "Test FAILED! Created pdb finished.pdb is different from tests/correct_outputs/024.finished.pdb"
#   	echo "Exit Code: 1"
#   	return 1
#   	exit 1
#fi
printf "Test passed.\n"
rm -r 024.output_wiggleToSite.txt wiggleToSite initial.pdb initial.off superimposed.off superimposed.pdb finished.pdb finished.off
echo "Exit Code: 0"
return 0