#!/bin/bash

printf "Testing 021.test.cdsSequence.cpp... "
g++ -std=c++17 -I $GEMSHOME/gmml/ -I $GEMSHOME/gmml/ -L$GEMSHOME/gmml/bin/ -Wl,-rpath,$GEMSHOME/gmml/bin/ tests/021.test.cdsSequence.cpp -lgmml -pthread -o cdsSequence
./cdsSequence > 021.output_cdsSequence.txt 2>&1
fileList=("021.output_cdsSequence.txt" "021.sequenceAsOffFile.off" "021.sequenceAsPdbFile.pdb")
for file in ${fileList[@]}; 
do
  	if [ -f $file ]; then
  	    if ! cmp $file tests/correct_outputs/$file  > /dev/null 2>&1; then
  	        printf "Test FAILED!\n $file is different from tests/correct_outputs/$file\n"
            echo "Exit Code: 1"
            return 1
        else
            rm $file    
        fi
    else
        printf "Test FAILED!\n $file does not exist\n"
        echo "Exit Code: 1"
        return 1
    fi      
done
printf "Test passed.\n"
rm cdsSequence
echo "Exit Code: 0"
return 0
