#!/bin/bash

printf "Testing 020.test.parameterFiles.cpp... "
g++ -std=c++17 -I $GEMSHOME/gmml/ -L$GEMSHOME/gmml/bin/ -Wl,-rpath,$GEMSHOME/gmml/bin/ tests/020.test.parameterFiles.cpp -lgmml -pthread -o loadParameters
./loadParameters

fileList=("prepAsPrepFile.prep" "prepAsPdbFile.pdb" "prepAsOffFile.off" "prepAsLibFile.lib")
for file in ${fileList[@]}; 
do
  	if [ ! -f $file ]; then
  	    printf "Test FAILED!\n $file does not exist\n"
        echo "Exit Code: 1"
        return 1
    fi
  	if ! cmp -s $file tests/correct_outputs/020.$file;  then
  	    printf "Test FAILED!\n $file is different from tests/correct_outputs/020.$file\n"
        echo "Exit Code: 1"
        return 1
    fi
    rm $file          
done
printf "Test passed.\n"
rm loadParameters
echo "Exit Code: 0"
return 0
