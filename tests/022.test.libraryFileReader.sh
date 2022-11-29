#!/bin/bash

printf "Testing 022.test.libraryFileReader.cpp... "
g++ -std=c++17 -I $GEMSHOME/gmml/ -L$GEMSHOME/gmml/bin/ -Wl,-rpath,$GEMSHOME/gmml/bin/ tests/022.test.cdsLibFile.cpp -lgmml -pthread -o libFileReader
./libFileReader > 022.output.txt

fileList=("libAsPdbFile.pdb" "libAsOffFile.off" "libAsLibFile.lib")
for file in ${fileList[@]}; 
do
  	if [ -f $file ]; then
  	    if ! cmp $file tests/correct_outputs/020.$file  > /dev/null 2>&1; then
  	        printf "Test FAILED!\n $file is different from tests/correct_outputs/022.$file\n"
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
rm libFileReader
echo "Exit Code: 0"
return 0
