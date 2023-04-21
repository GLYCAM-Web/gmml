#!/bin/bash

GMML_ROOT_DIR=$(git rev-parse --show-toplevel)

if [ "$(git config --get remote.origin.url)" != "https://github.com/GLYCAM-Web/gmml" ]; then
            exit 1
fi
printf "Testing PDBPreprocessor... ~15 seconds\n"
g++ -std=c++17 -I "${GMML_ROOT_DIR}" -L"${GMML_ROOT_DIR}"/bin -Wl,-rpath,"${GMML_ROOT_DIR}"/bin tests/004.PDBPreprocessor.cc -lgmml -pthread -o PDBPreprocessor
 
for filePath in $(ls tests/inputs/004.preprocessorInput_*.pdb)
do
    ./PDBPreprocessor $filePath > /dev/null 2>&1
    if [ -f Processed.pdb ]; then
        filename=$(basename $filePath)
        if ! cmp Processed.pdb tests/correct_outputs/$filename > /dev/null 2>&1; then
            printf "Test FAILED!.\nINFO: These files are different: Processed.pdb tests/correct_outputs/$filename\n"
        	echo "Exit Code: 1"
            return 1
        else
            printf "Test passed for $filePath.\n"
        fi
    else
        printf "Test FAILED!.\nProcessed.pdb was not created for $filePath.\n"
        echo "Exit Code: 1"
        return 1
	fi
done
rm Processed.pdb PDBPreprocessor > /dev/null 2>&1
echo "Exit Code: 0"	
return 0
