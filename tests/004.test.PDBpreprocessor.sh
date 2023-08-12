#!/bin/bash

GMML_ROOT_DIR=$(git rev-parse --show-toplevel)

if [[ "${GMML_ROOT_DIR}" != *"gmml" ]]; then
    echo -e "Test 004 failed, we think our GMML root directory is:\t${GMML_ROOT_DIR}\n"
    exit 1
fi

printf "Testing PDBPreprocessor... ~15 seconds\n"
g++ -std=c++17 -I "${GMML_ROOT_DIR}" -L"${GMML_ROOT_DIR}"/bin -Wl,-rpath,"${GMML_ROOT_DIR}"/bin tests/004.PDBPreprocessor.cc -lgmml -pthread -o PDBPreprocessor

for filePath in $(ls tests/inputs/004.preprocessorInput_*.pdb); do
    ./PDBPreprocessor "${filePath}" >/dev/null 2>&1
    if [ -f Processed.pdb ]; then
        filename=$(basename "${filePath}")
        if ! cmp Processed.pdb tests/correct_outputs/"${filename}" >/dev/null 2>&1; then
            echo -e "Test FAILED!.\nINFO: These files are different: Processed.pdb tests/correct_outputs/${filename}\n"
            echo "Exit Code: 1"
            return 1
        else
            echo -e "Test passed for ${filePath}.\n"
        fi
    else
        echo -e "Test FAILED!.\nProcessed.pdb was not created for ${filePath}.\n"
        echo "Exit Code: 1"
        return 1
    fi
done
rm Processed.pdb PDBPreprocessor >/dev/null 2>&1
echo "Exit Code: 0"
return 0
