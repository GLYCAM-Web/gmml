#!/bin/bash

printf "Testing 021.test.cdsSequence.cpp... "

GMML_ROOT_DIR=$(git rev-parse --show-toplevel)

if [[ "${GMML_ROOT_DIR}" != *"gmml" ]]; then
    echo -e "Test 021 failed, we think our GMML root directory is:\t${GMML_ROOT_DIR}\n"
    exit 1
fi

g++ -std=c++17 -I "${GMML_ROOT_DIR}"/ -I "${GMML_ROOT_DIR}"/ -L"${GMML_ROOT_DIR}"/bin/ -Wl,-rpath,"${GMML_ROOT_DIR}"/bin/ tests/021.test.cdsSequence.cpp -lgmml -pthread -o cdsSequence
./cdsSequence >021.output_cdsSequence.txt 2>&1
fileList=("021.output_cdsSequence.txt" "021.sequenceAsOffFile.off" "021.sequenceAsPdbFile.pdb")
for file in "${fileList[@]}"; do
    if [ -f "${file}" ]; then
        if ! cmp "${file}" tests/correct_outputs/"${file}" >/dev/null 2>&1; then
            echo -e "Test FAILED!\n ${file} is different from tests/correct_outputs/${file}\n"
            echo "Exit Code: 1"
            return 1
        else
            rm "${file}"
        fi
    else
        echo -e "Test FAILED!\n ${file} does not exist\n"
        echo "Exit Code: 1"
        return 1
    fi
done
printf "Test passed.\n"
rm cdsSequence
echo "Exit Code: 0"
return 0
