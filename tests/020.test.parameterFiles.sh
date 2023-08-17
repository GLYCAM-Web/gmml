#!/bin/bash

printf "Testing 020.test.parameterFiles.cpp... "

GMML_ROOT_DIR=$(git rev-parse --show-toplevel)

if [[ "${GMML_ROOT_DIR}" != *"gmml" ]]; then
    echo -e "Test 020 failed, we think our GMML root directory is:\t${GMML_ROOT_DIR}\n"
    exit 1
fi

g++ -std=c++17 -I "${GMML_ROOT_DIR}"/ -L"${GMML_ROOT_DIR}"/bin/ -Wl,-rpath,"${GMML_ROOT_DIR}"/bin/ tests/020.test.parameterFiles.cpp -lgmml -pthread -o loadParameters
./loadParameters

fileList=("prepAsPrepFile.prep" "prepAsPdbFile.pdb" "prepAsOffFile.off" "prepAsLibFile.lib")
for file in "${fileList[@]}"; do
    if [ ! -f "${file}" ]; then
        echo -e "Test FAILED!\n ${file} does not exist\n"
        echo "Exit Code: 1"
        return 1
    fi
    if ! cmp -s "${file}" tests/correct_outputs/020."${file}"; then
        echo -e "Test FAILED!\n ${file} is different from tests/correct_outputs/020.${file}\n"
        echo "Exit Code: 1"
        return 1
    fi
    rm "${file}"
done
echo -e "Test passed.\n"
rm loadParameters
echo "Exit Code: 0"

return 0
