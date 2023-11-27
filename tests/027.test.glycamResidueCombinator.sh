#!/bin/bash

printf "Testing 027.test.glycamResidueCombinator.cpp... "

GMML_ROOT_DIR=$(git rev-parse --show-toplevel)

if [[ "${GMML_ROOT_DIR}" != *"gmml" ]]; then
    echo -e "Test 027 failed, we think our GMML root directory is:\t${GMML_ROOT_DIR}\n"
    exit 1
fi

g++ -std=c++17 -I "${GMML_ROOT_DIR}"/ -L"${GMML_ROOT_DIR}"/bin/ -Wl,-rpath,"${GMML_ROOT_DIR}"/bin/ tests/027.test.glycamResidueCombinator.cpp -lgmml -pthread -o 027.glycamResiduecombinator.exe
./027.glycamResiduecombinator.exe ../dat/prep/GLYCAM_06j-1_GAGS_KDN.prep

fileList=("0GA.lib" "0YB.lib")
for file in "${fileList[@]}"; do
    if [ ! -f "${file}" ]; then
        echo -e "Test FAILED!\n ${file} does not exist\n"
        echo "Exit Code: 1"
        return 1
    fi
    if ! cmp -s "${file}" tests/correct_outputs/027."${file}"; then
        echo -e "Test FAILED!\n ${file} is different from tests/correct_outputs/027.${file}\n"
        echo "Exit Code: 1"
        return 1
    fi
    rm "${file}"
done
echo -e "Test passed.\n"
rm 027.glycamResiduecombinator.exe
echo "Exit Code: 0"

return 0
