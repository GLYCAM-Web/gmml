#!/bin/bash

printf "Testing 027.test.glycamResidueCombinator.cpp... "

GMML_ROOT_DIR=$(git rev-parse --show-toplevel)

if [[ "${GMML_ROOT_DIR}" != *"gmml" ]]; then
    echo -e "Test 027 failed, we think our GMML root directory is:\t${GMML_ROOT_DIR}\n"
    exit 1
fi

g++ -std=c++17 -I "${GMML_ROOT_DIR}"/ -L"${GMML_ROOT_DIR}"/bin/ -Wl,-rpath,"${GMML_ROOT_DIR}"/bin/ tests/027.test.glycamResidueCombinator.cpp -lgmml -pthread -o 027.glycamResiduecombinator.exe
outputFile=027.output.txt
./027.glycamResiduecombinator.exe ../dat/prep/GLYCAM_06j-1_GAGS_KDN.prep > $outputFile
fileList=("GLYCAM_06k.lib")
for file in "${fileList[@]}"; do
    if [ ! -f "${file}" ]; then
        echo -e "Test FAILED!\n ${file} does not exist\n"
        echo "Exit Code: 1"
        return 1
    fi
    if ! cmp -s "${file}" ../dat/CurrentParams/"${file}"; then
        echo -e "Test FAILED!\n ${file} is different from ../dat/CurrentParams/${file}\n"
        echo "Exit Code: 1"
        return 1
    fi
    rm "${file}"
done
if ! cmp -s $outputFile tests/correct_outputs/$outputFile; then
    echo -e "Test FAILED!\nOutput files different. Try:\ndiff $outputFile tests/correct_outputs/$outputFile"
    echo "Exit Code: 1"
    return 1    
fi
    
echo -e "Test passed.\n"
rm 027.glycamResiduecombinator.exe $outputFile 
echo "Exit Code: 0"

return 0
