#!/bin/bash

printf "Testing 026.test.editPDB.cpp... ~2 seconds. "

GMML_ROOT_DIR=$(git rev-parse --show-toplevel)

if [[ "${GMML_ROOT_DIR}" != *"gmml" ]]; then
    echo -e "Test 026 failed, we think our GMML root directory is:\t${GMML_ROOT_DIR}\n"
    exit 1
fi

g++ -std=c++17 -I "${GMML_ROOT_DIR}" -L"${GMML_ROOT_DIR}"/bin/ -Wl,-rpath,"${GMML_ROOT_DIR}"/bin/ tests/026.test.editPDB.cpp -lgmml -pthread -o editPdb
shopt -s nullglob
for filepath in tests/inputs/026.*.pdb; do
    file=$(basename "${filepath}")
    ./editPdb tests/inputs/"${file}" >026.output.txt
    if ! cmp 026.outputPdbFile.pdb tests/correct_outputs/"${file}"-output.pdb >/dev/null 2>&1; then
        echo -e "Test FAILED!. 026.outputPdbFile.pdb different from tests/correct_outputs/${file}-output.pdb\n Compare using diff or VMD\n"
        echo "Exit Code: 1"
        return 1
        #exit 1
    fi
done
printf "Test passed.\n"
rm ./editPdb 026.outputPdbFile.pdb 026.output.txt
echo "Exit Code: 0"
return 0
