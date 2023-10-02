#!/bin/bash

printf "Testing 019.test.newPDBClass.cpp... ~30 seconds. "

GMML_ROOT_DIR=$(git rev-parse --show-toplevel)

if [[ "${GMML_ROOT_DIR}" != *"gmml" ]]; then
    echo -e "Test 019 failed, we think our GMML root directory is:\t${GMML_ROOT_DIR}\n"
    exit 1
fi

g++ -std=c++17 -I "${GMML_ROOT_DIR}" -L"${GMML_ROOT_DIR}"/bin/ -Wl,-rpath,"${GMML_ROOT_DIR}"/bin/ tests/019.test.newPDBClass.cpp -lgmml -pthread -o newPdbClass
shopt -s nullglob
for filepath in tests/inputs/019.*.pdb; do
    file=$(basename "${filepath}")
    ./newPdbClass tests/inputs/"${file}" >019.output.txt
    if ! cmp 019.output.txt tests/correct_outputs/"${file}"-output.txt >/dev/null 2>&1; then
        echo -e "Test FAILED!. 019.output.txt different from tests/correct_outputs/${file}-output.txt\n Compare using diff\n"
        echo "Exit Code: 1"
        return 1
        #exit 1
    elif ! cmp 019.outputPdbFile.pdb tests/correct_outputs/"${file}"-output.pdb >/dev/null 2>&1; then
        echo -e "Test FAILED!. 019.outputPdbFile.pdb different from tests/correct_outputs/${file}-output.pdb\n Compare using diff or VMD\n"
        echo "Exit Code: 1"
        return 1
        #exit 1
    fi
done
printf "Test passed.\n"
rm ./newPdbClass 019.outputPdbFile.pdb 019.output.txt outputOffFile.off
echo "Exit Code: 0"
return 0
