#!/bin/bash

GMML_ROOT_DIR=$(git rev-parse --show-toplevel)

if [[ "${GMML_ROOT_DIR}" != *"gmml" ]]; then
    echo -e "Test 018 failed, we think our GMML root directory is:\t${GMML_ROOT_DIR}\n"
    exit 1
fi

printf "Testing 018.test.createGlycosylationTables.cpp... "
g++ -std=c++17 -I "${GMML_ROOT_DIR}"/ -L"${GMML_ROOT_DIR}"/bin/ -Wl,-rpath,"${GMML_ROOT_DIR}"/bin/ ../internalPrograms/createGlycosylationTables.cpp -lgmml -pthread -o gpBuilderTable
./gpBuilderTable tests/inputs/018.4mbzEdit.pdb >018.output_GlycoproteinBuilderTable.txt 2>&1
if ! cmp 018.output_GlycoproteinBuilderTable.txt tests/correct_outputs/018.output_GlycoproteinBuilderTable.txt >/dev/null 2>&1; then
    printf "Test FAILED!. tests/correct_outputs/018.output_GlycoproteinBuilderTable.txt different from 018.output_GlycoproteinBuilderTable.txt\ndiff tests/correct_outputs/018.output_GlycoproteinBuilderTable.txt 018.output_GlycoproteinBuilderTable.txt\n"
    diff tests/correct_outputs/018.output_GlycoproteinBuilderTable.txt 018.output_GlycoproteinBuilderTable.txt
    echo "Exit Code: 1"
    return 1
else
    printf "Test passed.\n"
    rm gpBuilderTable 018.output_GlycoproteinBuilderTable.txt
    echo "Exit Code: 0"
    return 0
fi
