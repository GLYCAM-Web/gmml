#!/bin/bash

GMML_ROOT_DIR=$(git rev-parse --show-toplevel)

if [[ "${GMML_ROOT_DIR}" != *"gmml" ]]; then
    echo -e "Test 017 failed, we think our GMML root directory is:\t${GMML_ROOT_DIR}\n"
    exit 1
fi

printf "Testing 017.test.GlycoproteinBuilder.cpp... "
g++ -std=c++17 -I "${GMML_ROOT_DIR}" -L"${GMML_ROOT_DIR}"/bin/ -Wl,-rpath,"${GMML_ROOT_DIR}"/bin/ "${GMML_ROOT_DIR}"/internalPrograms/GlycoproteinBuilder/main.cpp -lgmml -pthread -o gpBuilder
./gpBuilder tests/inputs/017.GlycoproteinBuilderInput.txt >output_GlycoproteinBuilder.txt 2>&1
fileList=("glycoprotein_initial.pdb" "glycoprotein.pdb" "glycoprotein.off" " glycoprotein_serialized.pdb" "output_GlycoproteinBuilder.txt")
for file in "${fileList[@]}"; do
    if [ -f "${file}" ]; then
        if ! cmp "${file}" tests/correct_outputs/017."${file}" >/dev/null 2>&1; then
            echo -e "Test FAILED!\n ${file} is different from tests/correct_outputs/017.${file}\n"
            echo "Exit Code: 1"
            return 1
            exit 1
        else
            rm "${file}"
        fi
    else
        echo -e "Test FAILED!\n ${file} does not exist\n"
        echo "Exit Code: 1"
        return 1
        exit 1
    fi
done
#These need only exist, as will have random data each time.
fileList=("0_glycoprotein.pdb" "1_glycoprotein.pdb")
for file in "${fileList[@]}"; do
    if [ ! -f "${file}" ]; then
        echo -e "Test FAILED!\n ${file} does not exist\n"
        echo "Exit Code: 1"
        return 1
        exit 1
    else
        rm "${file}"
    fi
done
printf "Test passed.\n"
rm gpBuilder
echo "Exit Code: 0"
return 0
exit 0
