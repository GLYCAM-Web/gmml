#!/bin/bash
printf "Testing 028.test.cdsCarbBuilderAll.cpp..."

GMML_ROOT_DIR=$(git rev-parse --show-toplevel)

if [[ "${GMML_ROOT_DIR}" != *"gmml" ]]; then
    echo -e "Test 028 failed, we think our GMML root directory is:\t${GMML_ROOT_DIR}\n"
    exit 1
fi

g++ -std=c++17 -I "${GMML_ROOT_DIR}"/ -L"${GMML_ROOT_DIR}"/bin/ -Wl,-rpath,"${GMML_ROOT_DIR}"/bin/ tests/028.test.cdsCarbBuilderAll.cpp -lgmml -pthread -o 028.carbBuilder.exe
rm -r 028.outputs/ >/dev/null 2>&1
mkdir 028.outputs/
./028.carbBuilder.exe tests/inputs/023.smallLibrary.txt _ 028.outputs ../dat/prep/GLYCAM_06j-1_GAGS_KDN.prep >028.output_carbohydrateBuilder.txt 2>&1

for i in $(cut -d _ -f1 tests/inputs/023.smallLibrary.txt); do
    if [ -f 028.outputs/"${i}".pdb ]; then
        echo "${i}.pdb succesfully created." >>028.output_carbohydrateBuilder.txt
        if ! cmp 028.outputs/"${i}".pdb tests/correct_outputs/023.outputs/"${i}".pdb >/dev/null 2>&1; then
            echo "Test FAILED! Created pdb file 028.outputs/${i}.pdb is different from tests/correct_outputs/023.outputs/${i}.pdb"
            echo "Exit Code: 1"
            return 1
        fi
    else
        echo "${i}.pdb not created." >>028.output_carbohydrateBuilder.txt
        if [ -f tests/correct_outputs/028.outputs/"${i}".pdb ]; then
            echo "Test FAILED! Did not create ${i}.pdb, yet it exists in tests/correct_outputs/028.outputs/${i}.pdb"
            echo "Exit Code: 1"
            return 1
        fi
    fi
done
if ! cmp 028.output_carbohydrateBuilder.txt tests/correct_outputs/028.output_carbohydrateBuilder.txt >/dev/null 2>&1; then
    printf "Test FAILED! Output file %s different from %s \n" 028.output_carbohydrateBuilder.txt tests/correct_outputs/028.output_carbohydrateBuilder.txt
    echo "Exit Code: 1"
    return 1
else
    printf "Test passed.\n"
    rm -r 028.outputs/  028.carbBuilder.exe 028.output_carbohydrateBuilder.txt
    echo "Exit Code: 0"
    return 0
fi
