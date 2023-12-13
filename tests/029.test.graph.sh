#!/bin/bash
printf "Testing 029.graph..."

GMML_ROOT_DIR=$(git rev-parse --show-toplevel)

if [[ "${GMML_ROOT_DIR}" != *"gmml" ]]; then
    echo -e "Test 029 failed, we think our GMML root directory is:\t${GMML_ROOT_DIR}\n"
    exit 1
fi

g++ -std=c++17 -I "${GMML_ROOT_DIR}"/ -L"${GMML_ROOT_DIR}"/bin/ -Wl,-rpath,"${GMML_ROOT_DIR}"/bin/ tests/029.test.graph.cpp -lgmml -pthread -o 029.graph.exe

./029.graph.exe >029.output_graph.txt 2>&1

if ! cmp 029.output_graph.txt tests/correct_outputs/029.output_graph.txt >/dev/null 2>&1; then
    printf "Test FAILED! Output file different. Try\ndiff %s %s\n" 029.output_graph.txt tests/correct_outputs/029.output_graph.txt
    echo "Exit Code: 1"
    return 1
    exit 1
fi
#if ! [ -f finished.pdb ]; then
#    echo "Test FAILED! Did not create finished.pdb"
#    echo "Exit Code: 1"
#    return 1
#    exit 1
#fi
printf "Test passed.\n"
rm -r 029.output_graph.txt ./029.graph.exe 
echo "Exit Code: 0"
return 0
