#!/bin/bash
## Note: Oliver was checking the functionality. It does not yet work as required, but it took a while to figure out how to run the code
## So this test is just a snapshot of how it's currently working and how I managed to get output.
printf "Testing 012.AddSolventNeutralize... "
g++ -std=c++0x -I $GEMSHOME/gmml/ -L$GEMSHOME/gmml/bin/ -Wl,-rpath,$GEMSHOME/gmml/bin/ tests/012.addSolventNeutralize.cpp -lgmml -pthread -o addSolventNeutralize
./addSolventNeutralize > 012.output_addSolventNeutralize.txt
if [ -f 012.addSolventNeutralize.pdb ] ; then
    if ! cmp 012.addSolventNeutralize.pdb tests/correct_outputs/012.addSolventNeutralize.pdb > /dev/null 2>&1; then
        printf "Test FAILED!. PDB file different\n"
        return 1;
    #elif ! cmp structure.off tests/correct_outputs/010.buildBySequenceRotamer.off > /dev/null 2>&1; then
    #    printf "Test FAILED!. Off file different.\n"
    #    return 1;
    elif ! cmp  012.output_addSolventNeutralize.txt tests/correct_outputs/012.output_addSolventNeutralize.txt > /dev/null 2>&1; then
        printf "Test FAILED!. Output file different\n"
        return 1;
    else
        printf "Test passed.\n"
        rm 012.addSolventNeutralize.pdb addSolventNeutralize 012.output_addSolventNeutralize.txt
        return 0;
    fi
else
    printf "Test FAILED!.\n"
    return 1;
fi

