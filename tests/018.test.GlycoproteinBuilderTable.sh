#!/bin/bash

printf "Testing 018.test.createGlycosylationTables.cpp... "
g++ -std=c++17 -I $GEMSHOME/gmml/ -L$GEMSHOME/gmml/bin/ -Wl,-rpath,$GEMSHOME/gmml/bin/ ../internalPrograms/createGlycosylationTables.cpp -lgmml -pthread -o gpBuilderTable
./gpBuilderTable tests/inputs/018.4mbzEdit.pdb > 018.output_GlycoproteinBuilderTable.txt 2>&1
if ! cmp  018.output_GlycoproteinBuilderTable.txt tests/correct_outputs/018.output_GlycoproteinBuilderTable.txt > /dev/null 2>&1; then
    printf "Test FAILED!. tests/correct_outputs/018.output_GlycoproteinBuilderTable.txt different from 018.output_GlycoproteinBuilderTable.txt\ndiff tests/correct_outputs/018.output_GlycoproteinBuilderTable.txt 018.output_GlycoproteinBuilderTable.txt\n"
    diff tests/correct_outputs/018.output_GlycoproteinBuilderTable.txt 018.output_GlycoproteinBuilderTable.txt
    echo "Exit Code: 1"
    return 1;
else
    printf "Test passed.\n"
    rm gpBuilderTable 018.output_GlycoproteinBuilderTable.txt
    echo "Exit Code: 0"	
    return 0;
fi
