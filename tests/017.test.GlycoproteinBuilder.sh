#!/bin/bash

printf "Testing 017.test.GlycoproteinBuilder.cpp... "
g++ -std=c++17 -I $GEMSHOME/gmml/ -L$GEMSHOME/gmml/bin/ -Wl,-rpath,$GEMSHOME/gmml/bin/ tests/017.test.GlycoproteinBuilder.cpp -lgmml -pthread -o gpBuilder
./gpBuilder 017.GlycoproteinBuilderInput.txt tests/inputs/ > 017.output_GlycoproteinBuilder.txt 2>&1
if ! cmp  tests/inputs/GlycoProtein_All_Resolved.pdb tests/correct_outputs/017.Glycoprotein_All_Resolved.pdb > /dev/null 2>&1; then
    printf "Test FAILED! tests/inputs/GlycoProtein_All_Resolved.pdb different from tests/correct_outputs/017.Glycoprotein_All_Resolved.pdb\n Compare using diff or VMD\n"
    echo "Exit Code: 1"
    return 1
else
    printf "Test passed.\n"
    rm gpBuilder 017.output_GlycoproteinBuilder.txt tests/inputs/GlycoProtein_All_Resolved.pdb
    echo "Exit Code: 0"
    return 0
fi
