#!/bin/bash

SourceCode="tests/reorderSequence.cc"
Executable="reorderSequence"
OutputFile="Reordered_Labeled_Sequences.txt"
CompareFile="tests/correct_outputs/Reordered_Labeled_Sequences.txt"
SequenceOne="DManpa1-3[DGalpb1-4DGalpb1-4DGalpb1-4DGalpb1-4]LRhapa1-OH"
SequenceTwo="DManp[2S,3Me]a1-6DManpa1-6[DGlcpNAcb1-4][DNeu5Aca2-6DGalpb1-4DGlcpNAc[3S]b1-2DManpa1-3]DManpb1-4DGlcpNAc[6Me]b1-4DGlcpNAcb1-OH";

printf "Testing Sequence reordering and labeling... "
g++ -std=c++0x -I $GEMSHOME/gmml/includes/ -L$GEMSHOME/gmml/bin/ -Wl,-rpath,$GEMSHOME/gmml/bin/ ${SourceCode}  -lgmml -o ${Executable}

./reorderSequence ${SequenceOne}  > "${OutputFile}"  2>&1
./reorderSequence ${SequenceTwo}  >> "${OutputFile}"  2>&1

if [ ! -f ${CompareFile} ] ; then
	printf "Test FAILED!  (cannot find comparison results)\n"
	exit 1
fi

if [ -f ${OutputFile} ]; then
    if ! cmp ${OutputFile} ${CompareFile} > /dev/null 2>&1; then
	printf "Test FAILED! (comparison of results).\n"
        exit 1;
    else
        printf "Test passed.\n"
        rm ${Executable} ${OutputFile}
        exit 0;
    fi
else
    printf "Test FAILED!  (no output file) .\n"
    exit 1;
fi
