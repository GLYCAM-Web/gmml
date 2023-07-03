#!/bin/bash

GMML_ROOT_DIR=$(git rev-parse --show-toplevel)

if [[ "${GMML_ROOT_DIR}" != *"gmml" ]] ; then
            echo "Test 009 failed, we think our GMML root directory is:\t${GMML_ROOT_DIR}\n"
            exit 1
fi

SourceCode="tests/009.reorderSequence.cc"
Executable="reorderSequence"
OutputFile="Reordered_Labeled_Sequences.txt"
CompareFile="tests/correct_outputs/Reordered_Labeled_Sequences.txt"
SequenceOne="DManpa1-3[DGalpb1-4DGalpb1-4DGalpb1-4DGalpb1-4]LRhapa1-OH"
SequenceTwo="DManp[2S,3Me]a1-6DManpa1-6[DGlcpNAcb1-4][DNeu5Aca2-6DGalpb1-4DGlcpNAc[3S]b1-2DManpa1-3]DManpb1-4DGlcpNAc[6Me]b1-4DGlcpNAcb1-OH";
SequenceThree="LFrupa2-6-DGlcpAa1-OME"

printf "Testing 009.reorderSequence.cc (Sequence reordering and labeling)... "
g++ -std=c++17 -I "${GMML_ROOT_DIR}"/ -L"${GMML_ROOT_DIR}"/bin/ -Wl,-rpath,"${GMML_ROOT_DIR}"/bin/ ${SourceCode} -lgmml -pthread -o ${Executable}

./reorderSequence ${SequenceOne}  > "${OutputFile}"  2>&1
./reorderSequence ${SequenceTwo}  >> "${OutputFile}"  2>&1
./reorderSequence ${SequenceThree}  >> "${OutputFile}"  2>&1

if [ ! -f ${CompareFile} ] ; then
	printf "Test FAILED!  (cannot find comparison results)\n"
	echo "Exit Code: 1"
	return 1
fi

if [ -f ${OutputFile} ]; then
    if ! cmp ${OutputFile} ${CompareFile} > /dev/null 2>&1; then
	printf "Test FAILED! (comparison of results)\n"
        echo "Exit Code: 1"
        return 1
    else
        printf "Test passed.\n"
        rm ${Executable} ${OutputFile}
        echo "Exit Code: 0"
        return 0
    fi
else
    printf "Test FAILED!  (no output file)\n"
    echo "Exit Code: 1"
    return 1
fi
