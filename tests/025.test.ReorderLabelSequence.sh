#!/bin/bash

GMML_ROOT_DIR=$(git rev-parse --show-toplevel)

if [[ "${GMML_ROOT_DIR}" != *"gmml" ]] ; then
            echo "Test 009 failed, we think our GMML root directory is:\t${GMML_ROOT_DIR}\n"
            exit 1
fi

SourceCode="tests/025.cdsSequence.cpp"
Executable="cdsReorderSequence"
OutputFile="025.output.reorderedLabeledSequences.txt"
CompareFile="tests/correct_outputs/025.output.reorderedLabeledSequences.txt"
SequenceOne="DManpa1-3[DGalpb1-4DGalpb1-4DGalpb1-4DGalpb1-4]LRhapa1-OH"
SequenceTwo="DManp[2S,3Me]a1-6DManpa1-6[DGlcpNAcb1-4][DNeu5Aca2-6DGalpb1-4DGlcpNAc[3S]b1-2DManpa1-3]DManpb1-4DGlcpNAc[6Me]b1-4DGlcpNAcb1-OH";
SequenceThree="LFrupa2-6-DGlcpAa1-OME"

printf "Testing $SourceCode (Sequence reordering and labeling)... "
g++ -std=c++17 -I "${GMML_ROOT_DIR}"/ -L"${GMML_ROOT_DIR}"/bin/ -Wl,-rpath,"${GMML_ROOT_DIR}"/bin/ ${SourceCode} -lgmml -pthread -o ${Executable}

./${Executable} ${SequenceOne}  > "${OutputFile}"  2>&1
printf "\n***********************************************************************\n\n" >> "${OutputFile}"
./${Executable} ${SequenceTwo}  >> "${OutputFile}"  2>&1
printf "\n***********************************************************************\n\n" >> "${OutputFile}"
./${Executable} ${SequenceThree}  >> "${OutputFile}"  2>&1

if [ ! -f ${CompareFile} ] ; then
	printf "Test FAILED! (cannot find correct output to compare against)\n"
	echo "Exit Code: 1"
	return 1
fi
if ! cmp ${OutputFile} ${CompareFile} > /dev/null 2>&1; then
	printf "\nTest FAILED! Output file differs from reference.\nTry:\ndiff ${OutputFile} ${CompareFile}\n"
    echo "Exit Code: 1"
    return 1
fi
printf "Test passed.\n"
rm ${Executable} ${OutputFile}
echo "Exit Code: 0"
return 0

