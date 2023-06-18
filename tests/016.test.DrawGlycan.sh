#!/bin/bash

GMML_ROOT_DIR=$(git rev-parse --show-toplevel)

if [[ "${GMML_ROOT_DIR}" != *"gmml" ]] ; then
            echo "Test 016 failed, we think our GMML root directory is:\t${GMML_ROOT_DIR}\n"
            exit 1
fi

printf "Testing 016.test.DrawGlycan.cc..."
g++ -std=c++17 -I "${GMML_ROOT_DIR}"/ -I "${GMML_ROOT_DIR}"/ -L"${GMML_ROOT_DIR}"/bin/ -Wl,-rpath,"${GMML_ROOT_DIR}"/bin/ tests/016.test.DrawGlycan.cc -lgmml -pthread -o drawGlycan
./drawGlycan
> 016.output_drawGlycan.txt
for dotFile in `ls *.dot`
do
	cat $dotFile >> 016.output_drawGlycan.txt
	dotFileName="${dotFile%.*}"
    dot -Tsvg:cairo:cairo $dotFile -o $dotFileName.svg > /dev/null 2>&1
    rm $dotFile
done

for svgFile in `ls *.svg`
do
	cmp $svgFile tests/correct_outputs/016.output_SVGs/$svgFile
	if ! cmp $svgFile tests/correct_outputs/016.output_SVGs/$svgFile > /dev/null 2>&1; then
		printf "Test FAILED! Output file %s different to tests/correct_outputs/016.output_SVGs/%s\n" $svgFile $svgFile
	  	echo "Exit Code: 1"
	  	return 1
	fi
	rm $svgFile
done

if ! cmp  016.output_drawGlycan.txt tests/correct_outputs/016.output_drawGlycan.txt > /dev/null 2>&1; then
    printf "Test FAILED! Output file different\n"
    echo "Exit Code: 1"
    return 1
else
    printf "Test passed.\n"
    rm drawGlycan 016.output_drawGlycan.txt
    echo "Exit Code: 0"
    return 0
fi
