#!/bin/bash

GMML_ROOT_DIR=$(git rev-parse --show-toplevel)

if [[ "${GMML_ROOT_DIR}" != *"gmml" ]]; then
    echo -e "Test 016 failed, we think our GMML root directory is:\t${GMML_ROOT_DIR}\n"
    exit 1
fi

#I am sorry for I have sinned, this is lazy and gross. I didnt want
#to deal with doing weird regex stuff within the diff file so I just
#back up the og test files then change the ones that will actually
#be changed so the pathing matches up to the bare metal paths
if [[ "${GMML_ROOT_DIR}" != "/programs/gems/gmml" ]]; then
    cp tests/correct_outputs/016.output_drawGlycan.txt tests/correct_outputs/016.output_drawGlycan.txt_backup
    sed -i -e "s#/programs/gems/gmml#${GMML_ROOT_DIR}#g" tests/correct_outputs/016.output_drawGlycan.txt
    #reverts our file back to normal once script exits and delets the backup
    trap 'cp -f tests/correct_outputs/016.output_drawGlycan.txt_backup tests/correct_outputs/016.output_drawGlycan.txt && rm tests/correct_outputs/016.output_drawGlycan.txt_backup' EXIT
fi

printf "Testing 016.test.DrawGlycan.cc..."
g++ -std=c++17 -I "${GMML_ROOT_DIR}"/ -I "${GMML_ROOT_DIR}"/ -L"${GMML_ROOT_DIR}"/bin/ -Wl,-rpath,"${GMML_ROOT_DIR}"/bin/ tests/016.test.DrawGlycan.cc -lgmml -pthread -o drawGlycan
./drawGlycan
>016.output_drawGlycan.txt
for dotFile in $(ls *.dot); do
    cat $dotFile >>016.output_drawGlycan.txt
    dotFileName="${dotFile%.*}"
    dot -Tsvg:cairo:cairo $dotFile -o $dotFileName.svg >/dev/null 2>&1
    rm $dotFile
done

for svgFile in $(ls *.svg); do
    cmp $svgFile tests/correct_outputs/016.output_SVGs/$svgFile
    if ! cmp $svgFile tests/correct_outputs/016.output_SVGs/$svgFile >/dev/null 2>&1; then
        printf "Test FAILED! Output file %s different to tests/correct_outputs/016.output_SVGs/%s\n" $svgFile $svgFile
        echo "Exit Code: 1"
        return 1
    fi
    rm $svgFile
done

if ! cmp 016.output_drawGlycan.txt tests/correct_outputs/016.output_drawGlycan.txt >/dev/null 2>&1; then
    printf "Test FAILED! Output file different\n"
    echo "Exit Code: 1"
    return 1
else
    printf "Test passed.\n"
    rm drawGlycan 016.output_drawGlycan.txt
    echo "Exit Code: 0"
    return 0
fi
