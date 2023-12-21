#!/bin/bash

printf "Testing 027.test.glycamResidueCombinator.cpp... "

GMML_ROOT_DIR=$(git rev-parse --show-toplevel)

if [[ "${GMML_ROOT_DIR}" != *"gmml" ]]; then
    echo -e "Test 027 failed, we think our GMML root directory is:\t${GMML_ROOT_DIR}\n"
    exit 1
fi

g++ -std=c++17 -I "${GMML_ROOT_DIR}"/ -L"${GMML_ROOT_DIR}"/bin/ -Wl,-rpath,"${GMML_ROOT_DIR}"/bin/ tests/027.test.glycamResidueCombinator.cpp -lgmml -pthread -o 027.glycamResiduecombinator.exe
outputFile=027.output.txt
./027.glycamResiduecombinator.exe ../dat/prep/GLYCAM_06j-1_GAGS_KDN.prep > $outputFile
file="GLYCAM_06k.lib"
if [ ! -f "${file}" ]; then
    echo -e "Test FAILED!\n ${file} does not exist\n"
    echo "Exit Code: 1"
     return 1
fi
if ! cmp -s "${file}" ../dat/CurrentParams/"${file}"; then
    echo -e "Test FAILED!\n ${file} is different from ../dat/CurrentParams/${file}\n"
    echo "Exit Code: 1"
     return 1
fi
  
sed -n '/\!entry/q;p' GLYCAM_06k.lib | sed 's/ "//g' | sed 's/"//g' | grep -v "index" >libEntries.txt
grep INT ../dat/prep/GLYCAM_06j-1_GAGS_KDN.prep | cut -d \  -f1 >prepEntriesList.txt
for id in `cat libEntries.txt`
do
    if ! grep -q "^$id" ../dat/prep/GLYCAM_06j-1_GAGS_KDN.prep;then
        echo "$id is new in lib file" >>$outputFile; 
    fi;
done

uniq -D libEntries.txt >>$outputFile

for id in `cat prepEntriesList.txt`; 
do 
    if ! grep -q " \"$id\"$" GLYCAM_06k.lib; 
    then 
        echo "$id not found in lib file" >>$outputFile  
    fi;     
done

if ! cmp -s $outputFile tests/correct_outputs/$outputFile; then
    echo -e "Test FAILED!\nOutput files different. Try:\ndiff $outputFile tests/correct_outputs/$outputFile"
    echo "Exit Code: 1"
    return 1    
fi
    
echo -e "Test passed.\n"
rm 027.glycamResiduecombinator.exe $outputFile libEntries.txt prepEntriesList.txt $file
echo "Exit Code: 0"

return 0
