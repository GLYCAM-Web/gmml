#!/bin/bash

GMML_ROOT_DIR=$(git rev-parse --show-toplevel)

if [ "$(git config --get remote.origin.url)" != "https://github.com/GLYCAM-Web/gmml" ]; then
            exit 1
fi


printf "Testing 017.test.GlycoproteinBuilder.cpp... "
g++ -std=c++17 -I "${GMML_ROOT_DIR}" -L"${GMML_ROOT_DIR}"/bin/ -Wl,-rpath,"${GMML_ROOT_DIR}"/bin/ "${GMML_ROOT_DIR}"/internalPrograms/GlycoproteinBuilder/main.cpp -lgmml -pthread -o gpBuilder
./gpBuilder tests/inputs/017.GlycoproteinBuilderInput.txt > output_GlycoproteinBuilder.txt 2>&1
fileList=("GlycoProtein_All_Resolved.pdb" "GlycoProtein_All_Resolved.off" "GlycoProtein_All_Resolved_Serialized.pdb" "output_GlycoproteinBuilder.txt")
for file in ${fileList[@]}; 
do
  	if [ -f $file ]; then
  	    if ! cmp $file tests/correct_outputs/017.$file  > /dev/null 2>&1; then
  	        printf "Test FAILED!\n $file is different from tests/correct_outputs/017.$file\n"
            echo "Exit Code: 1"
            return 1
            exit 1
        else
            rm $file    
        fi
    else
        printf "Test FAILED!\n $file does not exist\n"
        echo "Exit Code: 1"
        return 1
        exit 1
    fi      
done	
printf "Test passed.\n"
rm gpBuilder  
echo "Exit Code: 0"
return 0
exit 0
