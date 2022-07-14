#!/bin/bash

#NOTE: the -t flag is used to add our test files to the cmake file lists, this is important for linting.
INCLUDE_TEST_FILES=0
while getopts "t" option
do
        case "${option}" in
            t)
                INCLUDE_TEST_FILES=1
                ;;
			*)
                exit 1 
                ;;
        esac
done

FILES_UPDATED=0

echo "This script will walk through all the source and header files and put them into a file."
echo "This was done because it prevents us from having to adhere to some more annoying cmake"
echo "paradigms."
echo ""
#So I dont have to change up all teh paths after moving this dude into the scripts folder i just
#go up a dir. This does cause us to always have to call this script from within this directory,
#i.e. if we call this script in another script, the otehr scripts shell must cd into this scripts
#directory and then call this update script....
cd ..
#Grab all c++ files and chuck them into our file. We have the option to be able to also grab
#the tests c++ files in case we want to lint/run tests on em
echo "Hitting the list that contains all the paths to our source files"
if  [ "${INCLUDE_TEST_FILES}" == 1 ] || ! diff cmakeFileLists/cFileList.txt <(find "./src" -name "*.cc" -o -name "*.cpp" |  LC_ALL=C.UTF-8 sort) ; then
	if [ "${INCLUDE_TEST_FILES}" == 1 ]; then
        diff cmakeFileLists/cFileList.txt <(find "./src" "./tests" -name "*.cc" -o -name "*.cpp" |  LC_ALL=C.UTF-8 sort)
    fi
	if [ -f "./cmakeFileLists/cFileList.txt" ]; then
        echo "Previous data found, Removing old cc/cpp file list..."
        rm ./cmakeFileLists/cFileList.txt
    fi
    echo "Grabbing all *.cpp/*.cc files from /src/ and throwing into ./cmakeFileLists/cFileList.txt"
    if [ "${INCLUDE_TEST_FILES}" == 1 ]; then
        find "./src" "./tests" -name "*.cc" -o -name "*.cpp" |  LC_ALL=C.UTF-8 sort &> ./cmakeFileLists/cFileList.txt
    else
        find "./src" -name "*.cc" -o -name "*.cpp" |  LC_ALL=C.UTF-8 sort &> ./cmakeFileLists/cFileList.txt
    fi
    FILES_UPDATED=1
    echo "Finished fixing up the list that has all our source files"
else
    echo "No changes in our c file list."
fi
echo""

#Get all OUR header files, none of the external libraries. This is important because when cmake uses these file lists
#it will have all these files use the normal -I flag for gcc in the compile commands
echo "Hitting the list that contains all our paths to our directories that have our headers"
if ! diff cmakeFileLists/hDirectoryList.txt <( find "./includes" -type d ! -path "*External_Libraries*" |  LC_ALL=C.UTF-8 sort ) ; then
    if [ -f "./cmakeFileLists/hDirectoryList.txt" ]; then
        echo "Previous data found, Removing old header file directory path list..."
        rm ./cmakeFileLists/hDirectoryList.txt
    fi
    echo "Grabbing all dir paths to our h/hpp files and throwing into ./cmakeFileLists/hDirectoryList.txt"
    #Grabs all our header code dirs
    find "./includes" -type d ! -path "*External_Libraries*" |  LC_ALL=C.UTF-8 sort   &> ./cmakeFileLists/hDirectoryList.txt
    FILES_UPDATED=1
    echo "Finished fixing up the list that has all our include dirs"
else
    echo "No changes in our header directory list."
fi
echo ""

#Gets all our EXTERNAL lib headers. This is Eigen3 and pcg files. This is important because when cmake uses this file
#list to have all the files use the -isystem flag on the files which prevents linting on em and should help clean up
#our build outputs.
echo "Hitting the list that contains all our paths to the directories that have all our 3rd party libraries"
if ! diff cmakeFileLists/externalHDirectoryList.txt <(find "./includes" -type d -path "*External_Libraries*" |  LC_ALL=C.UTF-8 sort) ; then
    if [ -f "./cmakeFileLists/externalHDirectoryList.txt" ]; then
        echo "Previous data found, Removing old header file directory path list..."
        rm ./cmakeFileLists/externalHDirectoryList.txt
    fi
    echo "Grabbing all dir paths to our external libraries h/hpp files and throwing into ./cmakeFileLists/externalHDirectoryList.txt"
    #grabs only the header dirs that correspond to our 3rd party libs
    find . -type d -path "*External_Libraries*" |  LC_ALL=C.UTF-8 sort  &> ./cmakeFileLists/externalHDirectoryList.txt
    FILES_UPDATED=1
    echo "Finished fixing up the list that has all our external lib include dirs."
else
    echo "No changes in our external lib header directory list."
fi
echo ""


echo "Script completed."
echo ""
if [ 0 = "${FILES_UPDATED}" ]; then 
    echo "Files unchanged"
else
    echo "Files updated"
fi
