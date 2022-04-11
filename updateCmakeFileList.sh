#!/bin/bash

#okay this is kinda lazy but if anything goes bad with these commands it will b
#thrown in the given files

if [ -f "./cmakeFileLists/cFileList.txt" ]; then
        printf "\nPrevious data found, Removing old cc/cpp file list..."
        rm ./cmakeFileLists/cFileList.txt
fi
printf "\nGrabbing all *.cpp/*.cc files from /src/ and throwing into ./cmakeFileLists/cFileList.txt\n"
find ./src/ -name \*.cc -o -name \*.cpp &> ./cmakeFileLists/cFileList.txt
#lmk if we want to do a diff type thying to show what we are adding

if [ -f "./cmakeFileLists/hDirectoryList.txt" ]; then
	printf "\nPrevious data found, Removing old header file directory path list..."
	rm ./cmakeFileLists/hDirectoryList.txt
fi
printf "\nGrabbing all dir paths to our h/hpp files and throwing into ./cmakeFileLists/hDirectoryList.txt\n"
find ./includes/ -type d &> ./cmakeFileLists/hDirectoryList.txt

printf "\nFiles updated\n"
