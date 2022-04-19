#!/bin/bash

INCLUDE_TEST_FILES=0
while getopts "t" option
do
        case "$option" in
                        t)
                                INCLUDE_TEST_FILES=1
                                ;;
			*)
                                exit 1 
                                ;;
        esac
done

#okay this is kinda lazy but if anything goes bad with these commands it will b
#thrown in the given files
FILES_UPDATED=0
if  ! diff cmakeFileLists/cFileList.txt <(find ./src/ -name \*.cc -o -name \*.cpp | sort) || [ ${INCLUDE_TEST_FILES} == 1 ] ; then
	if [ -f "./cmakeFileLists/cFileList.txt" ]; then
        printf "\nPrevious data found, Removing old cc/cpp file list..."
        rm ./cmakeFileLists/cFileList.txt
    fi
    printf "\nGrabbing all *.cpp/*.cc files from /src/ and throwing into ./cmakeFileLists/cFileList.txt\n"
    if [ ${INCLUDE_TEST_FILES} == 1 ]; then
        find ./src/ ./tests/ -name *.cc -o -name \*.cpp | sort &> ./cmakeFileLists/cFileList.txt
    else
        find ./src/ -name *.cc -o -name *.cpp | sort &> ./cmakeFileLists/cFileList.txt
    fi
    FILES_UPDATED=1
fi

if ! diff cmakeFileLists/hDirectoryList.txt <(find ./includes/ -type d -path "*" | sort) ; then
    if [ -f "./cmakeFileLists/hDirectoryList.txt" ]; then
        printf "\nPrevious data found, Removing old header file directory path list..."
        rm ./cmakeFileLists/hDirectoryList.txt
    fi
    printf "\nGrabbing all dir paths to our h/hpp files and throwing into ./cmakeFileLists/hDirectoryList.txt\n"
    #for when TRUs fixed
    #find ./includes/ -type d ! -path "*External_Libraries*" | sort 
    find ./includes/ -type d -path "*"  | sort  &> ./cmakeFileLists/hDirectoryList.txt
    FILES_UPDATED=1
fi

#Will be used once we get our TRUs fixed. This is so we can denote everything in the external libs 
#folder as a system file, thus it wont be linted etc.
#if ! diff cmakeFileLists/externalHDirectoryList.txt <(find . -type d -path "*External_Libraries*" | sort) ; then
#    if [ -f "./cmakeFileLists/externalHDirectoryList.txt" ]; then
#        printf "\nPrevious data found, Removing old header file directory path list..."
#        rm ./cmakeFileLists/externalHDirectoryList.txt
#    fi
#    printf "\nGrabbing all dir paths to our external libraries h/hpp files and throwing into ./cmakeFileLists/externalHDirectoryList.txt\n"
#    #find ./includes/ -type d ! -path "*External_Libraries*" | sort 
#    find . -type d -path "*External_Libraries*" | sort  &> ./cmakeFileLists/externalHDirectoryList.txt
#    FILES_UPDATED=1
#fi

if [ 0 = ${FILES_UPDATED} ]; then 
    printf "\nFiles unchanged\n"
else
    printf "\nFiles updated\n"
fi
