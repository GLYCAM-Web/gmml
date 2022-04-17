#!/bin/bash

AUTO_RUN=0
while getopts "a" option
do
        case "$option" in
                        a)
                                AUTO_RUN=1
                                ;;
			*)
                                exit 1 
                                ;;
        esac
done


echo "I know its rough, but this will be used until we finish fixing up our includes"
if [ ${AUTO_RUN} == 1  ] && [ ! -d "./cmakeBuild" ]; then
	./make.sh -j $(grep -c ^processor /proc/cpuinfo) -c	
fi	


if [ -d "./cmakeBuild" ]; then
    #Dont really want to do a bunch of regex now. In the future I will.
    #gmml/src
    cd src/
    echo "Creating a cpp file that contains all of our source code headers"
    find ../includes -iname "*.hpp" -o -iname "*.h" ! -path "*External_Libraries*" | sort > allHeaders.cpp
    #have to delete a couple of dirs in there cause the find command doesnt work completely its wonky
    sed -i '/External_Libraries/d' ./allHeaders.cpp
    sed -i -e 's/^/\#include "/' ./allHeaders.cpp 
    sed -i -e 's/.hpp/\.hpp\"/' ./allHeaders.cpp
    #gmml
    cd ../
    echo "Updating the file list cmake uses to build"
    ./updateCmakeFileList.sh
    #now we have to remove all the eigen stuff from our include list cause we dont want traversed at all
    sed -i '/External_Libraries/d' ./cmakeFileLists/hDirectoryList.txt
    echo "Rebuilding our tool chain so it has the new allHeaders.cpp known"
    #cause cmake cache
    rm -rf ./cmakeBuild
    cmake ${CMAKE_BUILD_TYPE_FLAG} -S . -B ./cmakeBuild -DWRAP_GMML=${WRAP}
    #gmml/cmakeBuild
    cd cmakeBuild
	
    mv compile_commands.json backup_compile_commands.json
    #basically making the compile commands json only include the allheaders cpp file. This is so 
    #clang will use the all headers file and only that file to walk all of our headers :o) 
    #this is to patch bad clang behvior due to some bad includes that we will eventually find
    grep '\"file\":.*allHeaders.cpp' backup_compile_commands.json -B 3 -A 1 > compile_commands.json
    sed -i '1s/^/[\n/' compile_commands.json 
    sed -i '$a]' compile_commands.json
    sed -i "s/},/}/" compile_commands.json
    echo "Compile Commands workaround created, running clang so we can hit all our headers"
    # gmml
    cd ..
    #now actually run clang on our header files
    run-clang-tidy -p ./cmakeBuild -header-filter=.* -fix
        
    #Now to ask if we wanna keep the header walking compile commands, just in case. If not I just 
    #revert you back to your normal gmml status
    echo "This will delete the allHeaders.cpp and redo your make tooling and rebuild the lib. Basically this is if"
    echo"you just want to have the script be done. This will remake your code."
    echo ""
    read -p "Revert GMML back to normal? [y/n]" -n 1 -r
    echo ""
    if [ ${AUTO_RUN} == 1  ] || [[ $REPLY =~ ^[Yy]$ ]] ; then
        echo "Removing src/allHeaders.cpp file"
        rm ./src/allHeaders.hpp
        echo "Redoing filelists"
        ./updateCmakeFileList.sh
        echo "refreshing toolchain"
        cmake ${CMAKE_BUILD_TYPE_FLAG} -S . -B ./cmakeBuild -DWRAP_GMML=${WRAP}
        #./make.sh -j $(grep -c ^processor /proc/cpuinfo) -c 
        echo "Completed"
        exit 0
    fi
    echo "Okay header files have been run on, workaround complete, leaving everything set up to deal with headers. Just regen"
    echo "the build tooling using:"
    echo "cmake ${CMAKE_BUILD_TYPE_FLAG} -S . -B ./cmakeBuild -DWRAP_GMML=${WRAP}"
    #cmake ${CMAKE_BUILD_TYPE_FLAG} -S . -B ./cmakeBuild -DWRAP_GMML=${WRAP}
    echo "Now if you want to run clang on your src files go ahead. Good luck"
    echo "if you want to run clang tools on your c files just refresh the build tooling"
    echo "and then run:"
    echo "run-clang-tidy -p ./cmakeBuild -header-filter=.* -fix"
    echo "also maybe remove the build dir"
    AUTO_RUN=0
    exit 0    
fi
echo "run make.sh first you at least need the makefile"
exit 0



