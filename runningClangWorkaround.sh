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

#run-clang-tidy -p ./cmakeBuild/ -header-filter=\.hpp$ -fix

echo "I know its rough, but this will be used until we finish fixing up our includes"
if [ ${AUTO_RUN} == 1  ] && [ ! -d "./cmakeBuild" ]; then
	./make.sh -j $(grep -c ^processor /proc/cpuinfo) -c	
fi	

if [ -d "./cmakeBuild" ]; then
    #gmml/src
    cd src/
    echo "Creating a cpp file that contains all of our source code headers. This will not include any"
    echo "external library headers because we dont want to lint code that isnt ours."
    find ../includes -iname "*.hpp" -o -iname "*.h" ! -path "*External_Libraries*" | sort > allHeaders.cpp
    #have to delete a couple of dirs in there cause the find command doesnt work completely its wonky
    #delete all lines that include External_Libraries
    sed -i '/External_Libraries/d' ./allHeaders.cpp
    #add #include " to the beginning of each line
    sed -i -e 's/^/\#include "/' ./allHeaders.cpp
    #add closing quotes for each line
    sed -i -e 's/.hpp/\.hpp\"/' ./allHeaders.cpp
    #gmml
    cd ../
    echo "Updating the file list cmake uses to build. It is also including out test files so we can"
    echo "Lint our test files."
    ./updateCmakeFileList.sh -t
    echo "Rebuilding our tool chain so it has the new allHeaders.cpp known"
    #regens our makefile and the compile commands file. There is no need
    #to fully rebuild out codebase. I may be wrong but so far it seems like we are good.
    cmake ${CMAKE_BUILD_TYPE_FLAG} -S . -B ./cmakeBuild -DWRAP_GMML=${WRAP}
    
    #Below is no longer needed because we filter our lint by file extension. Keep in mind we currently have
    #all external library headers end with .h and all our source code end with .hpp
    #Going to keep the code cause it may be nice to use in the future.
    #gmml/cmakeBuild
    #cd cmakeBuild
    #mv compile_commands.json backup_compile_commands.json
    #basically making the compile commands json only include the allheaders cpp file. This is so 
    #clang will use the all headers file and only that file to walk all of our headers :o) 
    #this is to patch bad clang behvior due to some bad includes that we will eventually find
    #grep '\"file\":.*allHeaders.cpp' backup_compile_commands.json -B 3 -A 1 > compile_commands.json
    #sed -i '1s/^/[\n/' compile_commands.json 
    #sed -i '$a]' compile_commands.json
    #sed -i "s/},/}/" compile_commands.json
    #echo "Compile Commands workaround created, running clang so we can hit all our headers"
    # gmml
    #cd ..
    
    #now actually run clang on all our codes, this should run it on everything we have. This needs to run a
    #couple of times because tool limitations. Don't worry, this brute force wont be here forever. I will
    #eventually do a pattern like whats in the updateCmakeFileList script
    run-clang-tidy $(find ./src ./tests -type f -name "*.cpp") -p ./cmakeBuild -header-filter=\.hpp$ -fix
    run-clang-tidy $(find ./src ./tests -type f -name "*.cpp") -p ./cmakeBuild -header-filter=\.hpp$ -fix
    run-clang-tidy $(find ./src ./tests -type f -name "*.cpp") -p ./cmakeBuild -header-filter=\.hpp$ -fix
    run-clang-tidy $(find ./src ./tests -type f -name "*.cpp") -p ./cmakeBuild -header-filter=\.hpp$ -fix
    
    echo "Command that was run:"
    echo "run-clang-tidy \$(find ./src ./tests -type f -name \"*.cpp\") -p ./cmakeBuild -header-filter=\.hpp$ -fix"
                
    #Now check if we want to revert back to nromal gmml. 
    echo "This will delete the allHeaders.cpp and redo your make tooling and rebuild the lib. Basically this is if"
    echo "you just want to have the script be done. This will remake your code."
    echo ""
    read -p "Revert GMML back to normal? [y/n]" -n 1 -r
    echo ""
    if [ ${AUTO_RUN} == 1  ] || [[ $REPLY =~ ^[Yy]$ ]] ; then
        echo "Removing src/allHeaders.cpp file"
        rm ./src/allHeaders.hpp
        echo "Updating file lists cmake uses to walk our files."
        ./updateCmakeFileList.sh
        #Possibly do a clean compile? It honsestly isnt necessary from my understanding.
        echo "Refreshing our toolchain"
        cmake ${CMAKE_BUILD_TYPE_FLAG} -S . -B ./cmakeBuild -DWRAP_GMML=${WRAP}
        #./make.sh -j $(grep -c ^processor /proc/cpuinfo) -c 
        echo "Completed"
        exit 0
    fi
    echo "Okay header files have been run on, workaround complete, leaving everything set up to deal with headers. Once you are"
    echo "done with the header workaround file rebuild you make toolchain with:"
    echo "cmake ${CMAKE_BUILD_TYPE_FLAG} -S . -B ./cmakeBuild"
    echo "Now if you want to run clang on your src files go ahead. Good luck"
    AUTO_RUN=0
    exit 0    
fi
echo "run make.sh first you at least need the makefile"
exit 1



