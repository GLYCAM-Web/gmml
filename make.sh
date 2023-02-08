#!/bin/bash

##### PRETTY PRINT STUFF #####
#lazy and dont want to have to type all these color variables a bunch
RESET_STYLE='\033[0m'
PASSED_STYLE='\033[0;32m\033[1m'
INFO_STYLE='\033[0;33m\033[1m'
ERROR_STYLE='\033[0;31m\033[1m'
#variable that will be used to check if we are trying to run any linting
#on our test files
TESTIN_TIME=0

################################################################
#########                 FUNCTIONS                    #########
################################################################

get_numprocs() {
	if [ -z "${GEMSMAKEPROCS}" ]; then
		NMP=4
	elif ! [[ ${GEMSMAKEPROCS} =~ ^[0-9]+$ ]] ; then
		echo "Warning:  GEMSMAKEPROCS is not a valid integer; setting to 4"
		NMP=4
	elif [ "${GEMSMAKEPROCS}" -eq "0" ] ; then
		echo "Warning:  GEMSMAKEPROCS cannot be zero; setting to 4"
		NMP=4
	else
		NMP=${GEMSMAKEPROCS}
	fi
}

check_gemshome() {
	if [ -z "${GEMSHOME}" ]; then
		echo ""
		echo "Error:  GEMSHOME environment variable is not set! It should be set to"
		echo "$1"
		exit 1
	elif [ ! -d "${GEMSHOME}" ]; then
		echo ""
		echo "Error:  GEMSHOME environment variable is set to ${GEMSHOME} -- this does"
		echo "not appear to be a directory. It should be set to"
		echo "$1"
		exit 1
	fi
}

#this function is needed to help block people from running code if they have
#not updated the file lists that cmake uses and they are trying to run make
#this is important because it allows us to continue to use my workaround
#which allows us to use the ease of file globbing but without being bums
#and doing bad cmake practices where it runs ops we dont know
checkCMakeFileLists()
{
    echo -e "${INFO_STYLE}###### Checking our cmake file lists to make sure they match files in our codebase ######"
    BAD_FILE_LISTS=0
    #check our cfiles stuff
    echo -e "${RESET_STYLE}\tEnsuring the list of all c/cpp files, cFileList.txt, matches\n\tour current source code${INFO_STYLE}"
    if  [ "${TESTIN_TIME}" == 1 ] && ! diff cmakeFileLists/cFileList.txt <(find "./src" "./tests" -name "*.cc" -o -name "*.cpp" |  LC_ALL=C.UTF-8 sort); then
        echo -e "\tWarning: you want to build for testing but the cfiles isnt"
        echo -e "\tupdated with the needed files. Read the readme. Script will exit."
        BAD_FILE_LISTS=1
    elif ! diff cmakeFileLists/cFileList.txt <(find "./src" -name "*.cc" -o -name "*.cpp" |  LC_ALL=C.UTF-8 sort); then
        echo -e "\tWarning: you want to build but the cfiles isnt updated with"
        echo -e "\tthe needed files. Read the readme. Script will exit."
        BAD_FILE_LISTS=1
    fi
    #check our hdirectory stuff
    echo -e "\n${RESET_STYLE}\tEnsuring the list of all header directories, hDirectoryList.txt,\n\tmatches our current source code${INFO_STYLE}"
    if ! diff cmakeFileLists/hDirectoryList.txt <( find "./includes" -type d ! -path "*External_Libraries*" |  LC_ALL=C.UTF-8 sort) ; then
        echo -e "\tWarning: you want to build for testing but the hdirectory list file"
        echo -e "\tisnt updated with the needed files. Read the readme. Script will exit."
        BAD_FILE_LISTS=1
    fi
    #check our external lib stuff
    echo -e "${RESET_STYLE}\n\tEnsuring the list of all external lib header directories, \n\texternalHDirectoryList.txt, matches our current source code${INFO_STYLE}"
    if ! diff cmakeFileLists/externalHDirectoryList.txt <(find "./includes" -type d -path "*External_Libraries*" |  LC_ALL=C.UTF-8 sort) ; then
        echo -e "\tWarning: you want to build for testing but the list of directories"
        echo -e "\tfor our external libs isnt updated with the needed files. Read" 
        echo -e "\tthe readme. Script will exit."
        BAD_FILE_LISTS=1
    fi
    
    if [ "${BAD_FILE_LISTS}" == 1 ]; then
        echo -e "${ERROR_STYLE}\n==============================================================="
        echo "READ THE README, IT WILL EXPLAIN WHAT IS GOING ON."
        echo "Check under the \"compiling the library\" seciont. Be sure to read"
        echo "that section. This is a non-traditional workaround to help with tooling"
        echo "and to help with making CMAKE less horrid to deal with. That being"
        echo "said, we must know what this workaround is, how it works, and if we"
        echo "do not know this the purpose of the workaround is defeated thus we"
        echo "we will be forced to remove this workaround. Exiting."
        echo -e "===============================================================${RESET_STYLE}\n"
        exit 1
    fi
    echo -e "${PASSED_STYLE}CMake file lists are good, lets run it${RESET_STYLE}"
   
}

################################################################
##########               Print help message         ############
################################################################

printHelp()
{
echo -e "\n\t========= GMML MAKE.SH SCRIPT =========
If you are a user then this $0 is not for you:
go to the GEMS home dir; If needed GEMS $0 will 
spawn GMML building.

This is for building GMML in isolation! You can wrap GMML here
using swig but interfacing with the wrapped lib will be up to you

GEMSHOME should be set to the parent of the gmml directory.
After GMML is built the next step is testing, first
run cd ./tests/, then run ./compile_run_tests.bash
*************************************************************
Options are as follows:
\t-c\t\t\tClean all files from previous builds
\t-j <NUM_JOBS>\t\tBuild GMML with <NUM_JOBS>
\t-o <O0/O2/OG/debug>\tBuild GMML using no optimization, 2nd \n\t\t\t\t\tlevel optimization, or with debug symbols
\t-w\t\t\tWrap gmml in python using swig
\t-t\t\t\tChange up our lists that index our source code to also
\t\t\t\t\tindex our test files, needed for running autotooling
\t-h\t\t\tPrint this help message and exit
*************************************************************
Exiting."
exit 1
}

################################################################
#########                CHECK SETTINGS                #########
################################################################

gemshome=$(pwd)
check_gemshome "${gemshome}"
get_numprocs

################################################################
#########              CREATE CLIENT HOOKS             #########
################################################################

#Cannot use server side hooks when hosting on git-hub.
#Stuff in .git/hooks is ignored by git.
#Solution: The folder .hooks is tracked by git.
# Copy .hooks to .git/hooks during installation.
cp -r "${GEMSHOME}"/gmml/.hooks/* "${GEMSHOME}"/gmml/.git/hooks/
#I don't think this is ideal, and is perhaps silly. OG Apr 2017.

################################################################
#########                SET UP DEFAULTS               #########
################################################################

#Changing the following options are for developers.
CLEAN="0"
CMAKE_BUILD_TYPE_FLAG="-DCMAKE_BUILD_TYPE=Release"
#target gmml_wrapped = swig wrapped gmml
#target gmml = just gmml with no wrap
MAKE_TARGET="gmml"

################################################################
#########               COMMAND LINE INPUTS            #########
################################################################

#Proper cmake use is to use the build_type flag. We mess with specifics in the
# actual cmakelists.txt file. Listed below are the different build types in
# normal english. Ty stack overflow

#1. Release: high optimization level, no debug info, code or asserts.
#2. Debug: No optimization, asserts enabled, [custom debug (output) code enabled],
#   debug info included in executable (so you can step through the code with a
#   debugger and have address to source-file:line-number translation).
#3. RelWithDebInfo: optimized, *with* debug info, but no debug (output) code or asserts.
#4. MinSizeRel: same as Release but optimizing for size rather than speed.

#Below are the compiler flags passed for each build type before they were
# overloaded in our cmakelists.txt
#1. Release: `-O3 -DNDEBUG`
#2. Debug: `-O0 -g`
#3. RelWithDebInfo: `-O2 -g -DNDEBUG`
#4. MinSizeRel: `-Os -DNDEBUG`

#Below are what we changed the compiler flags to
#1. Release: `-O2 -DNDEBUG`
#2. Debug: `-O0 -g`
#3. RelWithDebInfo: `-O0 -DNDEBUG`	(This is our non-optimized build)
#4. MinSizeRel: `` 					(NOT USED)
#Note that our debug symbols are only available with non-optimized flags, this is
#intentional.

# Please refer to https://blog.feabhas.com/2021/07/cmake-part-1-the-dark-arts/
# Follow the paradigm
#check the compile_run_tests.sh file for a description of what this stuff is
while getopts "j:o:cwht" option
do
	case "${option}" in
			j)
				jIn="${OPTARG}"
				if [[ "${jIn}" =~ ^[1-9][0-9]*$ ]]; then
					NMP="${jIn}"
				else
					printHelp
				fi
				;;
			o)
				oIn="${OPTARG}"
				if [ "${oIn}" == "O0" ] || [ "${oIn}" == "no_optimize" ]; then
					CMAKE_BUILD_TYPE_FLAG="-DCMAKE_BUILD_TYPE=RelWithDebInfo"
				elif [ "${oIn}" == "O2" ] || [ "${oIn}" == "optimize" ]; then
					CMAKE_BUILD_TYPE_FLAG="-DCMAKE_BUILD_TYPE=Release"
				elif [ "${oIn}" == "debug" ] || [ "${oIn}" == "OG" ]; then
					CMAKE_BUILD_TYPE_FLAG="-DCMAKE_BUILD_TYPE=Debug"
				else
					printHelp
				fi
				;;
			c)
				CLEAN=1
				;;
			t)
                TESTIN_TIME=1
                ;;
			w)
				MAKE_TARGET="gmml_wrapped"
				;;
            h)
                printHelp
                ;;
			*)
				printHelp
				;;
	esac
done

echo -e "\n${INFO_STYLE}###### Starting GMML Build Script ######${RESET_STYLE}"
START_TIME=$(date +%s)
echo -e "Steps this script takes:
\t1) Check the indexed file lists present in cmakeFileLists are correct
\t2) If we want to aggressively clean all traces of GMML, do so
\t3) Generate the Makefile, and build directory, for GMML by using CMake
\t4) cd to generated dir, then use the generated Makefile with our specified target
\t5) Script completed
"
#check our file lists before we do anything 
checkCMakeFileLists

echo -e "\n${INFO_STYLE}###### Running GMML make.sh with these settings ######${RESET_STYLE}
Build Jobs:\t${NMP}
Build Type:\t${CMAKE_BUILD_TYPE_FLAG}
Make Target:\t${MAKE_TARGET}
Agro Clean:\t${CLEAN}
Hit Test Files:\t${TESTIN_TIME}
###############################################"
################################################################
#########                  COMPILE GMML                #########
################################################################

if [ "${CLEAN}" == "1" ]; then
	echo -e "\n${INFO_STYLE}###### Aggressive Cleaning GMML. Please note this is VERY aggressive! ######${RESET_STYLE}"
	if [ -d "./cmakeBuild" ]; then
		rm -rf ./cmakeBuild
	fi
	if [ -d "./lib" ]; then
		rm -rf ./lib
	fi
	echo -e "${PASSED_STYLE}Aggressive clean has been completed, all compiled/wrapped\n\tGMML code has been completely removed from this dir${RESET_STYLE}"
fi

#Note that we have to generate our makefile with cmake before we build
# The -S (source) flag is for where the cmakelists file is located
# and the -B (build dir) flag is for where everything will be built to
echo -e "\n${INFO_STYLE}###### Generating GMML makefile using cmake ######${RESET_STYLE}"
#The -D<VARIABLE> stuff passes the variable to cmake. This is done instead of
#having to deal with setting env vars to accomplish the same goal. Basically just tryna
#make life easier. Also this command ends up generating the Makefile into the cmakeBuild
#directory which is normal behavior for cmake. Also all our libraries etc. will be thrown
#in that directory.
cmake "${CMAKE_BUILD_TYPE_FLAG}" -S . -B ./cmakeBuild -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ \
|| { echo -e "${ERROR_STYLE}ERROR RUNNING CMAKE ON GMML: $0 FAILED, EXITING${RESET_STYLE}" ; exit 1; }
echo -e "${PASSED_STYLE}All GMML build files have been successfully generated by CMake${RESET_STYLE}"

#NOTE: All our build stuff is within the cmakeBuild dir
echo -e "\n${INFO_STYLE}###### Build files generated, now actually making GMML in ./gmml/cmakeBuild/ ######"
echo -e "###### COMMAND USED: ./gmml/cmakeBuild$ make -j${NMP} ${MAKE_TARGET} ######${RESET_STYLE}"

cd cmakeBuild || { echo -e "${ERROR_STYLE}ERROR BUILDING GMML: $0 FAILED, EXITING${RESET_STYLE}" ; exit 1; }

#actually use the makefile that cmake generated
make -j"${NMP}" "${MAKE_TARGET}" \
|| { echo -e "${ERROR_STYLE}ERROR BUILDING GMML: $0 FAILED, EXITING${RESET_STYLE}" ; exit 1; }
cd ..

################################################################
#########              WRAP UP TO GEMS                 #########
################################################################

STUPID_WRAPPED_MSG=""
#Wrapping is handled by cmake
if [ "${MAKE_TARGET}" == "gmml_wrapped" ]; then
    STUPID_WRAPPED_MSG=" and wrapping"
fi
echo -e "${PASSED_STYLE}
GMML compilation${STUPID_WRAPPED_MSG} is finished at $(date).
Time Taken: $(( $(date +%s) - START_TIME )) seconds${RESET_STYLE}"
exit 0

