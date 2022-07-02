#!/bin/bash

##### PRETTY PRINT STUFF #####
#lazy and dont want to have to type all these color variables a bunch
RESET_STYLE='\033[0m'
PASSED_STYLE='\033[0;32m\033[1m'
INFO_STYLE='\033[0;33m\033[1m'
ERROR_STYLE='\033[0;31m\033[1m'
#variable that will be used to check if we are trying to run any linting
#on our test files
AUTO_TESTIN_TIME=0

gemshome=$(pwd)

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
    #if we are doing the auto testing stuff just fix up our lists auto
    if [ "${AUTO_TESTIN_TIME}" == 1 ]; then
        echo -e "${INFO_STYLE}###### Auto running the cmake file list update to include test files ######${RESET_STYLE}"
        mkdir ./.tempCmakeLists
        cp ./cmakeFileLists/* ./.tempCmakeLists/
        trap 'cp -fv ./.tempCmakeLists/* ./cmakeFileLists/ && rm -rfv ${gemshome}/.tempCmakeLists && echo -e "${INFO_STYLE}###### Cmake file lists restored to state before script was run ######${RESET_STYLE}"' EXIT
        ./updateCmakeFileList.sh -t || { echo -e "${ERROR_STYLE}ERROR AUTO UPDATING CMAKEFILELISTS, EXITING${RESET_STYLE}" ; exit 1; }
        echo -e "${PASSED_STYLE}###### Auto updating cmakefilelists to include tests succeded ######${RESET_STYLE}\n"
    fi
    
    echo -e "${INFO_STYLE}###### Checking our cmake file lists to make sure they match files in our codebase ######"
    BAD_FILE_LISTS=0
    #check our cfiles stuff
    echo -e "${RESET_STYLE}\tEnsuring the list of all c/cpp files, cFileList.txt, matches\n\tour current source code${INFO_STYLE}"
    if  [ "${AUTO_TESTIN_TIME}" == 1 ] && ! diff cmakeFileLists/cFileList.txt <(find "./src" "./tests" -name "*.cc" -o -name "*.cpp" |  LC_ALL=C.UTF-8 sort); then
        echo -e "\n\tWarning: you want to build for testing but the cfiles isnt"
        echo -e "\tupdated with the needed files. Read the readme. Script will exit."
        BAD_FILE_LISTS=1
    elif  [ "${AUTO_TESTIN_TIME}" == 0 ] && ! diff cmakeFileLists/cFileList.txt <(find "./src" -name "*.cc" -o -name "*.cpp" |  LC_ALL=C.UTF-8 sort); then
        echo -e "\n\tWarning: you want to build but the cfiles isnt updated with"
        echo -e "\tthe needed files. Read the readme. Script will exit."
        BAD_FILE_LISTS=1
    fi
    #check our hdirectory stuff
    echo -e "\n${RESET_STYLE}\tEnsuring the list of all header directories, hDirectoryList.txt,\n\tmatches our current source code${INFO_STYLE}"
    if ! diff cmakeFileLists/hDirectoryList.txt <( find "./includes" -type d ! -path "*External_Libraries*" |  LC_ALL=C.UTF-8 sort) ; then
        echo -e "\n\tWarning: you want to build but the hdirectory list file"
        echo -e "\tisnt updated with the needed files. Read the readme. Script will exit."
        BAD_FILE_LISTS=1
    fi
    #check our external lib stuff
    echo -e "${RESET_STYLE}\n\tEnsuring the list of all external lib header directories, \n\texternalHDirectoryList.txt, matches our current source code${INFO_STYLE}"
    if ! diff cmakeFileLists/externalHDirectoryList.txt <(find "./includes" -type d -path "*External_Libraries*" |  LC_ALL=C.UTF-8 sort) ; then
        echo -e "\n\tWarning: you want to build for but the list of directories"
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
    echo -e "${PASSED_STYLE}###### CMake file lists are good, lets run it ######${RESET_STYLE}"
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
\t-o <O0/O2/OG/debug>\tBuild GMML using no optimization, 2nd
\t\t\t\t\tlevel optimization, or with debug symbols
\t-w\t\t\tWrap gmml in python using swig
\t-h\t\t\tPrint this help message and exit
*************************************************************
Exiting."
exit 1
}

################################################################
#########        SOME DUMB FUNCTIONS                   #########
################################################################
optionsBorked()
{
    [ "$#" == 2 ] || { echo -e "${ERROR_STYLE}ERROR OPTIONS BORKED NOT USED CORRECT${RESET_STYLE}"; exit 1; }
    case $2 in
        badArg) 
                echo -e "${ERROR_STYLE}ERROR: OPTION -${1} PASSED INCORRECT ARGUMENT${RESET_STYLE}"
                printHelp
                ;;
        notArg)
                echo -e "${ERROR_STYLE}ERROR: OPTION -${1} DOES NOT ALLOW AN ARGUMENT${RESET_STYLE}"
                printHelp
                ;;
        *)
                echo -e "${ERROR_STYLE}ERROR OPTIONS BORKED HIT WILD CARD ${RESET_STYLE}"
                exit 1
                ;;
    esac
}

################################################################
#########                CHECK SETTINGS                #########
################################################################
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
while getopts "j:o:cwht:" option
do
	case "${option}" in
			j)
				jIn="${OPTARG}"
				if [[ "${jIn}" =~ ^[1-9][0-9]*$ ]]; then
					NMP="${jIn}"
				else
					optionsBorked "${option}" "badArg"
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
					optionsBorked "${option}" "badArg"
				fi
				;;
			c)
                #to make sure we arent passing anything extra to this, we use a 
                #indirect expansion of optind, which contains next index of
                #something interesting aka character. This takes the index, gets
                #the variable it relates to, then 0 is the index we start at
                #aka offset then we have length of 1 cause we only care about next
                #character to ensure next arg is a flag. Remnant of trying to get optional multiflag passing
                #stuff where u can have a flag w/wout arguments
                { [ -z "${!OPTIND:0:1}" ] || [ "${!OPTIND:0:1}" == "-" ]; } || { optionsBorked "${option}" "notArg"; }
				CLEAN=1
				;;
			t)
                tIn="${OPTARG}"
                if [ "${tIn}" == "auto" ]; then
                    AUTO_TESTIN_TIME=1
                else
                    optionsBorked "${option}" "badArg"
                fi    
                ;;
			w)
                { [ -z "${!OPTIND:0:1}" ] || [ "${!OPTIND:0:1}" == "-" ]; } || { optionsBorked "${option}" "notArg"; }
				MAKE_TARGET="gmml_wrapped"
				;;
            h)
                printHelp
                ;;
			*)
                echo -e "${ERROR_STYLE}ERROR INCORRECT FLAG PASSED${RESET_STYLE}"
				printHelp
				;;
	esac
done

START_TIME=$(date +%s)

echo -e "\n${INFO_STYLE}###### Running GMML make.sh with these settings ######${RESET_STYLE}
Build Jobs:\t${NMP}
Build Type:\t${CMAKE_BUILD_TYPE_FLAG}
Make Target:\t${MAKE_TARGET}
Agro Clean:\t${CLEAN}"

if [ "${AUTO_TESTIN_TIME}" == 1 ]; then
    echo -e "Auto test prep: ${AUTO_TESTIN_TIME}"
fi

echo -e "Steps this script takes:
\t1) Check the indexed file lists present in cmakeFileLists are correct
\t2) If we want to aggressively clean all traces of GMML, do so
\t3) Generate the Makefile, and build directory, for GMML by using CMake
\t4) cd to generated dir, then use the generated Makefile with our specified target
\t5) Script completed
###############################################
"
#check our file lists before we do anything 
checkCMakeFileLists

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
	echo -e "${PASSED_STYLE}###### All compiled/wrapped GMML has been completely removed ######${RESET_STYLE}"
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
echo -e "${PASSED_STYLE}###### All GMML build files have been successfully generated by CMake ######${RESET_STYLE}"

#Now if we are just trying to run some linting etc. using tools on our codebase
#and want to also check out the testing source code files, we know that the whole
#thing will not compile, thus we reset the indexed file lists to their original state.
#and exit without compiling or anything fun
if [ "${AUTO_TESTIN_TIME}" == 1 ]; then
    echo -e "${INFO_STYLE}###### Restoring cmake file lists to their original state ######${RESET_STYLE}"
    #trap func takes care of the restore
    exit 0
fi

#NOTE: All our build stuff is within the cmakeBuild dir
echo -e "\n${INFO_STYLE}###### Now actually making GMML, COMMAND USED: ./gmml/cmakeBuild$ make -j${NMP} ${MAKE_TARGET} ######${RESET_STYLE}"
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

