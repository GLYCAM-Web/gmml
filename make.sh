#!/bin/bash

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
    echo "Checking our file lists to make sure you are good to compile"
    BAD_FILE_LISTS=0
    #check our cfiles stuff
    if  [ "${TESTIN_TIME}" == 1 ] && ! diff cmakeFileLists/cFileList.txt <(find "./src" "./tests" -name "*.cc" -o -name "*.cpp" |  LC_ALL=C.UTF-8 sort); then
        echo "Warning: you want to build for testing but the cfiles isnt"
        echo "updated with the needed files. Read the readme."
        BAD_FILE_LISTS=1
    elif ! diff cmakeFileLists/cFileList.txt <(find "./src" -name "*.cc" -o -name "*.cpp" |  LC_ALL=C.UTF-8 sort); then
        echo "Warning: you want to build but the cfiles isnt updated with"
        echo "the needed files. Read the readme."
        BAD_FILE_LISTS=1
    fi
    #check our hdirectory stuff
    if ! diff cmakeFileLists/hDirectoryList.txt <( find "./includes" -type d ! -path "*External_Libraries*" |  LC_ALL=C.UTF-8 sort) ; then
        echo "Warning: you want to build for testing but the hdirectory list file"
        echo "isnt updated with the needed files. Read the readme."
        BAD_FILE_LISTS=1
    fi
    #check our external lib stuff
    if ! diff cmakeFileLists/externalHDirectoryList.txt <(find "./includes" -type d -path "*External_Libraries*" |  LC_ALL=C.UTF-8 sort) ; then
        echo "Warning: you want to build for testing but the list of directories"
        echo "for our external libs isnt updated with the needed files. Read" 
        echo "the readme."
        BAD_FILE_LISTS=1
    fi
    
    if [ "${BAD_FILE_LISTS}" == 1 ]; then
        echo "==============================================================="
        echo "READ THE README, IT WILL EXPLAIN WHAT IS GOING ON. CHECK UNDER"
        echo "THE COMPILING THE LIBRARY SECTION. IF PEOPLE CANNOT FOLLOW THIS"
        echo "WORKAROUND I WILL IMPLEMENT CMAKE IN A MORE TRADIONAL WAY WHICH"
        echo "WILL MAKE DEVELOPMENT MORE ANNOYING"
        echo "==============================================================="
        exit "${BAD_FILE_LISTS}"
    fi
    echo "File lists are good lets run it"
   
}

################################################################
##########               Print help message         ############
################################################################

printHelp()
{
	echo "*************************************************************"
	echo "If you are a user then this $0 is not for you:"
	echo "go to the GEMS home dir; If needed GEMS $0 will spawn GMML building."
	echo ""
	echo "This is for building GMML in isolation! You can wrap GMML here"
	echo "using swig but interfacing with the wrapped lib will be up to you "
	echo ""
	echo "GEMSHOME should be set to the parent of the gmml directory."
	echo "After isolated GMML is built the next step is testing; do this:"
	echo "cd tests; compile_run_tests.bash"
	printf "*************************************************************\n"
	printf "Please note that once GMML is built, you can test it by running:\n"
	printf "cd tests; ./compile_run_tests.bash\n"
	printf "*************************************************************\n"
	printf "Options are as follows:\n"
	printf "\t-c\t\t\tClean all files from previous builds\n"
	printf "\t-j <NUM_JOBS>\t\tBuild GMML with <NUM_JOBS>\n"
	printf "\t-o <O0/O2/OG/debug>\tBuild GMML using no optimization, 2nd \n\t\t\t\tlevel optimization, or with debug symbols\n"
	printf "\t-w\t\t\tWrap gmml in python using swig\n"
	printf "\t-h\t\t\tPrint this help message and exit\n"
	printf "*************************************************************\n"
	echo "Exiting."
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
TARGET_MAKE_FILE="Makefile"
CMAKE_BUILD_TYPE_FLAG="-DCMAKE_BUILD_TYPE=Release"
#blank means all
MAKE_TARGET="all"

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

while getopts "j:o:cwh" option
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
				CLEAN="1"
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

echo "Starting installation of GMML at $(date)".

#check our file lists before we do anything 
checkCMakeFileLists

printf "\nBuilding with these settings:\n"
echo "GEMSHOME: ${GEMSHOME}"
echo "CMake build type: ${CMAKE_BUILD_TYPE_FLAG}"
echo "Build Target: ${MAKE_TARGET}"
echo "Clean: ${CLEAN}"

################################################################
#########                  COMPILE GMML                #########
################################################################

echo "Generating GMML ${TARGET_MAKE_FILE}."

if [ "${CLEAN}" == "1" ]; then
	echo ""
	echo "Cleaning GMML. Please note this is VERY aggressive!"
	if [ -d "./cmakeBuild" ]; then
		(cd ./cmakeBuild; make -f "${TARGET_MAKE_FILE}" clean)
		rm -rf ./cmakeBuild
	fi
	if [ -d "./lib" ]; then
		rm -rf ./lib
	fi
fi

#Note that we have to generate our makefile with cmake before we build
# The -S (source) flag is for where the cmakelists file is located
# and the -B (build dir) flag is for where everything will be built to
echo "Creating makefile using cmake"
echo ""
cmake "${CMAKE_BUILD_TYPE_FLAG}" -S . -B ./cmakeBuild || { echo "ERROR RUNNING CMAKE ON GMML: $0 FAILED, EXITING" ; exit 1; }


#NOTE: All our build stuff is within the cmakeBuild dir
echo "Making GMML."
cd cmakeBuild || { echo "ERROR BUILDING GMML: $0 FAILED, EXITING" ; exit 1; }
make -j"${NMP}" "${MAKE_TARGET}" || { echo "ERROR BUILDING GMML: $0 FAILED, EXITING" ; exit 1; }
cd ..

################################################################
#########              WRAP UP TO GEMS                 #########
################################################################

STUPID_WRAPPED_MSG=""
#Wrapping is handled by cmake
if [ "${MAKE_TARGET}" == "gmml_wrapped" ]; then
    STUPID_WRAPPED_MSG=" and wrapping"
fi
echo ""
echo "GMML compilation${STUPID_WRAPPED_MSG} is finished at $(date)".
exit

