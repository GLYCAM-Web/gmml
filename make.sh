#!/bin/bash

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
	elif [ ! -d ${GEMSHOME} ]; then
		echo ""
		echo "Error:  GEMSHOME environment variable is set to ${GEMSHOME} -- this does"
		echo "not appear to be a directory. It should be set to"
		echo "$1"
		exit 1
	fi
}

################################################################
##########               Print help message         ############
################################################################

printHelp()
{
	echo ""
	echo "If you are a user then this $0 is not for you:"
	echo "go to the GEMS home dir; If needed GEMS $0 will spawn GMML building."
	echo ""
	echo "This is for building GMML in isolation! You can wrap GMML here"
	echo "using swig but interfacing with the wrapped lib will be up to you "
	echo ""
	echo "GEMSHOME should be set to the parent of the gmml directory."
	echo "After isolated GMML is built the next step is testing; do this:"
	echo "cd tests; compile_run_tests.bash"
	echo ""

	printf "*************************************************************\n"
	printf "Please note that once GMML is built, you can test it by running:\n"
	printf "cd tests; ./compile_run_tests.bash\n"
	printf "*************************************************************\n"
	printf "Options are as follows:\n"
	printf "\t-c\t\t\tClean all files from previous builds\n"
	printf "\t-j <NUM_JOBS>\t\tBuild GMML with <NUM_JOBS>\n"
	printf "\t-o <O0/O2/OG/debug>\tBuild GMML using no optimization, 2nd \n\t\t\t\tlevel optimization, or with debug symbols\n"
	printf "\t-w\t\t\tWrap gmml in python using swig\n"
	printf "*************************************************************\n"
	echo "Exiting."
	exit 1
}

################################################################
#########                CHECK SETTINGS                #########
################################################################

echo "Starting installation of GMML at `date`".

gemshome=`pwd`
check_gemshome ${gemshome}
get_numprocs

################################################################
#########              CREATE CLIENT HOOKS             #########
################################################################

#Cannot use server side hooks when hosting on git-hub.
#Stuff in .git/hooks is ignored by git.
#Solution: The folder .hooks is tracked by git.
# Copy .hooks to .git/hooks during installation.
cp -r $GEMSHOME/gmml/.hooks/* $GEMSHOME/gmml/.git/hooks/
#I don't think this is ideal, and is perhaps silly. OG Apr 2017.

################################################################
#########                SET UP DEFAULTS               #########
################################################################

#Changing the following options are for developers.
CLEAN="0"
TARGET_MAKE_FILE="Makefile"
CMAKE_BUILD_TYPE_FLAG="-DCMAKE_BUILD_TYPE=Release"
WRAP="0"

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

while getopts "j:o:cw" option
do
	case "$option" in
			j)
				jIn=${OPTARG}
				if [[ ${jIn} =~ ^[0-9]+$ ]]; then
					NMP=${jIn}
				else
					printHelp
				fi
				;;
			o)
				oIn=${OPTARG}
				if [ ${oIn} == "O0" ] || [ ${oIn} == "no_optimize" ]; then
					CMAKE_BUILD_TYPE_FLAG="-DCMAKE_BUILD_TYPE=RelWithDebInfo"
				elif [ ${oIn} == "O2" ] || [ ${oIn} == "optimize" ]; then
					CMAKE_BUILD_TYPE_FLAG="-DCMAKE_BUILD_TYPE=Release"
				elif [ ${oIn} == "debug" ] || [ ${oIn} == "OG" ]; then
					CMAKE_BUILD_TYPE_FLAG="-DCMAKE_BUILD_TYPE=Debug"
				else
					printHelp
				fi
				;;
			c)
				CLEAN="1"
				;;
			w)
				WRAP="1"
				;;
			*)
				printHelp
				;;
	esac
done

printf "\nBuilding with these settings:\n"
printf "GEMSHOME: $GEMSHOME\n"
printf "TARGET_MAKE_FILE: $TARGET_MAKE_FILE\n"
printf "CMake build type: $CMAKE_BUILD_TYPE_FLAG\n"
printf "Wrap: $WRAP\n"
printf "Clean: $CLEAN\n"

#NOTE: DEFUNCT NEED TO REMOVE, HERE FOR REFERENCE ABOUT WHAT OLD FLAGS WERE
#i=1
#while [ ${i} -le $# ]; do
#	argument="${!i}"
#	if [ "$argument" = "clean" ]||[ "$argument" = "no_clean" ];then
#		CLEAN="${!i}"
#	elif [ "$argument" = "wrap" ]||[ "$argument" = "no_wrap" ];then
#		WRAP_GMML="${!i}"
#	elif [ "$argument" = "debug" ]||[ "$argument" = "no_debug" ];then
#		DEBUG="${!i}"
#	elif [ "$argument" = "optimize" ]||[ "$argument" = "no_optimize" ]||[ "$argument" = "O1" ]||[ "$argument" = "O2" ];then
#		OPTIMIZE="${!i}"
#	fi
#	i=$[$i+1]
#done

################################################################
#########                  COMPILE GMML                #########
################################################################

echo "Generating GMML $TARGET_MAKE_FILE."

if [ "$CLEAN" == "1" ]; then
	echo ""
	echo "Cleaning GMML. Please note this is VERY aggressive!"
	if [ -d "./cmakeBuild" ]; then
		cd ./cmakeBuild
		make -f $TARGET_MAKE_FILE clean
		cd ..
		rm -rf ./cmakeBuild
	fi
	if [ -d "./lib" ]; then
		rm -rf ./lib
	fi
#NOTE: The flattened build dir is soon to be removed
	if [ -d "./build" ]; then
		rm -rf ./build
	fi
fi

#Note that we have to generate our makefile with cmake before we build
# The -S (source) flag is for where the cmakelists file is located
# and the -B (build dir) flag is for where everything will be built to
echo "Creating makefile using cmake"
echo ""
cmake ${CMAKE_BUILD_TYPE_FLAG} -S . -B ./cmakeBuild -DWRAP_GMML=${WRAP}


#NOTE: All our build stuff is within the cmakeBuild dir
echo "Making GMML."
cd cmakeBuild
make -j${NMP}
cd ..
if [ -f ./cmakeBuild/libgmml.so ]; then
#	mkdir lib
#	cp ./cmakeBuild/libgmml.so lib/
else
	echo "Warning: libgmml was not created!"
fi


################################################################
#########              WRAP UP TO GEMS                 #########
################################################################

#Wrapping is handled by cmake

echo ""
echo "GMML compilation is finished at `date`".
exit

