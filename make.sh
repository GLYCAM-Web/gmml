#!/bin/bash

################################################################
#########                 FUNCTIONS                    #########
################################################################

get_numprocs() {
    if [ -z "$GEMSMAKEPROCS" ]; then
        NMP=4
    elif ! [[ $GEMSMAKEPROCS =~ ^[0-9]+$ ]] ; then
        echo "Warning:  GEMSMAKEPROCS is not a valid integer; setting to 4"
        NMP=4
    elif [ "$GEMSMAKEPROCS" -eq "0" ] ; then
        echo "Warning:  GEMSMAKEPROCS cannot be zero; setting to 4"
        NMP=4
    else
        NMP=$GEMSMAKEPROCS
    fi
}

check_gemshome() {
    if [ -z "$GEMSHOME" ]; then
        echo ""
        echo "Error:  GEMSHOME environment variable is not set! It should be set to"
        echo "$1"
        exit 1
    elif [ ! -d $GEMSHOME ]; then
        echo ""
        echo "Error:  GEMSHOME environment variable is set to $GEMSHOME -- this does"
        echo "not appear to be a directory. It should be set to"
        echo "$1"
        exit 1
# skip this since we are using GEMSHOME as a GMMLHOME.
#    elif [ ! "$GEMSHOME" = "$1" -a ! "$GEMSHOME" = "${1}/" ]; then
#        #try checking the inode incase there is a problem with symlinks
#        if [ `stat -c "%i" $GEMSHOME` != `stat -c "%i" ${1}` ]; then
#            echo ""
#            echo "Error:  GEMSHOME is expected to be $1 but it is currently"
#            echo "$GEMSHOME"
#            echo "        This will cause problems!"
#            exit 1
#        fi
    fi
}

################################################################
##########               Print help message         ############
################################################################

if [[ "$1" == "-help" ]] || [[ "$1" == "-h" ]] || [[ "$1" == "--help" ]]; then
    echo ""
    echo "If you are a user then this $0 is not for you:"
    echo "    go to the GEMS home dir; GEMS $0 will spawn GMML building."
    echo ""
    echo "This is for building GMML in isolation!"
    echo "    There is no wrapping and no connection to GEMS "
    echo "    although environment variables starting with GEMS are respected;"
    echo "    GEMSHOME should be set to the parent of the gmml directory."
    echo "After isolated GMML is built the next step is testing; do this:"
    echo "cd tests; compile_run_tests.bash"
    echo ""

    printf "*************************************************************\n"
    printf "Usage: $0 clean_gmml? debug_gmml?\n"
    printf "Example: $0 clean no_debug\n"
    printf "Default: $0 no_clean debug\n"
    printf "*************************************************************\n"
    printf "If selected the options do this:\n"
    printf "     1. Cleans gmml before making\n"
    printf "     2. Adds compiler options for debugging\n"
    printf "*************************************************************\n"
    echo "Exiting."
    exit 1
fi

################################################################
#########                CHECK SETTINGS                #########
################################################################

echo "Starting installation of GMML at `date`".

gemshome=`pwd`
check_gemshome $gemshome
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
CLEAN="no_clean"
DEBUG="no_debug"
TARGET_MAKE_FILE="Makefile"
OPTIMIZE="optimize"

################################################################
#########               COMMAND LINE INPUTS            #########
################################################################
i=1
while [ ${i} -le $# ]; do
    argument="${!i}"
    if [ "$argument" = "clean" ]||[ "$argument" = "no_clean" ];then
        CLEAN="${!i}"
    elif [ "$argument" = "wrap" ]||[ "$argument" = "no_wrap" ];then
        WRAP_GMML="${!i}"
    elif [ "$argument" = "debug" ]||[ "$argument" = "no_debug" ];then
        DEBUG="${!i}"
    elif [ "$argument" = "optimize" ]||[ "$argument" = "no_optimize" ]||[ "$argument" = "O1" ]||[ "$argument" = "O2" ];then
        OPTIMIZE="${!i}"
    fi
    i=$[$i+1]
done

printf "\nBuilding with these settings:\n"
printf "GEMSHOME: $GEMSHOME\n"
printf "TARGET_MAKE_FILE: $TARGET_MAKE_FILE\n"
printf "CLEAN: $CLEAN\n"
printf "DEBUG: $DEBUG\n"
printf "OPTIMIZE: $OPTIMIZE\n\n"

################################################################
#########                  COMPILE GMML                #########
################################################################

 if [ "$DEBUG" == "debug" ]; then
     DEBUGOPTIONS='-g'
 else
     DEBUGOPTIONS=''
 fi
 
 if [ "$OPTIMIZE" == "optimize" ]; then
     OPTIMIZE=""
     #Does -O2 by default
     NO_OPTIMIZE=""
 elif [ "$OPTIMIZE" == "no_optimize" ]; then
     OPTIMIZE=""
     NO_OPTIMIZE="QMAKE_CXXFLAGS_RELEASE -= -O2 -O1"
 elif [ "$OPTIMIZE" == "O1" ]; then
     OPTIMIZE="-O1"
     NO_OPTIMIZE="QMAKE_CXXFLAGS_RELEASE -= -O2"
 elif [ "$OPTIMIZE" == "O2" ]; then
     OPTIMIZE="-O2"
     NO_OPTIMIZE=""
 else
     OPTIMIZE=""
     NO_OPTIMIZE=""
 fi

 echo "Generating GMML $TARGET_MAKE_FILE."
 # Always create a new gmml.pro and makefile
 ## This is going to be broken up to variables instead of being this long command. Just wanted to get a working version pushed up.
 qmake -project -t lib -o gmml.pro "QMAKE_CXXFLAGS += -Wall -W -std=c++11 ${DEBUGOPTIONS} ${OPTIMIZE}" "QMAKE_CFLAGS += -Wall -W ${DEBUGOPTIONS}" "${NO_OPTIMIZE}" "DEFINES += _REENTRANT" "CONFIG = no_lflag_merge" "unix:LIBS = -L/usr/lib/x86_64-linux-gnu -lpthread" "OBJECTS_DIR = build" "DESTDIR = lib" -r src/ includes/ -nopwd
 qmake -o $TARGET_MAKE_FILE

 if [ "$CLEAN" == "clean" ]; then
     echo ""
     echo "Cleaning GMML."
     make -f $TARGET_MAKE_FILE distclean
     qmake -o $TARGET_MAKE_FILE
 fi

 echo "Making GMML."
  make -j ${NMP} -f $TARGET_MAKE_FILE 


################################################################
#########              WRAP UP TO GEMS                 #########
################################################################

# Wrapping is performed in GEMS make.sh

echo ""
echo "GMML compilation is finished at `date`".
exit

