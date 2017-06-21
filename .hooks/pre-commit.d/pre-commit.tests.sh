#!/bin/bash

check_gemshome() {
   if [ -z "$GEMSHOME" ]; then
      echo ""
      echo "Your GEMSHOME environment variable is not set! It should be set to"
      echo "$1"
      exit 1
   elif [ ! -d $GEMSHOME ]; then
      echo ""
      echo "Your GEMSHOME environment variable is set to $GEMSHOME -- this does"
      echo "not appear to be a directory. It should be set to"
      echo "$1"
      exit 1
   elif [ ! "$GEMSHOME" = "$1" -a ! "$GEMSHOME" = "${1}/" ]; then
      #try checking the inode incase there is a problem with symlinks
       if [ `stat -c "%i" $GEMSHOME` != `stat -c "%i" ${1}` ]; then
           echo ""
           echo "ERROR: GEMSHOME is expected to be $1 but it is currently"
           echo "$GEMSHOME    This will cause problems!"
           exit 1
       fi
   fi
}

cd ../
gemshome=`pwd`
cd -
check_gemshome $gemshome 

echo "Running mandatory tests..."
echo "Q1. What... is the air-speed velocity of an unladen swallow?"
cd $GEMSHOME/gmml/tests/
 source compile_run_tests.bash
cd -

if [ -f $GEMSHOME/gmml/tests/All_Tests_Passed ] ; then
    echo  "All tests have passed. Commits are allowed."
    exit 0
else
    echo "
         *****************************************************************
         You're a naughty boy, Fawlty!
         Tests have not been run, or have failed! 
         Commits are not allowed, and the police are on the way.
         Run tests/compile_run_tests.bash before commiting changes. 
         *****************************************************************
         "
    exit 1
fi
