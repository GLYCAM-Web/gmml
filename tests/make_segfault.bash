#!/bin/bash

#Manually change this number as you add tests:

# Required for compiling
# Can do this, or with root copy bin/* to /usr/lib/
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:../bin
export LD_LIBRARY_PATH

echo "Testing the segfault generator... "

COMMAND1='g++ -std=c++0x -I $GEMSHOME/gmml/includes/ -L$GEMSHOME/gmml/bin/ -Wl,-rpath,$GEMSHOME/gmml/bin/ tests/isegfault.cc -lgmml -o isegfault '
COMMAND2='./isegfault > /dev/null  2>&1'

echo "Compiling now with:"
echo ${COMMAND1}
eval ${COMMAND1}

echo "Running now with:"
echo ${COMMAND2}
##  The following will run the segfault text into a file.
##  { (eval ${COMMAND2})   2>&1 & } >segoutput.txt
##  If that is desired, use that instead of the following.
##  Note that changing from the following will frustrate 
##  attempts to capture the segfault exit code.
eval ${COMMAND2}
exitcode=$?
if [ "${exitcode}" -ne "139" ] && [ "${exitcode}" -ne "-11" ] ; then
	echo "The exit code is ${exitcode}, so there was no segfault.  Test failed."
else
	echo "Test Passed."
fi
rm isegfault


