#!/bin/bash

if [ -z "$GEMSHOME" ]; then
      echo ""
      echo "Your GEMSHOME environment variable is not set! It should be set something like
            export GEMSHOME=/yourpath/gems/"
      exit 1
fi
 
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
