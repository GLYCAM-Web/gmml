#!/bin/bash
 
if [ -f ../../tests/All_Tests_Passed ] ; then
    echo  "All tests have passed. Commits are allowed."
    exit 1
else
    echo "
         *****************************************************************
         You're a naughty boy, Fawlty!
         Tests have not been run, or have failed! 
         Commits are not allowed, and the police are on the way.
         Run tests/compile_run_tests.bash before commiting changes. 
         *****************************************************************
         "
    exit 0
fi
