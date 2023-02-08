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

check_dir_exists() {
    if [ ! -d "$1" ]; then
        echo ""
        echo "Your $1 directory does not exist."
    fi
}

cd ../
 gemshome=`pwd`
cd -
check_gemshome $gemshome 

## OG Oct 2021 have the hooks update themselves.
cp -r $GEMSHOME/gmml/.hooks/* $GEMSHOME/gmml/.git/hooks/

#### Allow skipping tests ####
branch=`git rev-parse --abbrev-ref HEAD`
if [[ "$branch" != "gmml-dev" ]] && [[ "$branch" != "gmml-test" ]] && [[ "$branch" != "stable" ]]; then
    printf "Branch is %s\nSkipping tests is allowed.\nDo you want to skip them?\ns=skip\na=abort\nEnter anything to run tests.\n" $branch
    read -p "Enter response: " response < /dev/tty
    if [[ $response == [sS] ]]; then
        printf "Skipping tests!\n"
        exit 0;
    elif [[ $response == [aA] ]]; then
        printf "Abort!\n"
        exit 1;
    else
        printf "Running tests.\n"
    fi
fi
#### Allow skipping tests ####



#Compile gmml if not compiled:
echo "Pulling all changes"
git pull
result=$? # record the exit status of previous command
if [ $result -eq 1 ] ; then
    echo "Could not pull gmml"
    exit 1
fi
cd $GEMSHOME/
 git pull
 result=$? # record the exit status of previous command
 if [ $result -eq 1 ] ; then
     echo "Could not pull gems"
     exit 1
 fi
 #Add these removes so the tests don't pass on an old version of the library
 rm -f gmml.py _gmml.so
 rm -rf ./gmml/lib
if [ -d "./gmml/cmakeBuild" ]; then
	 echo "Removing the libgmml.so from our cmakeBuild directory"
	rm ./gmml/cmakeBuild/libgmml.so
fi
echo "Compiling gmml using GEMS ./make.sh, no wrap flag cause it auto wraps"
 ./make.sh 
 cd $GEMSHOME/gmml

echo "Running mandatory tests..."
cd $GEMSHOME/gmml/tests/
 bash compile_run_tests.bash
 result=$? # record the exit status from compile_run_tests.bash
cd -
if [ $result -eq 0 ] ; then
    echo  "GMML level tests have passed. Doing gems level tests."
    cd $GEMSHOME/tests/
     bash run_tests.sh
     gems_tests_result=$? # record the exit status of previous command
     if [ $gems_tests_result -ne 0 ]; then
         echo "GEMS level tests have failed. Make sure you have pulled the latest version and are on the appropriate branch. " 
         echo "If you are up-to-date, this failure indicates that you have caused the outputs of $GEMSHOME/tests to change. You can open the $GEMSHOME/tests/run_tests.sh file and run the test line by line to get an output file. Compare it to the saved \"correct\" version in $GEMSHOME/tests/correct_outputs." 
         echo "Sometimes the changes you make are fine, and you just need to update what the correct output is by overwriting the old output. Make sure it is ok though, or you will be mur-didely-urdered."
         exit 1
     else
         echo "All tests have passed. Pushing allowed"
         exit 0
     fi
else
    echo "
         *****************************************************************
         The GMML level tests have failed! 
         Push cancelled.
         *****************************************************************
         "
    exit 1
fi
