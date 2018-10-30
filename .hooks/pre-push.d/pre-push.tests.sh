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

#Compile gmml if not compiled:
echo "Pulling all changes"
git pull
echo "Compiling gmml with ./make.sh no_clean no_wrap"
cd $GEMSHOME/
 git pull
 #Add these removes so the tests don't pass on an old version of the library
 rm -f ./gmml/bin/libgmml.so.1.0.0
 rm -f ./gmml/bin/libgmml.so
 rm -f ./gmml/bin/libgmml.so.1
 rm -f ./gmml/bin/libgmml.so.1.0
 rm -rf gmml_wrap.cxx gmml_wrap.o gmml.py gmml.pyc _gmml.so
 ./make.sh no_clean wrap
cd -

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
         echo "GEMS level tests have failed. Make sure you have pulled the latest version. At time of writing (July 2018) you probably need the gems-dev branch if using the gmml-dev branch." 
         echo "If you are up-to-date, this failure indicates that you have caused the outputs of $GEMSHOME/tests to change. You can open the $GEMSHOME/tests/run_tests.sh file and run the test line by line to get an output file. Compare it to the saved \"correct\" version in $GEMSHOME/tests/correct_outputs." 
         echo "Sometimes the changes you make are fine, and you just need to update what the correct output is by overwriting the old output. Make sure it is ok though, or you will be mur-didely-urdered."
         exit 1
     else
         echo "GEMS level tests have passed. Checking glycoprotein builder."
         if [ ! -d "$GEMSHOME/gmml/programs/GlycoproteinBuilder" ]; then
             cd $GEMSHOME/gmml/programs/
             git clone https://github.com/gitoliver/GlycoProteinBuilder.git GlycoproteinBuilder
         fi
         cd $GEMSHOME/gmml/programs/GlycoproteinBuilder
         git pull
         make clean
         make
         ./run_tests.sh
         gpbuilder_tests_result=$? # record the exit status of previous command
         if [ $gpbuilder_tests_result -eq 0 ]; then
             echo "All tests have passed. Pushing allowed"
             exit 0
         else
             echo "The tests in $GEMSHOME/gmml/programs/GlycoproteinBuilder have failed. Check $GEMSHOME/gmml/programs/GlycoproteinBuilder/run_tests.sh" 
             echo "Push cancelled."
         fi
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
