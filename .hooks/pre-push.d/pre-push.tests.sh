#!/bin/bash

check_gemshome()
{
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
        if [ $(stat -c "%i" $GEMSHOME) != $(stat -c "%i" ${1}) ]; then
            echo ""
            echo "ERROR: GEMSHOME is expected to be $1 but it is currently"
            echo "$GEMSHOME    This will cause problems!"
            exit 1
        fi
    fi
}

check_dir_exists()
{
    if [ ! -d "$1" ]; then
        echo ""
        echo "Your $1 directory does not exist."
    fi
}

cd ../
gemshome=$(pwd)
cd -
check_gemshome $gemshome

## OG Oct 2021 have the hooks update themselves.
cp -r $GEMSHOME/gmml/.hooks/* $GEMSHOME/gmml/.git/hooks/

TEST_SKIP=0
#### Allow skipping tests ####
branch=$(git rev-parse --abbrev-ref HEAD)
if [[ "$branch" != "gmml-dev" ]] && [[ "$branch" != "gmml-test" ]] && [[ "$branch" != "stable" ]]; then
    printf "Branch is %s\nSkipping tests is allowed.\nDo you want to skip them?\ns=skip\na=abort\nEnter anything to run tests.\n" $branch
    read -p "Enter response: " response </dev/tty
    if [[ $response == [sS] ]]; then
        echo -e "Skipping tests!\n"
        TEST_SKIP=1
    elif [[ $response == [aA] ]]; then
        printf "Abort!\n"
        exit 1
    else
        printf "Running tests.\n"
    fi
fi

#sane git checking
echo "Checking if our current branch is on remote, if the branch status below is empty"
echo -e "then we know that the branch is not on remote.\n"
echo "Branch hash and name: $(git ls-remote --heads origin "$(git rev-parse --abbrev-ref HEAD)")"
echo ""
if [ -n "$(git ls-remote --heads origin "$(git rev-parse --abbrev-ref HEAD)")" ]; then
    #we hit here if our branch is actually on remote, thus we must check
    #that the current branch is up to date on remote
    echo "Branch is on remote, now to check if local is behind remote"
    if [ "$(git rev-list --left-only --count origin/"$(git rev-parse --abbrev-ref HEAD)"..."$(git rev-parse --abbrev-ref HEAD)")" != 0 ]; then
        #here the given value is non zero thus remote is ahead of our local so we want to go ahead and stop everything and just exit with an error
        echo "ERROR: REMOTE IS AHEAD OF YOUR LOCAL BRANCH, PULL BEFORE YOU TRY TO PUSH"
        #since remote is ahead of local we know we want to stop the push and make the user pull the new code
        exit 1
    fi
    echo "Local branch is not behind the remote branch, proceeding"
else
    echo "Branch is not on remote, so no need to check if local is behind remote, proceeding"
fi
#if we dont wanna do tests we just exit here cause we know we are kosher to push
if [ "${TEST_SKIP}" == 1 ]; then
    echo "Skipping tests, you are good to push."
    exit 0
else
    echo "Beginning tests"
fi

#we only hit here if we know that we arent skipping tests and the local branch is not behind the remote branch (if the
#remote branch even exists so lets run it

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
if [ $result -eq 0 ]; then
    echo "GMML level tests have passed. Doing gems level tests."
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
