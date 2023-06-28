#!/bin/bash

##### PRETTY PRINT STUFF #####
#lazy and dont want to have to type all these color variables a bunch, also improves readability
RESET_STYLE='\033[0m'
GREEN_BOLD='\033[0;32m\033[1m'
YELLOW_BOLD='\033[0;33m\033[1m'
RED_BOLD='\033[0;31m\033[1m'

#To replace our gemshome dealio. Not great but egh should be good enough...... hopefully.....
if [[ $(pwd) == */gems/gmml ]]; then
    GMML_DIR=$(pwd)
    GEMS_DIR=$(cd .. && pwd)
else
    echo -e "${RED_BOLD}ERROR: Trying to run from incorrect directory being: $(pwd)\nRun from the base GMML directory.\nABORTING${RESET_STYLE}"
    exit 1
fi

#Just a global var to check if GEMS or GMML is out of date. Here we are using string existence
#instead of a wannabe bool type flag
BRANCH_IS_BEHIND=""

#How many commits a feature branch can be missing from the gmml-test branch
MAX_FEATURE_BEHIND_TEST=15

#this isnt coded for babies, follow our git patterns and we are good, BORK ON MERGE????
ensure_feature_close()
{
    CURRENT_BRANCH=$(git rev-parse --abbrev-ref HEAD)
    echo -e "\nEnsuring feature branch ${YELLOW_BOLD}${CURRENT_BRANCH}${RESET_STYLE} in repo ${YELLOW_BOLD}GMML${RESET_STYLE} is\nless than ${MAX_FEATURE_BEHIND_TEST} commits behind gmml-test branch"

    if [ "$(git rev-list --left-only --count origin/gmml-test..."${CURRENT_BRANCH}")" -gt "${MAX_FEATURE_BEHIND_TEST}" ]; then
        echo -e "${RED_BOLD}ERROR:${RESET_STYLE} MISSING TOO MANY COMMITS FROM GMML-TEST, INCORPERATE CURRENT GMML-TEST CODE INTO YOUR FEATURE BRANCH THEN TRY AGAIN.\nABORTING\n"
        exit 1
    else
        echo -e "${GREEN_BOLD}passed...${RESET_STYLE} Feature branch is not too far from gmml-test, proceding\n"
    fi
}

if [ -z "${GEMSHOME}" ]; then
    echo -e "${YELLOW_BOLD}WARNING: Your GEMSHOME environment variable is not set! It should be set to the GEMS directory\nthat is the parent of the GMML directory. This can cause some issues as some of the codebase\nstill relies upon the GEMSHOME variable. Continuing but you have been warned.${RESET_STYLE}"
fi

## OG Oct 2021 have the hooks update themselves.
#TODO: Do this more auto like, if this script is updated the first run the next time will not
#reflect the made changes due to the old script calling the copy then continuing.
cp -r "${GMML_DIR}"/.hooks/* "${GMML_DIR}"/.git/hooks/

#imagine this means "check if current branch is behind origin of the same branch". Basically all
#we are doing are checking either the GEMS or GMML repo to ensure there are no commits that the
#user has not pulled.
# @param $1 is either GEMS or GMML, whichever defined is what will be checked and pathing
# will be based off of the gemsdir variable. I picked caps in order to just make the check more
# apparent about how its checking an actual repo status
check_if_branch_behind()
{
    #we check we have a single arg, then if the single arg is what we want.
    if [ "$#" == 1 ] && {
        [ "$1" == "GEMS" ] || [ "$1" == "GMML" ]
    }; then
        #TODO: More legit error checking in this
        if [ "$1" == "GEMS" ]; then
            cd "${GEMS_DIR}" || {
                echo -e "${RED_BOLD}failed...${RESET_STYLE} We could not change directory to the following:\n\t ${GEMS_DIR}"
                echo "Exiting..."
                exit 1
            }
        elif [ "$1" == "GMML" ]; then
            cd "${GEMS_DIR}/gmml" || {
                echo -e "${RED_BOLD}failed...${RESET_STYLE} We could not change directory to the following:\n\t ${GEMS_DIR}/gmml"
                echo "Exiting..."
                exit 1
            }
        else
            echo -e "${RED_BOLD}Error checking if branch was behind when changing dir!${RESET_STYLE}"
            echo "Given param: $1"
            exit 1
        fi
        CURRENT_BRANCH=$(git rev-parse --abbrev-ref HEAD)
        #let the user know what is going on, i.e. what branch is being checked and in what repo
        echo -e "\nEnsuring branch ${YELLOW_BOLD}${CURRENT_BRANCH}${RESET_STYLE} in repo ${YELLOW_BOLD}$1${RESET_STYLE} is up to date...."

        #Done to check if we actually need to hit up origin to see if we are behind, the "branch hash and name on origin"
        #part is mostly there to bug check
        echo -e "First checking if ${YELLOW_BOLD}${CURRENT_BRANCH}${RESET_STYLE} is on remote"
        #if the branch status below is empty then we know that the branch is not on remote
        #this was just for deboogin, keep it cause it can help - P
        #echo "Branch hash and name on origin: $(git ls-remote --heads origin "${CURRENT_BRANCH}")"
        #check if we get a non-empty return aka branch is on remote, if we are on a head branch we want to "force" acceptance to go ahead
        #and add if our branch name is HEAD then we fail
        if [ -n "$(git ls-remote --heads origin "${CURRENT_BRANCH}")" ] || [ "${CURRENT_BRANCH}" == "HEAD" ]; then
            #we hit here if our branch is actually on remote, thus we must check
            #that the current branch is up to date on remote
            echo -e "Branch is present on remote, now to check if local is behind remote"
            #the left only part of this command shows how many commits behind the repo on the right is from
            #the left repo that has origin tacked onto it., the right side
            if [ "$(git rev-list --left-only --count origin/"${CURRENT_BRANCH}"..."${CURRENT_BRANCH}")" != 0 ] || [ "${CURRENT_BRANCH}" == "HEAD" ]; then
                #here the given value is non zero thus remote is ahead of our local so we want to go ahead and stop everything and just exit with an error
                echo -e "${RED_BOLD}ERROR:${RESET_STYLE} $1 REMOTE IS AHEAD OF YOUR LOCAL BRANCH, PULL BEFORE YOU TRY TO PUSH"
                #since remote is ahead of local we know we want to stop the push and make the user pull the new code,
                #BUT we want to go ahead and check both repos first and THEN exit so the user can know if they
                #need to pull both GEMS and GMML. Basically this var is only populated if local branch is behind origin branch
                if [ -n "${BRANCH_IS_BEHIND}" ]; then
                    BRANCH_IS_BEHIND="${BRANCH_IS_BEHIND}\n\t$1"
                else
                    BRANCH_IS_BEHIND="\t$1"
                fi
            else
                echo -e "${GREEN_BOLD}passed...${RESET_STYLE}Local branch ${YELLOW_BOLD}${CURRENT_BRANCH}${RESET_STYLE} on repo ${YELLOW_BOLD}$1${RESET_STYLE} is not behind the remote branch, proceeding"
            fi
        else
            echo -e "${GREEN_BOLD}passed...${RESET_STYLE} Branch is not on remote, so no need to check if local is behind remote, proceeding"
        fi
        cd "${GEMS_DIR}"/gmml || {
            echo -e "${RED_BOLD}failed...${RESET_STYLE} We could not change directory to the following:\n\t ${GEMS_DIR}/gmml"
            echo "Exiting..."
            exit 1
        }
    #bad input, so end script
    else
        echo -e "${RED_BOLD}ERROR${RESET_STYLE} Incorrect param given to check_if_branch_behind function. Exiting."
        echo -e "\tInput given:$1"
        exit 1
    fi
} #End ensuring branches not behind function

#before we try anything major we first figure out if any branches are behind.
check_if_branch_behind "GMML"
check_if_branch_behind "GEMS"
#echo -e "Checking if branches are behind remote completed\n"
#if they are behind then the string we are checking wont be empty, thus we should exit.
if [ -n "${BRANCH_IS_BEHIND}" ]; then
    echo -e "${RED_BOLD}failed...${RESET_STYLE} At least one of you branches are behind."
    echo -e "Pull the following repos:\n${RED_BOLD}${BRANCH_IS_BEHIND}${RESET_STYLE}\n"
    echo "Exiting..."
    exit 1
fi

TEST_SKIP=0
#### Allow skipping tests ####
branch=$(git rev-parse --abbrev-ref HEAD)
if [[ "${branch}" != "gmml-dev" ]] && [[ "${branch}" != "gmml-test" ]] && [[ "${branch}" != "stable" ]]; then

    if [[ "${branch}" != hotfix* ]]; then
        #Ensuring feature branch is not too far from gmml-test
        ensure_feature_close
    else
        echo -e "You are applying making a hotfix for one of our main branches EXCEPT gmml-test,\nif this is not the case ABORT AND BUG PRESTON\n"
    fi

    echo -e "Branch is ${branch}\nSkipping tests is allowed.\nDo you want to skip them?\ns=skip\na=abort\nEnter anything to run tests.\n"
    read -p "Enter response: " response </dev/tty
    if [[ "${response}" == [sS] ]]; then
        echo -e "Skipping tests!\n"
        TEST_SKIP=1
    elif [[ "${response}" == [aA] ]]; then
        printf "Abort!\n"
        exit 1
    else
        printf "Running tests.\n"
    fi
fi

if [ "${TEST_SKIP}" == 1 ]; then
    echo "Skipping tests, you are good to push."
    exit 0
else
    echo "Beginning tests"
fi

#Add these removes so the tests don't pass on an old version of the library
rm -f gmml.py _gmml.so
rm -rf ./gmml/lib
if [ -d "./gmml/cmakeBuild" ]; then
    echo "Removing the libgmml.so from our cmakeBuild directory"
    rm ./gmml/cmakeBuild/libgmml.so
fi
echo "Compiling gmml using GEMS ./make.sh, no wrap flag cause it auto wraps"

cd "${GEMS_DIR}" || {
    echo -e "${RED_BOLD}failed...${RESET_STYLE} We could not change directory to the following:\n\t ${GEMS_DIR}"
    echo "Exiting..."
    exit 1
}

./make.sh || {
    echo -e "\n${RED_BOLD}failed... Could not build GEMS and GMML...."
    echo -e "Exiting... ${RESET_STYLE}"
    exit 1
}

cd "${GEMS_DIR}"/gmml || {
    echo -e "${RED_BOLD}failed...${RESET_STYLE} We could not change directory to the following:\n\t ${GEMS_DIR}/gmml"
    echo "Exiting..."
    exit 1
}

echo "Running mandatory tests..."
cd "${GEMS_DIR}"/gmml/tests/ || {
    echo -e "${RED_BOLD}failed...${RESET_STYLE} We could not change directory to the following:\n\t ${GEMS_DIR}/gmml/tests/"
    echo "Exiting..."
    exit 1
}

nice -7 ./compile_run_tests.bash -j "$(nproc --all --ignore=2)"
result=$? # record the exit status from compile_run_tests.bash
cd "${GEMS_DIR}"/gmml || {
    echo -e "${RED_BOLD}failed...${RESET_STYLE} We could not change directory to the following:\n\t ${GEMS_DIR}/gmml"
    echo "Exiting..."
    exit 1
}
if [ "${result}" -eq 0 ]; then
    echo "GMML level tests have passed. Doing gems level tests."
    cd "${GEMS_DIR}"/tests/ || {
        echo -e "${RED_BOLD}failed...${RESET_STYLE} We could not change directory to the following:\n\t ${GEMS_DIR}/tests"
        echo "Exiting..."
        exit 1
    }
    bash run_tests.sh
    gems_tests_result=$? # record the exit status of previous command
    if [ "${gems_tests_result}" -ne 0 ]; then
        echo "GEMS level tests have failed. Make sure you have pulled the latest version and are on the appropriate branch. "
        echo "If you are up-to-date, this failure indicates that you have caused the outputs of ${GEMS_DIR}/tests to change. You can open the ${GEMS_DIR}/tests/run_tests.sh file and run the test line by line to get an output file. Compare it to the saved \"correct\" version in ${GEMS_DIR}/tests/correct_outputs."
        echo "Sometimes the changes you make are fine, and you just need to update what the correct output is by overwriting the old output. Make sure it is ok though, or you will be mur-didely-urdered."
        exit 1
    else
        echo -e "${GREEN_BOLD}All tests have passed. Pushing allowed.${RESET_STYLE}"
        exit 0
    fi
else
    echo -e "${RED_BOLD}
*****************************************************************
The GMML level tests have failed!
Push cancelled.
*****************************************************************
         ${RESET_STYLE}"
    exit 1
fi
