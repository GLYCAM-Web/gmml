#!/bin/bash

################################################################
##########               Cool variables             ############
################################################################
GMML_TEST_JOBS=1
#Get base list of all files we want so we dont need to deal with
#figuring out file list more than once, we want to take advantage of output splitting 
#here so we can hit the list ez pz
# shellcheck disable=SC2207
GMML_TEST_FILE_LIST=($(ls ./*.test.*.sh))
#This is mostly to make sure that we actually keep track of our failed tests. 
GMML_FAILED_TESTS=0
#Just to keep track for our cool output at the end.
GMML_PASSED_TESTS=0
#amount of tests we are actually running, thus amount of tests we need to pass
#GMML_REQUIRED_PASSING=$(find . -type f -name "*.test.*.sh" | wc -l)

##### PRETTY PRINT STUFF #####
#lazy and dont want to have to type all these color variables a bunch
RED='\033[0;31m' 
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
BOLD_STYLE='\033[1m'
#Done by combo of these 2: NORMAL_STYLE='\033[0m' & NO_COLOR='\033[0m'
RESET_STYLE='\033[0m\033[0m'

################################################################
##########               Print help message         ############
################################################################

printHelp()
{
        echo ""
        echo "==== GMML TEST RUNNING SCRIPT ===="
        echo "$0 is used to allow us to run multiple of our test scripts at"
        echo "once. This is not to be confused with compiling the test code,"
        echo "we will be running multiple of our scripts at once. There are"
        echo "some light GNU utils that could be used to make this cleaner but"
        echo "we want to keep these bash scripts as bare bones as possible."
        echo ""
        echo "We find tests to run by globbing \"*.test.*.sh\", to disable a"
        echo "test just change the file name of test you would like to disable"
        echo "so it does not match that glob."
        echo "*************************************************************"
        printf "Options are as follows:\n"
        printf "\t-j <NUM_JOBS>\t\tRun <NUM_JOBS> scripts at once\n"
        printf "\t-h \t\t\tPrint this msg\n"
        echo "*************************************************************"
        echo ""
        echo "Exiting"
        exit 1
}

################################################################
#########               COMMAND LINE INPUTS            #########
################################################################

while getopts "j:h" option
do
    case "${option}" in
        j)
            jIn=${OPTARG}
            if [[ ${jIn} =~ ^[0-9]+$ ]]; then
                GMML_TEST_JOBS=${jIn}
            else
                printHelp
            fi
            ;;
        h)
            printHelp
            ;;
        *)
            printHelp
            ;;
    esac
done

echo ""
echo -e "${YELLOW}${BOLD_STYLE}#### Beginning GMML tests ####${RESET_STYLE}"
echo "Number of tests found: ${#GMML_TEST_FILE_LIST[@]}"
echo "Begining testing using ${GMML_TEST_JOBS} jobs."
echo ""

mkdir -v ./tempTestOutputs
echo ""
#always delete this dir before we exit for any reason, get rid of em cause no useful info in em
trap "rm -rf ./tempTestOutputs" EXIT

#Gonna be used in the main loop, here so its close to first usage.
JOB_PIDS=()
JOB_OUTPUT_FILES=()

cleaningUpJobs()
{
    #first check if we be in sync, if not, we abort and freak out
    if  [ ${#JOB_PIDS[@]} != ${#JOB_OUTPUT_FILES[@]} ] ; then
        echo -e "${RED}${BOLD_STYLE}!!!! WARNING PID LIST AND FILE LIST OUT OF SYNC ABORTING !!!!${RESET_STYLE}"
        #nuke all subshells
        kill 0
        exit 1
    fi
    #now loop through our array, check the PIDs in our job array, if job
    #is done we nuke the data in the array, note we start at the end of the
    #array cause 0 is our base line index, so mutating array size wont cause
    #a bork if we do it this way (hopefully)
    for (( CURR_INDEX=${#JOB_PIDS[@]}-1; CURR_INDEX>=0; CURR_INDEX-- )) ; do
        #if our process is NOT running then we go ahead and get rid
        #of both the PID in the array and the file being held onto in
        #our output file array
        if ! ps -p "${JOB_PIDS[${CURR_INDEX}]}" > /dev/null ; then
            #now we need to extract the exit code from the file itself, this will 
            #be the last line of the file. I store as a variable cause we know we
            #gonna have to run the tail command at least once, and possibly more,
            #so just take the hit once and we are good. Gotta use 2 cause include newline
            DUMB_EXIT_CODE=$(tail -c -2 "${JOB_OUTPUT_FILES[${CURR_INDEX}]}")            
            if [[ ${DUMB_EXIT_CODE} == 0 ]]; then
                ((GMML_PASSED_TESTS++))
            elif [[ ${DUMB_EXIT_CODE} == 1 ]]; then
                ((GMML_FAILED_TESTS++))
            else
                echo -e "${RED}${BOLD_STYLE}!!!! WARNING: EXIT CODE LINE INCORRECT, EXITING WHOLE SCRIPT !!!!${RESET_STYLE}"
                #nuke all subshells
                kill 0
                exit 1
            fi
            echo ""
            #color and output the jobs output
            GREP_COLOR="1;31" grep --color=always '.*FAILED.*\|$' "${JOB_OUTPUT_FILES[${CURR_INDEX}]}" | GREP_COLOR="1;32" grep --color '.*passed.*\|$'
            
            #Now remove the PID from our "scheduler" array, ngl more of a tracker
            unset JOB_PIDS["${CURR_INDEX}"]
            #Now remove the job output from list
            unset JOB_OUTPUT_FILES["${CURR_INDEX}"]
            
            #Now we have to clean up our arrays, this kinda sucks and am unsure better way to fix
            JOB_PIDS=("${JOB_PIDS[@]}")
            JOB_OUTPUT_FILES=("${JOB_OUTPUT_FILES[@]}")
        fi
    done
}
                   

for CURRENT_TEST in "${GMML_TEST_FILE_LIST[@]}" ; do 
    #When we have all our wanted jobs running, we wait for 1 or more to exit then clean up our finished jobs
    if [ ${#JOB_PIDS[@]} == "${GMML_TEST_JOBS}" ] && [ ${#JOB_OUTPUT_FILES[@]} == "${GMML_TEST_JOBS}" ] ; then
        #This waits for 1, or more, job to complete
        #Once a job completes, we then need to go ahead and clean up our dirty job scheduler/tracker
        wait -n
        cleaningUpJobs
        echo ""
    fi
    #let user know whats going on
    echo -e "${YELLOW}${BOLD_STYLE}Beginning test:${RESET_STYLE} ${CURRENT_TEST}"
    
    #Wanted to do below but couldnt figure out force check that the only .sh
    #we would be replacing would be at the very end like I can with sed....
    #CURRENT_TEST_OUTPUT_FILENAME="${CURRENT_TEST//.sh/.output.txt}"
    
    #do some sed stuff to make our filename, this is a POSIX compliant method
    CURRENT_TEST_OUTPUT_FILENAME="./tempTestOutputs/"$(echo "${CURRENT_TEST}" | sed -e "s/\.sh$/.output.txt/")
  
    #actually run the script we currently want to run
    source "${CURRENT_TEST}" &> "${CURRENT_TEST_OUTPUT_FILENAME}" &
    
    #NOTE: Both below MUST ALWAYS BE IN SYNC, lazy way of doing some
    #instability to try and catch bad behavior..... bad idea? who knows!
    #Add the most recent job PID to our PID array
    JOB_PIDS+=($!)
    #Add corresponding output file to our output file array
    JOB_OUTPUT_FILES+=("${CURRENT_TEST_OUTPUT_FILENAME}")
done
#wait for all remaining files in pid to complete then clean....
#this has to be done because the wait -n just waits for one or more jobs
#to complete, at the end of the for loop there is a possibility of us having
#other jobs in the process of completing so we just wait for all the remaining
#jobs then we run the clean up then we done.
wait
cleaningUpJobs

case "${GMML_PASSED_TESTS}" in
    ${#GMML_TEST_FILE_LIST[@]})
        RESULT_COLOR=${GREEN}${BOLD_STYLE}
        ;;
    *)
        RESULT_COLOR=${RED}${BOLD_STYLE}
        ;;
esac

echo -e """
${RESULT_COLOR}
#### GMML TESTS COMPLETED ####
Required tests: ${#GMML_TEST_FILE_LIST[@]}
Passed tests:   ${GMML_PASSED_TESTS}
Failed tests:   ${GMML_FAILED_TESTS}${RESET_STYLE}"""

#Ps paranoid programming
if [ "${GMML_PASSED_TESTS}" == "${#GMML_TEST_FILE_LIST[@]}" ] ; then 
    echo -e "${GREEN}${BOLD_STYLE}GMML TESTS PASSED${RESET_STYLE}\n"
    #delete output files here
    exit 0
elif [[ $(( GMML_PASSED_TESTS + GMML_FAILED_TESTS )) != "${#GMML_TEST_FILE_LIST[@]}" ]] ; then
    echo -e "${RED}${BOLD_STYLE}!!! ERROR WE DIDNT GET EXIT CODES FROM ALL NEEDED SCRIPTS !!!${RESET_STYLE}\n"
    exit 1
elif [ "${GMML_FAILED_TESTS}" != 0 ] ; then
    echo -e "${RED}${BOLD_STYLE}!!! ${GMML_FAILED_TESTS} GMML TESTS FAILED !!!${RESET_STYLE}\n"
    exit 1
else
    echo -e "${RED}${BOLD_STYLE}!!! SOMETHING BORKED OH NO NOT GOOD !!!${RESET_STYLE}\n"
    exit 1
fi
