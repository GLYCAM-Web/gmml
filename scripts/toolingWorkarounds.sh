#!/bin/bash

################################################################
##########               Cool variables             ############
################################################################
#this is to just run all of the needed workarounds automagically, first we are gonna
#need to check if all our needed programs are installed :o), all the clang-tidy etc.
#stuff needs to be built off of clang-14
REQUIRED_PROGRAMS=("clang-14" "clang" "clang-tidy" "clang-format")
##### PRETTY PRINT STUFF #####
#lazy and dont want to have to type all these color variables a bunch
RESET_STYLE='\033[0m'
PASSED_STYLE='\033[0;32m\033[1m'
INFO_STYLE='\033[0;33m\033[1m'
ERROR_STYLE='\033[0;31m\033[1m'


################################################################
##########           Helper function(s/ality)       ############
################################################################
#first go back up to the gmml directory and check where we be
cd ..
if [ "$(basename "${PWD}")" != "gmml" ]; then
    echo -e "${ERROR_STYLE}ERROR SCRIPT RUNNING IN WRONG DIRECTORY"
    echo -e "current directory below:${RESET_STYLE}\n${PWD}"
    exit 1
fi

#restore compile commands to what they used to be
restoreCompileCommands()
{
    echo -e "${INFO_STYLE}###### RESTORING COMPILE_COMMANDS.json TO ORIGINAL STATE ######${RESET_STYLE}"
    if [ ! -f "./compile_commands_BACKUP.json" ]; then
         echo -e "${ERROR_STYLE}ERROR COMPILE_COMMANDS_BACKUP.JSON FILE WAS NOT FOUND${RESET_STYLE}"
         exit 1
    fi
    mv -fv ./compile_commands_BACKUP.json ./cmakeBuild/compile_commands.json || { echo -e "${ERROR_STYLE}ERROR COULDNT RESTORE THE COMPILE_COMMANDS.JSON TO PREVIOUS STATE${RESET_STYLE}"; exit 1; }
    echo -e "${PASSED_STYLE}###### compile_commands.json restored to original state ######${RESET_STYLE}"
    return 0
}

#backup compile commands.json
backupCompileCommands()
{
    ERROR_HAPPENED=0
    if [ ! -d "./cmakeBuild" ]; then
        echo -e "${ERROR_STYLE}ERROR CMAKEBUILD DIR WAS NOT FOUND${RESET_STYLE}"
    fi
    if [ ! -f "./cmakeBuild/compile_commands.json" ]; then
         echo -e "${ERROR_STYLE}ERROR COMPILE_COMMANDS.JSON FILE WAS NOT FOUND${RESET_STYLE}"
    fi
    if [ "${ERROR_HAPPENED}" != 0 ]; then
        echo -e "${ERROR_STYLE}ERROR COULDNT BACK UP NEEDED FILES, EXITING${RESET_STYLE}"
        exit 1
    fi
    echo -e "${INFO_STYLE}###### CREATING BACKUP OF THE COMPILE_COMMANDS.json FILE ######${RESET_STYLE}"
    mv -v ./cmakeBuild/compile_commands.json ./compile_commands_BACKUP.json
    return 0
}

#create a file that contains all headers in our codebase, this is needed for
#some tests due to our include structure being wonky. Keep in mind this does break some tests tho.
#This will 
autoAllHeadersFile()
{
    #switch to src dir so the generated output of the find command be correct
    cd ./src || { echo -e "${ERROR_STYLE}ERROR CANNOT GET INTO SRC DIR${RESET_STYLE}" ; exit 1; }
    #actually get paths of all the header files
    find ../includes -iname "*.hpp" -o -iname "*.h" ! -path "*External_Libraries*" | LC_ALL=C.UTF-8 sort > ./allHeaders.cpp
    #go back to normal gmml dir
    cd ..
    #now include the '#include "' dealio at begining of each line
    sed -i -e 's/^/\#include "/' ./src/allHeaders.cpp
    #now include closing quote, i only do this to hpp endings cause if we index
    #.h files we know something is wrong so i want it to nuke the build if thats
    #the case
    sed -i -e 's/.hpp/\.hpp\"/' ./src/allHeaders.cpp

    #update our cmake file lists and ensure we regenerate our compile_commands.json
    #NOTE: Probably gonna change to like auto_testing or something idk, need to figure out
    #if i want to keep the compile commands + lists with all the test files indexed
    ./make.sh -t "auto"
    trap 'rm -v ./src/allHeaders.cpp && echo -e "${INFO_STYLE}###### Allheaders.cpp workaround file removed.######${INFO_STYLE}"' EXIT
    return 0
}

#beat up the compile_commands.json so it only has the the allHeaders.cpp file
#commands in it, this is for when we are running 
mangleCompileCommands()
{
    echo -e "${INFO_STYLE}###### Making allHeaders.cpp the only thing in compile_commands.json ######${INFO_STYLE}"
    #generate the actual all header file & gen the new compile commands file
    autoAllHeadersFile
    #back up the compile commands
    backupCompileCommands
    
    #first only rip out the json block from the full compile_commands.json and chuck it
    #to where we expect the compile_commands.json to be
    grep '\"file\":.*allHeaders.cpp' ./compile_commands_BACKUP.json -B 3 -A 1 > ./cmakeBuild/compile_commands.json || { echo -e "${ERROR_STYLE}ERROR COULDNT RIP NEEDED LINES OUT OF COMPILE_COMMANDS.JSON${RESET_STYLE}" ; exit 1; }
    #add the opening bracket and newline at the beginning of the file
    sed -i '1s/^/[\n/' ./cmakeBuild/compile_commands.json
    #add closing bracket
    sed -i '$a]' ./cmakeBuild/compile_commands.json
    #remove the comma at the end of the bracket from the lines we ripped out
    sed -i "s/},/}/" ./cmakeBuild/compile_commands.json
    
    echo -e "${PASSED_STYLE}###### allHeaders.cpp now the only thing in compile_commands.json ######${INFO_STYLE}"
    return 0
}

#repair all header guards and take care of the needed workarounds. This will
#check for and remove duplicate includes, sort said includes, repair the header guards
#themselves, add in comments to the header guard macros
#NOTE: NEED A BETTER NAME!!!
repairHeaders()
{
    #go ahead and get whats needed to run on all our header files
    mangleCompileCommands
    
    #get all OUR header files and when we find em we add a temp comment ad the end of
    #any endif macro directive. This is a little sloppy, and does cause a couple extra
    #comments to be made in places they shouldnt but we will leave that be for now.....
    find ./includes -iname "*.hpp" -print0 | xargs -0 sed -i "s/\#endif$/\#endif \/\/ TEMP_COMMENT/g"
    #First we have to make sure the header guards are actually proper
    run-clang-tidy -checks='-*, llvm-header-guard' -p ./cmakeBuild/ -header-filter=.hpp -fix || \
        { echo -e "${ERROR_STYLE}ERROR COULDNT APPLY REPAIR HEADERS CHANGES${RESET_STYLE}" ; exit 1; }
    #Now, since our directory structure is not something that LLVM auto stuff
    #expects we gotta do some gross stuff
    BASE_PATH_AS_GUARD_STYLE="$(cd .. && echo "${PWD}" && cd ./gmml/)" || { echo -e "${ERROR_STYLE}ERROR CANNOT GET BASE LINE PATH${RESET_STYLE}" ; exit 1; }
    #Make all capital letters
    BASE_PATH_AS_GUARD_STYLE="${BASE_PATH_AS_GUARD_STYLE^^}"
    #now replace all forward slashs with underscores, using # as the operator
    #seperator in sed so its easier to look at
    BASE_PATH_AS_GUARD_STYLE=$(echo "${BASE_PATH_AS_GUARD_STYLE}" | sed -e "s#\/#\_#g" | sed -e "1s#\_##")
    
    #this whole block of find then sedding basically just removes the extra long
    #total system path from the pattern. Has to be done cause our code structure
    #is not what the auto tooling expects
    find . -type f -iname "*.hpp" -print0 | xargs -0 sed -i "s/\#ifndef.*${BASE_PATH_AS_GUARD_STYLE}_/\#ifndef /"
    find . -type f -iname "*.hpp" -print0 | xargs -0 sed -i "s/\#define.*${BASE_PATH_AS_GUARD_STYLE}_/\#define /"
    find . -type f -iname "*.hpp" -print0 | xargs -0 sed -i "s/\#endif.*\/\/.*${BASE_PATH_AS_GUARD_STYLE}_/\#endif \/\/ /"
    
    #once our header guards are good we then ensure all duplicate includes are removed
    #it is only the llvm dupe remover and it is pretty dumb, eventuall we wanna
    #use IWYU
     run-clang-tidy -checks='-*, readability-duplicate-include' -p ./cmakeBuild/ -header-filter=.hpp -fix || \
        { echo -e "${ERROR_STYLE}ERROR COULDNT APPLY REPAIR HEADERS CHANGES${RESET_STYLE}" ; exit 1; }
    
    #Now we sort the order of includes just so it looks nicer and more consistent
     run-clang-tidy -checks='-*, llvm-include-order, readability-braces-around-statements, llvm-namespace-comment' -p ./cmakeBuild/ -header-filter=.hpp -fix || \
        { echo -e "${ERROR_STYLE}ERROR COULDNT APPLY REPAIR HEADERS CHANGES${RESET_STYLE}" ; exit 1; }
    
    
    #Now remove all the temp comments left over lol 
    find ./includes -iname "*.hpp" -print0 | xargs -0 sed -i "s/\#endif \/\/ TEMP_COMMENT$/\#endif/g"
    
    restoreCompileCommands
    return 0
}

#now this does some gentle fixes for our actual source files
repairSources()
{
    #first make sure that our compile_commands.json etc. is nice and goode
    ./make.sh -t "auto"
    
    run-clang-tidy -checks='-*, readability-duplicate-include' -p ./cmakeBuild/ -fix || \
        { echo -e "${ERROR_STYLE}ERROR COULDNT APPLY REPAIR SOURCES CHANGES${RESET_STYLE}" ; exit 1; }
    
    #Now we sort the order of includes just so it looks nicer and more consistent
     run-clang-tidy -checks='-*, llvm-include-order, readability-braces-around-statements, llvm-namespace-comment' -p ./cmakeBuild/ -fix || \
        { echo -e "${ERROR_STYLE}ERROR COULDNT APPLY REPAIR SOURCES CHANGES${RESET_STYLE}" ; exit 1; }
    
    run-clang-tidy -p ./cmakeBuild/ -fix
}


#run the actual clang formatter on all code, since its a simple command its 
#basically a wrapper. Does keep code more understandable tho
formatAsOne()
{
    #we need to take advantage of word splitting so the clang prog is given a list
    #of files instead of one long file name
    # shellcheck disable=SC2046
    clang-format -i $(find ./src ./includes ./tests -type f -iname "*.cc" -o -iname "*.cpp" -o -iname "*.hpp")
}

#loop through all reqs and check if we got em, if not we exit
MISSING_PROGS=0
for CURR_PROG_CHECK in "${REQUIRED_PROGRAMS[@]}"; do
    command -v "${CURR_PROG_CHECK}" | grep -o "${CURR_PROG_CHECK}" > /dev/null || { echo -e "${ERROR_STYLE}ERROR MISSING PROGRAM:${CURR_PROG_CHECK}${RESET_STYLE}" ; MISSING_PROGS=1; }
done
if [ "${MISSING_PROGS}" -gt 0 ]; then
    echo -e "\n${ERROR_STYLE}ERROR MISSING REQUIRED PROGRAMS, GO INSTALL THE"
    echo -e "PROGRAMS LISTED ABOVE IN ORDER TO RUN THIS CODE."
    exit 1
fi

################################################################
##########               Print help message         ############
################################################################

printHelp()
{
echo -e "
===== UPDATE THIS LOLE =====
$0 is used to allow us to run multiple of our
test scripts at once. This is not to be confused with compiling
the test code, we will be running multiple of our scripts at once.
There exists some GNU utils that could be used to accomplish this
but we want to keep these bash scripts as bare bones as possible.

We find tests to run by globbing \"*.test.*.sh\", to disable a
test just change the file name of test you would like to disable
so it does not match that glob.
*************************************************************
Options are as follows:
\t-j <NUM_JOBS>\t\tRun <NUM_JOBS> scripts at once
\t-h \t\t\tPrint this msg
*************************************************************
Exiting"
exit 1
}

################################################################
#########               COMMAND LINE INPUTS            #########
################################################################
#get opts allows for users to pass flags to the script, i.e. calling
#./compile_run_tests.bash -j 11 will allow us to run the script with 11 jobs.
#The colon after j (in the while decleration) means that said flag expects an "input"
#The lack of a colon after h means that said flag does not accept any "input"
while getopts "r:h" option
do
    case "${option}" in
        #NOTE: GONNA NEED TO CHANGE THIS PATTERN, SINCE SO FRAGILE AS IS
        #CANT JUST HAVE TESTS RUN OUT OF A SPECIFIC ORDER. ENFORCE STEPS
        #WE TAKE!!!
        #to RUN a specific tool blurb
        r)
            rIn="${OPTARG}"
            rIn="${rIn,,}"
            case "${rIn}" in
                headerstime)
                    repairHeaders
                    ;;
                sourcestime)
                    repairSources
                    ;;
                formattime)
                    formatAsOne
                    ;;
                *)
                    printHelp
                    ;;
            esac
            ;;
        h)
            printHelp
            ;;
        #we need a wildcard case so that our script will print out the help msg
        #if a user passes an unexpected flag
        *)
            printHelp
            ;;
    esac
done


echo "Cool script started"



echo "cool script ended"
exit 0
