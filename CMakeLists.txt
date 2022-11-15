#There is a little cruft in here, but it is things that were done in
#order to make the previous build logic work. The cruft is depricated
#but all the wonky workarounds will hopefully be somewhat useful in the
#future.
cmake_minimum_required(VERSION 3.13.4)

#you always have to create a project. This is what your files
#are named, how you create objects, etc.
project(gmml VERSION 1.0.0 DESCRIPTION "GLYCAM Molecular Modeling Library")

################## BEGIN GETTING TOOLING TOGETHER ###################

#Exit if user doesnt have a new enough compiler to use c++17
#Need to check its at least version 7.0 cause that is supposed to
#be when c++17 was full supported by g++
if (NOT CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 7.0)
    message(WARNING "Compiler used: ${CMAKE_CXX_COMPILER_ID}")
    message(WARNING "Compiler version used: ${CMAKE_CXX_COMPILER_VERSION}")
    message(WARNING "Error, your compiler version is a high enough version")
    #fatal error msg causes it to stop and exit
    message(FATAL_ERROR "Check gcc version is at least 7.0")
endif()

#Basically we always hit here UNLESS we pass the DEBUG_NO_WRAP=True. This allows us to produce
#a makefile and other assoicaited files that have zero reference to swig, this helps for
#debugging purposes
if (NOT  ${DEBUG_NO_WRAP})
    #Get swig package
    find_package(SWIG 4.0 COMPONENTS python REQUIRED)
    #Get python lib package
    find_package(Python3 3.9 COMPONENTS Interpreter Development EXACT REQUIRED)
endif()

#Show user what compiler they using
message("################## TOOLS CMAKE IS USING ###################")
message(STATUS "CMake Version:\t${CMAKE_VERSION}")
message("")
#https://cmake.org/cmake/help/latest/variable/CMAKE_LANG_COMPILER_ID.html
message(STATUS "C++ Compiler Used:\t${CMAKE_CXX_COMPILER_ID}")
message(STATUS "C++ Compiler Ver:\t${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "C++ Compiler Path:\t${CMAKE_CXX_COMPILER}")
message("")
message(STATUS "C Compiler Used:\t${CMAKE_C_COMPILER_ID}")
message(STATUS "C Compiler Version:\t${CMAKE_C_COMPILER_VERSION}")
message(STATUS "C Compiler Path:\t${CMAKE_C_COMPILER}")
#Since we dont check for these packages if are are doing the total
#wrap skip thing, we also want to skip over places that try to reference
#packages we are not sure we have
if (NOT ${DEBUG_NO_WRAP})
    message("")
    message(STATUS "Swig Executable:\t${SWIG_EXECUTABLE}")
    message(STATUS "Swig Version:\t${SWIG_VERSION}")
    message(STATUS "Swig Library Dir:\t${SWIG_DIR}")
    message("")
    message(STATUS "Python Interpreter:\t${Python3_EXECUTABLE}")
    message(STATUS "Python Version:\t${Python3_VERSION}")
    message(STATUS "Python Include Dir:\t${Python3_INCLUDE_DIRS}")
    message(STATUS "Python Library Dir:\t${Python3_LIBRARY_DIRS}")
    #actually let us use swig tool
    include(${SWIG_USE_FILE})
endif()
message("###########################################################\n")
################## END GETTING TOOLING TOGETHER ###################


################## BEGIN SETTING BASE COMPILE FLAGS ###################
#export compile commands for sourcetrail, etc.
#it is just a json of how commands interact, most tools use this
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
#the find package that is commented out a bit below just checks your
#system for the pthreads package then chucks it in. It kinda works
#but setting the flag is easier, it just passes the -pthred flag
#to gcc, nice and easy. You can add more custom flags here
set(CMAKE_CXX_FLAGS -pthread)

#in order to have 2 seperate values to keep track of.
#done in order to move our *.o files to a directory that will keep
#the swig wrapping in gems happy
#set(PROJECT_OBJ "${PROJECT_NAME}_obj")

#Now we want to override the typical build type flags. Original ones are as follows:
#1. Release: `-O3 -DNDEBUG`
#2. Debug: `-O0 -g`
#3. RelWithDebInfo: `-O2 -g -DNDEBUG`
#4. MinSizeRel: `-Os -DNDEBUG`
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g")
set(CMAKE_CXX_FLAGS_RELEASE "-O2 -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O0 -DNDEBUG")
################## END SETTING BASE COMPILE FLAGS ###################

#since we want to output .so we need to add lib
#this creates a library with the name lib<project obj name>.so
#well more so a bunch of objects that are made to create a library
#Please note that all *.cc/cpp files must be placed here. Run command
#below in the . directory to make life less painful.
# $find ./src/ -name \*.cc -o -name \*.cpp

################## BEGIN READING THE FILE LISTS ###################
#Now if we want to add any fresh files to our system run teh included
#fileListUpdate.sh script
file(READ ${CMAKE_SOURCE_DIR}/cmakeFileLists/cFileList.txt CFILES)
#This changes our delimeter from a newline or spaces to just a semicolon
STRING(REGEX REPLACE "\\./" "${CMAKE_SOURCE_DIR}/" CFILES ${CFILES})
STRING(REGEX REPLACE ";" "\\\\;" CFILES "${CFILES}")
STRING(REGEX REPLACE "\n" ";" CFILES "${CFILES}")
STRING(REGEX REPLACE "[\n;]$" "" CFILES "${CFILES}")

#same as the -I flag for gcc.
file(READ ${CMAKE_SOURCE_DIR}/cmakeFileLists/hDirectoryList.txt HDIRECTORIES)
#this regex is needed because the relative pathing does not work when inputing from the file
#this allows us to mimic relative pathing, also prevents the use of globbing in
#cmake which is a big no-no. We do not glob, it invites anomalous behavior.
STRING(REGEX REPLACE "\\./" "${CMAKE_SOURCE_DIR}/" HDIRECTORIES "${HDIRECTORIES}")
STRING(REGEX REPLACE ";" "\\\\;" HDIRECTORIES "${HDIRECTORIES}")
STRING(REGEX REPLACE "\n" ";" HDIRECTORIES "${HDIRECTORIES}")
STRING(REGEX REPLACE "[\n;]$" "" HDIRECTORIES "${HDIRECTORIES}")

#same as the -isystem flag for gcc. This is for all our external header file dirs
#important cause tooling etc. 
file(READ ${CMAKE_SOURCE_DIR}/cmakeFileLists/externalHDirectoryList.txt EXHDIRECTORIES)
STRING(REGEX REPLACE "\\./" "${CMAKE_SOURCE_DIR}/" EXHDIRECTORIES "${EXHDIRECTORIES}")
STRING(REGEX REPLACE ";" "\\\\;" EXHDIRECTORIES "${EXHDIRECTORIES}")
STRING(REGEX REPLACE "\n" ";" EXHDIRECTORIES "${EXHDIRECTORIES}")
STRING(REGEX REPLACE "[\n;]$" "" EXHDIRECTORIES "${EXHDIRECTORIES}")
################## END READING THE FILE LISTS ###################


################## BEGIN SETTING UP GMML TARGET ###################
#add a library that will be an object libary built out of all our cfiles, this will
#be our baseline normal gmml that isnt wrapped
add_library(${PROJECT_NAME} SHARED
	${CFILES})

#Add the headers to our specific gmml target
target_include_directories(${PROJECT_NAME} PUBLIC
	"${CMAKE_SOURCE_DIR}" "${HDIRECTORIES}")
#Declaring the system will cause the compile commands to have the 
#prefix -isystem which should cause the autotooling to ignore these files
target_include_directories(${PROJECT_NAME} SYSTEM PUBLIC
    "${EXHDIRECTORIES}"
    )
	
#Specifies the header file we are gonna make into a shared object file,
#since our gmml.hpp is a god file we just include that
set_target_properties(${PROJECT_NAME} PROPERTIES
	PUBLIC_HEADER "includes/gmml.hpp")
    
#for the project object value this passes -fpic which is needed
#to make our files into a shared libary
target_compile_options(${PROJECT_NAME} PUBLIC "-fPIC" -Wall -W)
################## END SETTING UP GMML TARGET ###################

#we only skip here if we are trying to gen a makefile that has 0 wrap
#acapabilities. This is needed to help linting, etc. run because we dont generate
#a compile commands that includes the wrapping commands. Also if we dont even
#have the swig stuff included in our cmake list file none of these commands
#will work.
if (NOT ${DEBUG_NO_WRAP})
    ################## BEGIN SETTING UP GMML_WRAPPED TARGET ###################
    #have cpp property set on for file gmml.i
    #This property is associated with the gmml.i file itself
    set_property(SOURCE gmml.i PROPERTY CPLUSPLUS ON)
    #add a swig library with name ${PROJECT_NAME}_wrapped, wrapping the code
    #in python, using our gmml.i as our source file for swig
    swig_add_library(${PROJECT_NAME}_wrapped LANGUAGE python SOURCES gmml.i)
    #makes sure the actual module produced by swig for wrapping is good, cmake 3.14
    #forward renames the module to whatever the target name is i.e. gmml_wrapped
    set_target_properties(${PROJECT_NAME}_wrapped PROPERTIES OUTPUT_NAME gmml)
    #add the python headers so we can actually use them and wrap gmml in python
    target_include_directories(${PROJECT_NAME}_wrapped PUBLIC ${Python3_INCLUDE_DIRS})

    #Now we link the swig library to our gmml libary and our python library
    #so it can actually be useful
    swig_link_libraries(${PROJECT_NAME}_wrapped ${PROJECT_NAME} ${Python3_LIBRARIES})
    ################## END SETTING UP GMML_WRAPPED TARGET ###################
endif()

#links a library out of our group of objects. Recall, done in this manner
#to easily allow for the preservation of *.o files which is utilized
#by gems during swig wrapping
#add_library(${PROJECT_NAME} SHARED $<TARGET_OBJECTS:${PROJECT_OBJ}>)

################## BEING COPYING BUILT FILES TO NEEDED DIRS ###################
#NOTE THAT THIS WILL NOT RUN IF WE DO NOT REHIT OUR TARGET, i.e. ALL BUILD IS EXACT SAME
add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
	COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_SOURCE_DIR}/lib/)
#Same issue as above. If we do not populate our lib dir or dont make it that
# means that we did not run any targets for cmake
add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
	COMMAND ${CMAKE_COMMAND} -E copy
	$<TARGET_FILE:${PROJECT_NAME}> ${CMAKE_SOURCE_DIR}/lib/)
################## END COPYING BUILT FILES TO NEEDED DIRS ###################


################## BEGIN ALL BUILT FILES INCLUDED IN `make clean` ###################
#to properly remove swig data. As with other swig stuff we skip this when we 
#dont wanna wrap.
if (NOT ${DEBUG_NO_WRAP})
    set_property(TARGET ${PROJECT_NAME}_wrapped PROPERTY 
        ADDITIONAL_MAKE_CLEAN_FILES)
endif()

#To allow make clean to remove the new lib dir
set_property(DIRECTORY PROPERTY ADDITIONAL_MAKE_CLEAN_FILES
	"../lib")
################## END ALL BUILT FILES INCLUDED IN `make clean` ###################

#TO REMOVE, but kinda good reference
#chuck all .o files to /build dir.
#add_custom_command(TARGET ${PROJECT_NAME}
#  POST_BUILD
#  COMMAND ${CMAKE_COMMAND} -E make_directory ../build/
#  COMMAND ${CMAKE_COMMAND} -E copy_if_different
#  $<TARGET_OBJECTS:${PROJECT_OBJ}> ../build/
#  COMMAND_EXPAND_LISTS)
