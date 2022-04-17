cmake_minimum_required(VERSION 3.13.4)

#you always have to create a project. This is what your files
#are named, how you create objects, etc.
project(gmml VERSION 1.0.0 DESCRIPTION "GLYCAM Molecular Modeling Library")

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
set(PROJECT_OBJ "${PROJECT_NAME}_obj")

#Now we want to override the typical build type flags. Original ones are as follows:
#1. Release: `-O3 -DNDEBUG`
#2. Debug: `-O0 -g`
#3. RelWithDebInfo: `-O2 -g -DNDEBUG`
#4. MinSizeRel: `-Os -DNDEBUG`
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g")
set(CMAKE_CXX_FLAGS_RELEASE "-O2 -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O0 -DNDEBUG")

#since we want to output .so we need to add lib
#this creates a library with the name lib<project obj name>.so
#well more so a bunch of objects that are made to create a library
#Please note that all *.cc/cpp files must be placed here. Run command
#below in the . directory to make life less painful.
# $find ./src/ -name \*.cc -o -name \*.cpp

#Now if we want to add any fresh files to our system run teh included
#fileListUpdate.sh script
file(READ ${CMAKE_SOURCE_DIR}/cmakeFileLists/cFileList.txt CFILES)
#This changes our delimeter from a newline or spaces to just a semicolon
STRING(REGEX REPLACE ";" "\\\\;" CFILES "${CFILES}")
STRING(REGEX REPLACE "\n" ";" CFILES "${CFILES}")
#add a library that will be an object libary built out of all our cfiles
add_library(${PROJECT_OBJ} OBJECT
	${CFILES}
	)

#Sets whatever will be made from the PROJECT_OBJ to our project version and to
#use the denoted header
set_target_properties(${PROJECT_OBJ} PROPERTIES
	VERSION ${PROJECT_VERSION}
	PUBLIC_HEADER includes/gmml.hpp)

#same as the -I flag for qmake.
file(READ ${CMAKE_SOURCE_DIR}/cmakeFileLists/hDirectoryList.txt HDIRECTORIES)
#this regex is needed because the relative pathing does not work when inputing from the file
#this allows us to mimic relative pathing, also prevents the use of globbing in
#cmake which is a big no-no. We do not glob, it invites anomalous behavior.
STRING(REGEX REPLACE "./" "${CMAKE_SOURCE_DIR}/" HDIRECTORIES "${HDIRECTORIES}")
STRING(REGEX REPLACE ";" "\\\\;" HDIRECTORIES "${HDIRECTORIES}")
STRING(REGEX REPLACE "\n" ";" HDIRECTORIES "${HDIRECTORIES}")
#Now let cmake know that our list is where all the include directories are
#Done this way so we dont have to have individual cmake files cascading through
#Our filesystem, that is annoying and extra work
include_directories(
    "${CMAKE_SOURCE_DIR}"
	"${HDIRECTORIES}"
	)

#Will be used once we get TRUs fixed so we dont run checks on the eigen library
#file(READ ${CMAKE_SOURCE_DIR}/cmakeFileLists/externalHDirectoryList.txt EXDIRECTORIES)
#STRING(REGEX REPLACE "./" "${CMAKE_SOURCE_DIR}/" EXDIRECTORIES "${EXDIRECTORIES}")
#STRING(REGEX REPLACE ";" "\\\\;" EXDIRECTORIES "${EXDIRECTORIES}")
#STRING(REGEX REPLACE "\n" ";" EXDIRECTORIES "${EXDIRECTORIES}")
#include_directories(SYSTEM 
#    "${CMAKE_SOURCE_DIR}"
#    "${EXDIRECTORIES}"
#    )

#for the project object value this passes -fpic which is needed
#to make our files into a shared libary
target_compile_options(${PROJECT_OBJ} PUBLIC "-fPIC" -Wall -W)

#links a library out of our group of objects. Recall, done in this manner
#to easily allow for the preservation of *.o files which is utilized
#by gems during swig wrapping
add_library(${PROJECT_NAME} SHARED $<TARGET_OBJECTS:${PROJECT_OBJ}>)

#Trying to get swig working....

#conditional for user input. To set this variable so we wrap when calling cmake
#we have -DWRAP_GMML(=1? hopefully dont need =1)
if(WRAP_GMML)
	#Get swig package
	find_package(SWIG 4.0 COMPONENTS python REQUIRED)
	#Get python lib package
	find_package(Python3 3.9 COMPONENTS Interpreter Development REQUIRED)

	if (SWIG_FOUND AND Python3_FOUND)
		message(STATUS "Wrapping GMML library in python...")
		#Include ability to actually use our swig
		include(${SWIG_USE_FILE})
		#Include actual dir where python headers reside
		include_directories(${Python3_INCLUDE_DIRS})
		message(STATUS "Swig Version: ${SWIG_VERSION}")
		message(STATUS "Python Version: ${Python3_VERSION}")
		message(STATUS "Python path: ${Python3_INCLUDE_DIRS}\n")

		#have cpp property set on for file gmml.i
		#This property is associated with the gmml.i file itself
		set_property(SOURCE gmml.i PROPERTY CPLUSPLUS ON)
		#add a swig library with name ${PROJECT_NAME}_wrapped, wrapping the code
		#in python, using our gmml.i as our source file for swig
		swig_add_library(${PROJECT_NAME}_wrapped LANGUAGE python SOURCES gmml.i)
		#Now we link the swig library to our gmml libary and our python library
		#so it can actually be useful
		swig_link_libraries(${PROJECT_NAME}_wrapped ${PROJECT_NAME} ${Python3_LIBRARIES})

		#to properly remove swig data
		set_property(TARGET ${PROJECT_NAME}_wrapped PROPERTY 
			ADDITIONAL_MAKE_CLEAN_FILES)

	else()
	#We couldnt find swig or python libs
		if (NOT SWIG_FOUND)
			message(FATAL_ERROR "SWIG WAS NOT FOUND, YOU CANNOT WRAP")
		elseif (NOT Python3_FOUND)
			message(FATAL_ERROR "PYTHON LIBRARY WAS NOT FOUND, YOU CANNOT WRAP")
		else()
			message(FATAL_ERROR "NO IDEA HOW YOU GOT HERE, THIS ISNT GOOD")
		endif()
	endif()
endif()

#NOTE THAT THIS WILL NOT RUN IF WE DO NOT REHIT OUR TARGET, i.e. ALL BUILD IS EXACT SAME
add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
	COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_SOURCE_DIR}/lib/)
#Same issue as above. If we do not populate our lib dir or dont make it that
# means that we did not run any targets for cmake
add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
	COMMAND ${CMAKE_COMMAND} -E copy
	$<TARGET_FILE:${PROJECT_NAME}> ${CMAKE_SOURCE_DIR}/lib/)

#To allow make clean to remove the new lib dir
set_property(DIRECTORY PROPERTY ADDITIONAL_MAKE_CLEAN_FILES
	"../lib")

#TO REMOVE, but kinda good reference
#chuck all .o files to /build dir.
#add_custom_command(TARGET ${PROJECT_NAME}
#  POST_BUILD
#  COMMAND ${CMAKE_COMMAND} -E make_directory ../build/
#  COMMAND ${CMAKE_COMMAND} -E copy_if_different
#  $<TARGET_OBJECTS:${PROJECT_OBJ}> ../build/
#  COMMAND_EXPAND_LISTS)
