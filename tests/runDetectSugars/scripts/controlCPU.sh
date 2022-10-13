#!/bin/bash

## This script is used to run all of the PDBs in chunks of the value of the CPUS variable below.
## For example: if CPUS is 16, it will divide all the directories up into 16 lists of directories
##  to send the the ontology.sh script to iterate through the directories to run one at a time per
##  ontology.sh call. Which will be the size of the CPUS variable below.

USAGE="""
Usage: controlCPU.sh DIRECTORY [CPUS]
DIRECTORY is the path to a directory with one or more subdirectory of PDBs
Example: /programs/repos/PDB/all contains many subdirectories of PDBs in the form of pdb1atn.ent.gz
"""

## Save first argument to SRC
SRC="$(readlink -f $1)"

## Check if SRC is a directory.
if [ ! -d "${SRC}" ]; then
	echo -e "Error: '${SRC}': Not an existing directory:\n${USAGE}" >&2
	exit 1
fi

## Save second argument to CPUS
CPUS="$2"

## Check if CPUS is a valid unsigned integer. If not, set to default value of 1.
if [ "${CPUS}" == "" ] || ! [[ ${CPUS} =~ ^[0-9]+$ ]] ; then
	# 1 is the default because my script to run the archive uses 20 CPU cores to build by distance
	CPUS=1
fi

if [ ! -e $(pwd)/ontology.sh ]; then
	echo "It appears you don't have the partner script, ontology.sh. It is needed to run the rest."
	exit 1
fi
## Initialize an array to store all the subdirectories into.
DIRECTORIES_ARRAY=()

count=0
for DIRECTORY in $(ls ${SRC}); do
	DIRECTORIES_ARRAY[${count}]="${SRC}/${DIRECTORY}"
	# echo "${SRC}/${DIRECTORY}"
	((count++))
done

total_directories=${count}
processes=$((${total_directories}/${CPUS}))
((processes++))

for ((i=0; i<${#DIRECTORIES_ARRAY[@]}; )); do
	for ((j=0; j<${processes} && i<${#DIRECTORIES_ARRAY[@]}; i++, j++)); do
		DIRECTORIES+="${DIRECTORIES_ARRAY[${i}]} "
	done
	## Pass DIRECTORIES to sub-script. DIRECTORIES should be a string of ${processes} directories to then run 1 PDB at a time.
	$(pwd)/ontology.sh ${DIRECTORIES} &
	#echo ${DIRECTORIES}
	unset DIRECTORIES
done
