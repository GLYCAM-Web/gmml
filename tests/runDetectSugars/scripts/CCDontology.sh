#!/bin/bash

## Make sure the Output Directory is available and if not, make it.
if [ ! -d "${GEMSHOME}/PDB" ]; then
	mkdir "${GEMSHOME}/PDB"
fi

## Set some variables
DEST="${GEMSHOME}/PDB/unprocessed-extract"
RDEST="${GEMSHOME}/PDB/all-reports"
PDBDEST="${GEMSHOME}/PDB/processed"
## Make sure the Output Directories are available and if not, make them.
if [ ! -d "${DEST}" ]; then
	mkdir "${DEST}"
fi
if [ ! -d "${RDEST}" ]; then
	mkdir "${RDEST}"
fi
if [ ! -d "${PDBDEST}" ]; then
	mkdir "${PDBDEST}"
fi

## Create a log file. (Note: I don't think this does anything. It was from a previous script, so I kept it.)
UNPROCESSED="${GEMSHOME}/PDB/all-reports/unprocessed.log"
if [ ! -e "${UNPROCESSED}" ]; then
	touch "${UNPROCESSED}"
fi

## Amino Lib files found in the gmml directory.
AMINO_LIBS="${GEMSHOME}/gmml/dat/CurrentParams/leaprc.ff12SB_2014-04-24/amino12.lib,${GEMSHOME}/gmml/dat/CurrentParams/leaprc.ff12SB_2014-04-24/aminont12.lib,${GEMSHOME}/gmml/dat/CurrentParams/leaprc.ff12SB_2014-04-24/aminoct12.lib,${GEMSHOME}/gmml/dat/CurrentParams/leaprc_GLYCAM_06j-1_2014-03-14/GLYCAM_aminont_06j_12SB.lib,${GEMSHOME}/gmml/dat/CurrentParams/leaprc_GLYCAM_06j-1_2014-03-14/GLYCAM_aminoct_06j_12SB.lib,${GEMSHOME}/gmml/dat/CurrentParams/leaprc_GLYCAM_06j-1_2014-03-14/GLYCAM_amino_06j_12SB.lib"

## Declaring the function to run each PDB.
run_PDB() {

	## Capture the parameter/argument to the function as a variable.
	FILEPATH="$1"
	## Copy the file to the unproccessed-extract directory.
	#cp "${FILEPATH}" "${DEST}"

	## Printing for debugging purposes.
	#echo "FILE: ${FILEPATH}"

	## Set a variable to cut up the FILEPATH.
	##CUT_LENGTH=$(( ${#FILEPATH} - 13 ))

	## Set the ZIPFILE variable and print it for debugging purposes.
	FILENAME="$(basename -- ${FILEPATH})"
	#echo "ZIPFILE: ${ZIPFILE}"

	## Set the FILENAME variable and print it for debugging purposes.
	#FILENAME="$(echo "${ZIPFILE}" | cut -c4-7)"
#if grep -q "$FILENAME" 12-13-20/12-13-20LoopKilled.txt; then
	echo "FILE: ${FILEPATH}"
	#echo "ZIPFILE: ${ZIPFILE}"
	#echo "FILENAME: ${FILENAME}"
	#cp "${FILEPATH}" "${DEST}"

	## Set the output variables.
	#ZIPDEST="$(echo "${DEST}/${ZIPFILE}")"
	#ENTDEST="$(echo "${DEST}/pdb${FILENAME}.ent")"
	#FILENAME="$(echo "${FILENAME}" | tr [:lower:] [:upper:])"
	#FILEDEST="$(echo "${PDBDEST}/${FILENAME}.pdb")"

	## Unzip he gzipped file and move it to it's destination.
	#gunzip "${ZIPDEST}"
	#mv "${ENTDEST}" "${FILEDEST}"
	REPDEST="$(echo "${RDEST}/${FILENAME}.txt")"

	## Print some variables for debugging purposes.
	#echo "FILEDEST: ${FILEDEST}"
	echo "REPDEST: ${REPDEST}"

	#python3 ${GEMSHOME}/testbin/cycle.py -amino_libs ${AMINO_LIBS} -pdb ${FILEDEST} &> ${REPDEST} & pid=$!
        
	
#IF GREP $FILENAME is not in 12-7-20.log then run script	
#-q stops grep from writing output
#	if ! grep -q "$FILENAME" 12-7-20.log; then

		${GEMSHOME}/gmml/tests/detect_sugars ${FILENAME} &> ${REPDEST} 2>&1 & pid=$!
		#&& rm -f ${FILEDEST} >> /dev/null 2&>1 & pid=$!

	#( sleep 60 && kill -HUP $pid ) 2>/dev/null & watcher=$!
	#if wait $pid 2>/dev/null; then
        #   pkill -HUP -P $watcher
        #   wait $watcher
        #else
        #   printf "Failed\n"
        #   echo ${FILEDEST} >> killed.txt
        #fi

	## Some of the PDBs cause our code to get stuck in forever loops.
	## This will kill any straying processes, so the overall purpose of the
	## script can continue. This will log an Error, so we can go back and figure
	## out why.
		time=1
		while kill -0 ${pid} > /dev/null 2>&1; do
			if [ ${time} -gt 2400000 ]; then
				kill -9 ${pid}
				sleep 5
				echo "${FILENAME}" >> killed.txt
 			fi
			((time++))
		done
#fi
} ## End function run_PDB

## MAIN
## Iterate through each argument passed to this script.
for DIRECTORY in "$@"; do

	## Get all the files in the directory and run the function on it.
	for FOLDER in $(ls ${DIRECTORY}); do
		for FILE in $(ls ${DIRECTORY}/${FOLDER}/*.pdb); do
			## Call the function with the PDB file.ent.gz
			run_PDB "${DIRECTORY}/${FOLDER}/${FILE}"
	        done
	done

done

## EXIT_SUCCESS
exit 0
