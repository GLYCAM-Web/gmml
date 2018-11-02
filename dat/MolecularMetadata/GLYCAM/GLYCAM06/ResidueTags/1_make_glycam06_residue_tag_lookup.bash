#!/bin/bash
##
## File GLYCAM06_make_glycam06_residue_name-2-type_lookup.bash
##
## Begun on 2018-06-18 by BLFoley
##
## Use this script to generate a C++ data representation that can be used to
## map residue types to GLYCAM06 residue names.
##
## It is cast as a bash script so that it will be easy to change it if
## the desired C++ data representation changes

########################################
##
## Declarations / Preliminaries
##
########################################
## Before getting to docs for the main contents, we will make some
## declarations and do some preliminary stuff in a convenient location
##
GMMLOUTPATH='src/MolecularMetadata/GLYCAM'
OUTPATH="../../../../../${GMMLOUTPATH}"
OUTFILE='glycam06residueinfo.cpp' ## name of c++ file to write
if [ -e ${OUTPATH}/${OUTFILE} ] ; then
   echo "

The specified output file:

    ${OUTFILE}

at path:

    ${OUTPATH}

already exists.  Please re/move it and then re-run this script."
fi

########################################
##
## Current representation:  Multimap
##
########################################

########################################
##
##  Types contained herein
##
##      To learn more about the types written by this
##      code, see the file GLYCAM06_res-name-mappings.bash.
##
########################################

##
## Load in the residue tags.
##
. 2_general_residue-tags.bash

##
##  Sanity checking:
##
NumSets=${#TYPES[@]}
if [ "${NumSets}" -ne "${#NAMES[@]}" ] ; then
	echo "

The number of sets of types is:  ${#TYPES[@]}
The number of sets of names is:  ${#NAMES[@]}

These numbers are not the same, but should be.

Exiting.

"
	exit 1
fi

##
## Start by writing header-type info to the hpp file
##
echo "#include \"../../../../includes/MolecularMetadata/GLYCAM/glycam06residueinfo.hpp\"

/** \\file:  ${GMMLOUTPATH}/${OUTFILE}
 * GLYCAM06 metadata for residues
 *
 * This file was generated automatically on:
 *     $(date)
 *
 * by a script named:
 *     ${0}
 *
 * The script was begun on 16 June 2018 by BLFoley and
 * can be found in:
 *     dat/MolecularMetadata/scripts/GLYCAM/GLYCAM06/ResidueTags
 *
 * See that and associated scripts for more information.
 */

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////

using gmml::MolecularMetadata::GLYCAM::Glycam06NamesToTypesLookupContainer;

Glycam06NamesToTypesLookupContainer::Glycam06NamesToTypesLookupContainer()
{
    glycam06NamesToTypesLookupMap_ =
    {" > ${OUTPATH}/${OUTFILE}

##
## Add the data for the map
##
##   Set a nice format:
format="        { %-8s  ,  %-40s } , \n"
NumSets=$((NumSets-1))
i=0
## print out all but the last one
while [ "${i}" -lt "${NumSets}" ] ; do
        if [ "x`printf '%s' "${NAMES[${i}]}" | tr -d "$IFS"`" = x ] ; then
                i=$((i+1))
                continue
        fi
	##
	##   Ensure that the data doesn't contain newlines (for the comments to be right)
	thisNames="${NAMES[${i}]//$'\n'/ }"
	thisTypes="${TYPES[${i}]//$'\n'/ }"
	echo "//    Names:   ${thisNames} " >> ${OUTPATH}/${OUTFILE}
	echo "//    Types:   ${thisTypes} " >> ${OUTPATH}/${OUTFILE}
	for nam in ${NAMES[${i}]} ; do
		for typ in ${TYPES[${i}]} ; do
		printf "${format}" \"${nam}\" \"${typ}\" >> ${OUTPATH}/${OUTFILE}
		done
	done
	i=$((i+1))
done
## Close this data set
echo "
    }; // close Glycam06NamesToTypesLookupMap
}" >> ${OUTPATH}/${OUTFILE}
