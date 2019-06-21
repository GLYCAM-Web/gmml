#!/bin/bash
##
########################################
## File GLYCAM06_res-name-mappings.bash
##
## Begun on 2018-06-19 by BLFoley
## 
## To be used with 
##
##    GLYCAM06_make_glycam06_residue_name-2-type_lookup.bash
##
## This file contains the actual mappings to keep the other code
## easier to read.
##
## Want to DRY this code out?  Oh, thank you!!!!  Do it!!!!
##
########################################

## 
## This function generates these numbered arrays:
##
declare -a TYPES
declare -a NAMES
##
## The indexing in these arrays starts with zero:
##
i=0
##

########################################
##
##  About the Types 
##
########################################
##
##  Lists of types used:
##
##         Before each section, the most common types used in 
##         that section are listed.
##
##         Notes:  
##
##         *  These lists will certainly change for future needs.  
##            Consider any current lists to be approximations only.
##         *  The lists might not be complete in any sense of the
##            term.  That is.
##            -- The data used to make the code might contains 
##               elements not on the lists.
##            -- The data needed for full representation miht
##               not be in the list or in the code.
##
##  Regarding the data representation in the C++ code:
##
##         *  Currently, there is no hierarchy of types.  That is, 
##            'C-linking' and 'amino-acid' are just part of the list.
##
##         *  Groupings below are for the convenience of humans to read.
##
##  Types are not exclusive:
##
##         The presence of any one type does not keep another 
##         type from also being present.  
##
##         For example:
##         ZOLT and ZOLS are both aglycon and amino-acid.
##
##  List format:
##
##  Type
##              Sub-types
##
###################################

###################################
##  Amino Acids
###################################
##
##  amino-acid
##              C-linking
##              N-linking
##              O-linking
##              mid-chain
##              C-termial
##              N-terminal
##              zwitterionic
##              C-terminal-cap
##              N-terminal-cap
##
###################################
##
## HYP and variants
TYPES[${i}]=" amino-acid modified mid-chain "
NAMES[${i}]=" HYP "
i=$((i+1))
TYPES[${i}]=" amino-acid modified N-terminal "
NAMES[${i}]=" CHYP "
i=$((i+1))
TYPES[${i}]=" amino-acid modified C-terminal "
NAMES[${i}]=" NHYP "
##
## Zwitterionic aglycons
i=$((i+1))
TYPES[${i}]=" amino-acid O-linking zwitterionic aglycon"
NAMES[${i}]=" ZOLS  ZOLT "
##
## N-linking to peptides/proteins
i=$((i+1))
TYPES[${i}]=" amino-acid N-linking mid-chain "
NAMES[${i}]=" NLN   "
i=$((i+1))
TYPES[${i}]=" amino-acid N-linking C-terminal "
NAMES[${i}]=" CNLN  "
i=$((i+1))
TYPES[${i}]=" amino-acid N-linking N-terminal "
NAMES[${i}]=" NNLN  "
##
## O-linking to peptides/proteins
i=$((i+1))
TYPES[${i}]=" amino-acid O-linking mid-chain "
NAMES[${i}]=" OLP  OLS  OLT  OLY   "
i=$((i+1))
TYPES[${i}]=" amino-acid O-linking C-terminal "
NAMES[${i}]=" COLP COLS COLT COLY "
i=$((i+1))
TYPES[${i}]=" amino-acid O-linking N-terminal "
NAMES[${i}]=" NOLP  NOLS  NOLT  NOLY "

###################################
##  Aglycons
###################################
##
##  aglycon
##
##  monosaccharides that can be aglycons are added elsewhere in
##      these scripts
##
###################################

i=$((i+1))
TYPES[${i}]=" aglycon "
NAMES[${i}]=" OME ROH TBT "

###################################
##  Counter Ions
###################################
##
##  counter-ion
##              formal-charge=fc
##                     fc is the integer representing the formal charge
##
###################################

i=$((i+1))
TYPES[${i}]=" counter-ion formal-charge=+2  "
NAMES[${i}]=" CA2 "

###################################
##  Derivatives
###################################
##
##  derivative
##              formal-charge=fc
##                     fc is the integer representing the formal charge
##              head-hybridization=sp2
##              head-hybridization=sp3
##
###################################

i=$((i+1))
TYPES[${i}]=" derivative head-hybridization=sp2 "
NAMES[${i}]=" ACX "
i=$((i+1))
TYPES[${i}]=" derivative head-hybridization=sp3 "
NAMES[${i}]=" MEX "
i=$((i+1))
TYPES[${i}]=" derivative head-hybridization=sp3 formal-charge=-1 "
NAMES[${i}]=" SO3 "

###################################
##  Monosaccharides
###################################
##
##  carbohydrate
##  monosaccharide
##              alpha
##              beta
##              D-isomer
##              L-isomer
##              pyranose
##              furanose
##              ketose
##              aldose
##              ulose
##              uronate
##              n-carbon=N
##                 ... where N is an integer corresponding to the
##                     number of carbons, e.g., n-carbon=3, n-carbon=8
##                     
##                     Leaving this as n-carbon so it will be general,
##                     applying also to lipid tails, etc.
##              modified
##                     has derivatives or other modification
##              unsaturated-uronate
##              deoxy
##              formal-charge=fc
##                     fc is the integer representing the formal charge
##              amine
##                     has an amine group that is not modified (sulfated, N-acetylated, etc), e.g. GlcN
##              N-acetyl
##                     There is one or more N-acetyl group that is
##                     normally part of this carb, so both GlcNAc and
##                     bacillosamine will match
##              N-sulfo
##              protonated
##              BFMP=aXb
##                     ...where aXb is the BFMP code for the ring pucker
##              non-standard-name
##
###################################

## These definitions are in a separate file that is generated by a script
. 3_carbohydrate_residue_tags.bash
