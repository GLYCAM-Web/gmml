# GlycoProteinBuilder
Uses GEMS/GMML to add and adapt 3D structures of N-glycans and O-glycans onto glycoproteins. It can do this for Asn, Ser, Thr and Tyr.

## Schematic
![schematic](schematic/schematic.png)

### Notes
Project is under development, contact olivercgrant "at" gmail.com with queries. 
This code will replace the glycoprotein builder currently available on glycam.org/gp.
Has been tested on linux, but should install on both Mac and Windows with appropriate C++ complilers.

### Prerequisites

You'll need GMML. See here for installation instructions: http://glycam.org/docs/gems/download-and-install/.

### Installation of GlycoProteinBuilder
export GEMSHOME=<Your Path To Gems > # eg: export GEMSHOME=/home/oliver/Programs/gems

#### Makefile:
Some commands defined in the Makefile are:
* $ make
* $ make all
* $ make bin/gp_builder
* $ make build/<file_name>.o
* $ make clean

    So you can type "make" to compile the program.
    
#### Comand line (skip the Makefile):
g++ -std=c++17 -I $GEMSHOME/gmml/includes/ -I includes/ -L$GEMSHOME/gmml/bin/ -Wl,-rpath,$GEMSHOME/gmml/bin/ src/*.cpp -lgmml -o gp_builder

### Setup
Edit or create an input.txt file and place in a folder called tests/. See tests/simple/input.txt for an example.

You must provide:

    a protein 3D structure

    glycan 3D structure(s) or sequences in glycam condensed nomenclature (just like the carb builder here: glycam.org/cb)

    input.txt, which contains:

        protein file name

        the protein residue numbers you want to attach to (no automatic detection of sequons)

        the glycan you want to attach. Either the name in a library of PDB files or the glycan sequence
    
        example here: https://github.com/gitoliver/GlycoProteinBuilder/blob/stable/tests/tough/input.txt

### Bead based overlap calculation
In order to speed up the overlap calculation, certain atoms in the glycan and protein are replaced with large spheres that encompass the neighbouring atoms. The overlap calculation only looks at the beads, and as there are much fewer of them, it will be faster. The downside is that it is not as accurate and may be unncessarily optimizing. However, as the beads will encompass all atoms in the protein/glycan, if the bead overlap reaches zero, the per atom overlap will be zero.
Here is a figure showing the atoms being replaced by beads:
![bead replacment](schematic/beads.png)

## Algorithm
This is regularly being tweaked. Right now it's a combination of random movements and "wiggling" (grid searching within the known ranges of the glycosidic bonds).

### Known problems
It's perfect. Except when you mess up the inputs, then it segfaults, but that's on you. It's also slow (minutes to hours) and doesn't adjust any part of the protein except the Asn when trying to get the glycan to fit. 

