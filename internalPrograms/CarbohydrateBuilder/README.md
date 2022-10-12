# CarbohydrateBuilder
Uses GMML to create 3D structures of carbohydrates based on the GLYCAM forcefield.

### Notes
Project is under development, contact olivercgrant "at" gmail.com with queries. 
The underlying code is used by the builder currently available on glycam.org/cb.
Has been tested on linux, but should install on both Mac and Windows with appropriate C++ complilers.

### Prerequisites
You'll need GMML. See here for installation instructions: http://glycam.org/docs/gems/download-and-install/.

### Installation
export GEMSHOME=<Your Path To Gems > # eg: export GEMSHOME=/home/oliver/Programs/gems

#### Makefile:
Some commands defined in the Makefile are:
* $ make
* $ make all
* $ make clean

    So you can just type "make" to compile the program.
    
#### Comand line (skip the Makefile):
g++ -std=c++17 -I $GEMSHOME/gmml/includes/ -I includes/ -L$GEMSHOME/gmml/bin/ -Wl,-rpath,$GEMSHOME/gmml/bin/ src/*.cpp -lgmml -o bin/carbBuilder

### Testing
Once compiled, you can run:
./bin/carbBuilder
to get a usage statement and example inputs. 
You must create an output folder e.g:
mkdir outputs/

### Known issues
Very large (40-ish residues), branched oligosaccharides will take seemingly forever to resolve clashes. This will be fixed in upcoming releases.
