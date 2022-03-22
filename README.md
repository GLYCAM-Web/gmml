# GMML
The GLYCAM Molecular Modeling Library (GMML) is typically used as a library accessed by GEMS (GLYCAM Extensible Modeling Script).

[Prerequisites](#prerequisites)

[Obtaining the software](#obtaining-the-software)

[Compiling the Library](#compiling-the-library)

[Testing the Library](#testing-the-library)

[Documentation](#documentation)

---
## Prerequisites

In order to build GMML, you are required to have the following software available on your system:

* `libssl1.1`
* `libssl-dev`
* `cmake` (Version >= 3.13.4)
* `boost`
* `g++` (or a c++ compiler of your choice)
* `make`
* `python` (version >= 3.9.5)
* `python-dev` (version >= 3.9.5, possibly not needed)
* `git`

Installation instructions will vary according to what package manager your distro uses. If you are using apt as a package manager on a linux system, you should be able to use a command like this:

```bash
sudo apt-get update &&\
sudo apt-get install libssl1.1 libssl-dev git python3.9 python3.9-dev libboost-all-dev cmake g++ git-all
```
For other linux distros, please follow the instructions for the package managment software included with your system.

---
## Obtaining the software
The following does not require root access, but it does require one has `git` installed.

1. Navigate to the directory that you would like to have `gmml` live. Please note that in order to use the produced library with [`gems`](https://github.com/glycam-web/gems/) the `gmml` directory must be placed within the `gems` directory.

2. Clone `gmml` from the git repo and place into a folder named `gmml`
```bash
git clone https://github.com/GLYCAM-Web/gmml.git gmml
```

**NOTE:** There are non-git ways to obtain gems and gmml.  Please see the [`gems repo`](https://github.com/GLYCAM-Web/GEMS) and the [`gmml repo`](https://github.com/GLYCAM-Web/GMML) on github for more information.

---
## Compiling the Library

To compile the library first make sure you are still in the `gmml` directory.
```bash
pwd
./gmml
```
To control the number of processors used during the *`make`* process, use the `-j` flag for our `make.sh`, so to run with 8 cores we would run `./make.sh -j 8`

Now all one must do is run the make script.

```bash
$./make.sh
```

This will create the needed `cmake` files and will add the following directories within the `gmml` directory:

* `lib` (Contains: the `gmml` shared object libary, `libgmml.so`)
* `cmakeBuild` (Contains: all files produced by the `cmake` command, a `compile_commands.json` file to be used with tools of your choice, and all files contained within the directories listed above)

You can either use the `libgmml.so` file within the `lib` directory or the `libgmml.so` file within the `cmakeBuild` directory. They are the exact same.

Both the `build` and `lib` directories must remain untouched because `gems` utilizes them both and expects both to be in the state that `./make.sh` leaves them.

Please enter `./make.sh -h` for help regarding the make script.

---
## Testing the Library

From within the `gmml` directory, you must change your current working directory to the `gmml/tests` directory.

```bash
gmml$ cd tests/
gmml/tests$ ./compile_run_tests.bash
```

The output will tell you whether or not the library is behaving appropriately and if all tests are passed the output will look similar to the following:

```bash
gmml/tests$ ./compile_run_tests.bash

Number of tests found: 14
Beginning testing.


Using test file:  000.test.buildBySequenceOldWay.sh
Testing buildBySequence... Test passed.

Using test file:  001.test.buildBySequenceMetaWay.sh
Testing buildBySequenceMeta... Test passed.

Using test file:  002.test.createAssemblyWritePDB.sh
Testing create_Assembly_WritePDB... Test passed.

Using test file:  003.test.SuperimpositionEigen.sh
Testing superimposition_Eigen... Test passed.

Using test file:  004.test.PDBpreprocessor.sh
Testing PDBpreprocessor... Test passed.

Using test file:  005.test.Overlaps.sh
Testing Overlaps function... Test passed.

Using test file:  006.test.BFMP-RingShapeCalculation.sh
Testing BFMP Ring Shape Calculation... Test passed.

Using test file:  007.test.DetectSugars.sh
Testing detectSugars... Test passed.

Using test file:  008.test.PDB2GlycamAndSubgraphMatching.sh
Testing pdb2glycam and molecule subgraph matching... 1 matches found.
Test passed.

Using test file:  009.test.Reorder_and_Label_Sequence.sh
Testing Sequence reordering and labeling... Test passed.

Using test file:  010.test.buildBySequenceRotamer.sh
Testing buildBySequenceRotamer... Test passed.

Using test file:  011.test.writeResNumbers.sh
Testing writing original and new residue numbers into a PDB file... Test passed.

Using test file:  012.test.AddSolventNeutralize.sh
Testing 012.AddSolventNeutralize... Test passed.

Using test file:  013.test.buildOligoaccharideLibrary.sh
Testing buildOligosaccharide library... Test passed.

14 tests were attempted
14 tests passed
14 were required
```
---
## Documentation

The official documentation for both GEMS and GMML can be found on the main GLYCAM website:

* GEMS - [http://glycam.org/gems](http://glycam.org/gems "GEMS")
* GMML - [http://glycam.org/gmml](http://glycam.org/gmml "GMML")

---
## Depreciated Instructions

The GLYCAM Molecular Modeling Library (GMML) was designed to be used as a library accessed by GEMS (GLYCAM Extensible Modeling Script), but can be used as a standalone library.

More information about GEMS can be found here:

Website:  [http://glycam.org/gems](http://glycam.org/gems "GLYCAM GEMS")  
Github:  [https://github.com/GLYCAM-Web/gems](https://github.com/GLYCAM-Web/gems "GLYCAM GEMS - Github")

To get started, follow the [Download and Install](http://glycam.org/docs/gems/download-and-install/ "Download and Install") instructions. These instructions will walk you through the steps to obtain and configure the software, and also test the installation.

To compile and use the programs that are based on gmml (e.g. the carbohydrate or glycoprotein builders) go to their subfolders (e.g. internalPrograms/GlycoproteinBuilder/) and follow the compilation instructions in the readme there.
