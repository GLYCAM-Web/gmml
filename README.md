# GMML
The GLYCAM Molecular Modeling Library (GMML) is typically used as a library accessed by GEMS (GLYCAM Extensible Modeling Script).

[Overview](#overview)

[Prerequisites](#prerequisites)

[Obtaining the software](#obtaining-the-software)

[Compiling the Library](#compiling-the-library)

[Testing the Library](#testing-the-library)

[Documentation](#documentation)

---
## Overview

GMML provides a library for common molecular modeling tasks.  It is 
particularly well-tuned for models of carbohydrates and systems that
contain carbohydrates.

### Used by [GLYCAM-Web](https/glycam.org)

This code also serves as the main molecular modeling engine for GLYCAM-Web.  

### Funding Sources

We are very grateful to our funders.  
[Please check them out!](https://github.com/GLYCAM-Web/website/blob/master/funding.md)


## Prerequisites

In order to build GMML, you are required to have the following software available on your system:

* `libssl1.1`
* `libssl-dev`
* `cmake` (Version >= `3.13.4`)
* `boost`
* `g++` (Version >= `7.0`)
* `make`
* `python3.9` (Version `3.9.12`)
* `python3.9-dev` (Version `3.9.12`)
* `git`
* `swig` (Version `4.0.2`)

Installation instructions will vary according to what package manager your distro uses. If you are using apt as a package manager on a linux system, you should be able to use a command like this:

```bash
sudo apt-get update &&\
sudo apt-get install libssl1.1 libssl-dev git python3.9 python3.9-dev libboost-all-dev cmake g++ git-all
```
For other linux distros, please follow the instructions for the package managment software included with your system.

Please note that swig 4.0.2 must be installed from [their website](https://www.swig.org/download.html)

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
To control the number of processors used during the *`make`* process, use the `-j` flag for our `make.sh`, so to run with 8 cores we would run `./make.sh -j 8`.

Also, we have the option to wrap our code using `swig`. There are two methods to do this.

1. Once the makefile is generated using `cmake`, you can go into the `cmakeBuild` directory (or wherever you threw the makefile) and use the `gmml_wrapped` make target.

2. You can just call, from the base `GMML` directory, `./make.sh -w` 

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

### Updating file lists and whatnot

First off, **DO NOT JUST FIRE THE SCRIPT AND NOT KNOW WHAT IS GOING ON. THIS METHOD IS DONE SO I DONT HAVE TO FORCE PEOPLE TO USE A FULL CMAKE PATTERN, IF THESE GOOD GRACES ARE VIOLATED YOU WILL RUIN IT FOR EVERYONE BECAUSE I WILL BE FORCED TO REMOVE THIS WORKAROUND**. The reason I am making such a big deal about this is because no one should just auto grab what files we need to build and compile and run, we would have less knowledge about what our code is doing and it greatly increases the chances of introducing unknown behavior. So what is the usual method people use for cmake? Well, it is annoying and I am lazy so we do this workaround. Typically you have to have a `CMakeLists.txt` in every single directory that you add with the `add_subdirectory(<DIR HERE>)` cmake command. That is super annoying, but it has some merrits; mostly being that we know what is going on with out code. Other method that is **EXTREMELY FROWNED** is using cmake globbing (calling commands in the `CMakeLists.txt` file) to get all our `.cc/.cpp/.h/hpp` files. This is bad because its auto grabbing files to build and we should know whats going on with our code. So I did a middle of the ground method, we run a diff on the command thats used to grab the data thats in the file lists and on the file lists themselves. If the diff is different we update the lists, if not we dont update the lists. Thats about it. You can also pass a `-t` flag to the `updateCmakeFiles.sh` so it grabs all the test code files which we would want to do when we are running code analysis on our code base so we can lint the code and code. When you run the script, make sure you know whats going on and why each file is either being removed or added to file lists.

The `cmakeFileLists` directory contains the ouput from our `updateCmakeFiles.sh` script. This script goes through and grabs all our files that we want to compile. There are 3 types:

* `cFileList.txt` - this contains all of our cpp files and where they be

* `hDirectoryList.txt` - this contains all of the directories that OUR source headers are. In the compiler this gets passed `-I` flag

* `externalHDirectoryList.txt` - this contains all the directoires of the EXTERNAL source code headers. For example, the eigen library. This is done so our `compile_commands.json` will use all these files with the `-isystem` flag which makes running tools much easier. 

---
## Testing the Library

From within the `gmml` directory, you must change your current working directory to the `gmml/tests` directory.

```bash
gmml$ cd tests/
gmml/tests$ ./compile_run_tests.bash
```

Please note that running GMML bare metal will cause some tests to fail. This is of no concern because these tests need some extra things running to check, but those are internal for now. 

The output will tell you whether or not the library is behaving appropriately and if all tests are passed the output will look similar to the following:

```bash
Number of tests found: 18
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
Testing PDBPreprocessor... Test passed.

Using test file:  005.test.Overlaps.sh 
Testing Overlaps function... Test passed.

Using test file:  006.test.BFMP-RingShapeCalculation.sh 
Testing BFMP Ring Shape Calculation... Test passed.

Using test file:  007.test.DetectSugars.sh 
Testing detectSugars... Test passed.

Using test file:  008.test.PDB2GlycamAndSubgraphMatching.sh 
Testing pdb2glycam and molecule subgraph matching... Iupac name: DGalpb1-4DGlcpNAcb1-3DGalpb1-4DGlcpb1-ROH
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

Using test file:  014.test.SequenceParser.sh 
Testing 014.test.SequenceParser.cc... Test passed.

Using test file:  015.test.SequenceToAssembly.sh 
Testing 015.test.SequenceAssembly.cc... Test FAILED!. Output file different

Using test file:  016.test.DrawGlycan.sh 
Testing 016.test.DrawGlycan.cc...
ls: cannot access '*.svg': No such file or directory
Test FAILED!. Output file different

Using test file:  017.test.GlycoproteinBuilder.sh 
Testing 017.test.GlycoproteinBuilder.cpp... Test passed.

18 tests were attempted
16 tests passed 
18 were required
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
