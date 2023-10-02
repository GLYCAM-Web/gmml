# GMML
The GLYCAM Molecular Modeling Library (GMML) is typically used as a library accessed by GEMS (GLYCAM Extensible Modeling Script).  

[Overview](#overview)

[Prerequisites](#prerequisites)

[Obtaining the software](#obtaining-the-software)

[Compiling the Library](#compiling-the-library)

[Testing the Library](#testing-the-library)

[Coding Standards](#coding-standards)

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

### Building GMML

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
* `libeigen3-dev` (Version >= `3.3.7`)

Installation instructions will vary according to what package manager your distro uses. If you are using apt as a package manager on a linux system, you should be able to use a command like this:

```bash
sudo apt-get update &&\
sudo apt-get install libssl1.1 libssl-dev git python3.9 python3.9-dev libboost-all-dev cmake g++ git-all libeigen3-dev
```
For other linux distros, please follow the instructions for the package managment software included with your system.

Please note that swig 4.0.2 must be installed from [their website](https://www.swig.org/download.html). Just kidding, most current linux distros already have `swig4.0` available right out of the box.

### Contributing to GMML

If you want to contribute to `gmml` you will also need to install the following packages:

* `clang-tidy-15`
* `shellcheck`
* `shfmt`

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

As of now, we rely upon an shell enviroment variable named `GEMSHOME` but this need is in the process of being removed. The code expects `GEMSHOME` to be the parent directory of where `gmml` is being built. An easy one liner, from within the `gmml` directory, is as follows:

```bash
user@host:.../gmml$ cd .. && export GEMSHOME=$(pwd) && cd -
```

Now to compile the library first make sure you are still in the `gmml` directory.
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

**DO NOT JUST FIRE THE `updateCmakeFileList.sh` SCRIPT AND NOT KNOW WHAT IS GOING ON. The method implemented is done in order to avoid a taxing typical cmake pattern; if the script is just fired off too many times we will have to remove it in order to avoid possible undefined behavior**. Please note that not only for cmake, but for all compilers, one should not just grab every file present and compile; these type of things must have some thought to them. The reason why one should never just glob files that one *thinks* are what one needs to compile is due to the huge increase in chances of introducing unknown behavior.

Basically treat this the same way as one treats using `git add --all` as bad practice due to priming the code base to have a bunch of random files (that should not be pushed) added to the repo; instead of being able to directly avoid `git add --all` and using `git add <YOUR_SPECIFIC_FILES>` instead, **YOU** must be the difference between that logic if you call the script check the git.

The `cmakeFileLists` directory contains the ouput from our `./updateCmakeFileList.sh` script. This script goes through and grabs all our files that we want to compile. There are 3 types:

* `cFileList.txt` - this contains all of our cpp files and where they be

* `hDirectoryList.txt` - this contains all of the directories that OUR source headers are. In the compiler this gets passed `-I` flag

* `externalHDirectoryList.txt` - this contains all the directoires of the EXTERNAL source code headers. For example, the eigen library. This is done so our `compile_commands.json` will use all these files with the `-isystem` flag which makes running tools much easier. 

---
## Testing the Library

From within the `gmml` directory, you must change your current working directory to the `gmml/tests` directory. Note that `<NUM_JOBS>` is however many tests you want to run at once.

```bash
gmml$ cd tests/
gmml/tests$ ./compile_run_tests.bash -j <NUM_JOBS>
```

Please note that running GMML bare metal will cause test 016 (svg drawing) to fail, this is due to not setting the created svgs correctly and will eventually be fixed but for now don't worry if `016.test.DrawGlycan.sh` fails while running on bare metal; if you are utilizing the dev enviroment all tests are expected to pass. This is of no concern because these tests need some extra things running to check, but those are internal for now.

The output will tell you whether or not the library is behaving appropriately and if all tests are passed the output will look similar to the following:

```bash
#### Beginning GMML tests ####
Number of tests found:  23
Number of testing jobs: 4

mkdir: created directory './tempTestOutputs'

Beginning test: ./000.test.buildBySequenceOldWay.sh
Beginning test: ./002.test.createAssemblyWritePDB.sh
Beginning test: ./003.test.SuperimpositionEigen.sh
Beginning test: ./004.test.PDBpreprocessor.sh

Testing 003.superimpositionEigen.cc... Test passed
Exit Code: 0

Beginning test: ./005.test.Overlaps.sh

Testing 002.create_Assembly_WritePDB.cc... Test passed
Exit Code: 0

Beginning test: ./006.test.BFMP-RingShapeCalculation.sh

Testing 000.buildBySequence.cc... Test passed
Exit Code: 0

Beginning test: ./007.test.DetectSugars.sh

Testing 006.ringShapeDetection.cc (BFMP Ring Shape Calculation)... Test passed.
Exit Code: 0

Beginning test: ./008.test.PDB2GlycamAndSubgraphMatching.sh

Testing 005.overlaps.cc... Test passed
Exit Code: 0

Beginning test: ./009.test.Reorder_and_Label_Sequence.sh

Testing 009.reorderSequence.cc (Sequence reordering and labeling)... Test passed.
Exit Code: 0

Beginning test: ./010.test.buildBySequenceRotamer.sh

Testing 008.pdb2glycam.cc and molecule subgraph matching... Test passed.
Exit Code: 0

Beginning test: ./011.test.writeResNumbers.sh

Testing 010.buildBySequenceRotamer.cc... Test passed.
Exit Code: 0

Beginning test: ./012.test.AddSolventNeutralize.sh

Testing 011.writeResNumbers.cc (write original and new residue numbers into a PDB file)... Test passed.
Exit Code: 0

Beginning test: ./014.test.SequenceParser.sh

Testing 014.test.SequenceParser.cc... Test passed.
Exit Code: 0

Beginning test: ./015.test.SequenceToAssembly.sh

Testing 012.AddSolventNeutralize... Test passed.
Exit Code: 0

Beginning test: ./016.test.DrawGlycan.sh

Testing PDBPreprocessor... ~15 seconds
Test passed for tests/inputs/004.preprocessorInput_4mbz.pdb.
Test passed for tests/inputs/004.preprocessorInput_original.pdb.
Test passed for tests/inputs/004.preprocessorInput_reduce.pdb.
Exit Code: 0

Beginning test: ./017.test.GlycoproteinBuilder.sh

Testing 016.test.DrawGlycan.cc...0.svg tests/correct_outputs/016.output_SVGs/0.svg differ: byte 132, line 2
Test FAILED! Output file 0.svg different to tests/correct_outputs/016.output_SVGs/0.svg
Exit Code: 1

Beginning test: ./018.test.GlycoproteinBuilderTable.sh

Testing 018.test.createGlycosylationTables.cpp... Test passed.
Exit Code: 0

Beginning test: ./019.test.newPDBClass.sh

Testing 007.detectSugars.cc... Test passed.
Exit Code: 0

Beginning test: ./020.test.parameterFiles.sh

Testing 015.test.SequenceAssembly.cc... Test passed.
Exit Code: 0

Beginning test: ./021.test.cdsSequence.sh

Testing 020.test.parameterFiles.cpp... Test passed.
Exit Code: 0

Beginning test: ./022.test.libraryFileReader.sh

Testing 022.test.libraryFileReader.cpp... Test passed.
Exit Code: 0

Beginning test: ./023.test.carbohydrateBuilder.sh

Testing 021.test.cdsSequence.cpp... Test passed.
Exit Code: 0

Beginning test: ./024.test.wiggleToSite.sh

Testing 024.wiggleToSite...Test passed.
Exit Code: 0

Testing 023.carbohydrateBuilder... Test passed.
Exit Code: 0

Testing 019.test.newPDBClass.cpp... ~30 seconds. Test passed.
Exit Code: 0

Testing 017.test.GlycoproteinBuilder.cpp... Test passed.
Exit Code: 0

######## GMML TESTS COMPLETED ########
Required tests: 23
Passed tests:   22
Failed tests:   1
Time taken:     34 seconds
######################################

!!! OUTPUT OF THE 1 GMML TEST(S) THAT FAILED !!!

Testing 016.test.DrawGlycan.cc...0.svg tests/correct_outputs/016.output_SVGs/0.svg differ: byte 132, line 2
Test FAILED! Output file 0.svg different to tests/correct_outputs/016.output_SVGs/0.svg
Exit Code: 1

!!! FINISHED PRINTING FAILED TESTS !!!
```

---
## Documentation

The official documentation for both GEMS and GMML can be found on the main GLYCAM website:

* GEMS - [http://glycam.org/gems](http://glycam.org/gems "GEMS")
* GMML - [http://glycam.org/gmml](http://glycam.org/gmml "GMML")

---
## Coding Standards

In order to make deving on the library consistent, we must enforce coding standards. They will be added piecewise, including the appropriate tests (be them pre-commit/push hooks, ci/cd hooks, etc.) and will be outlined below.

### Branch Naming

We use a gitflow based workflow, due to us not being able to think of all the structures that could break things we slowly move our changes up in our production branches where each one should be increasingly stable. When developing on `gmml` create feature branches off of `gmml-test`.

All branch names must take the form of `<branchType>_<descriptiveName>`. Be sure that you have a good descriptive name. The branch types are as follows:

- `feature` - any new feature that is being created. This will be the most commonly used branch type used.
- `bugfix` - branches that are fixing bugs
- `hotfix` - This is to fix critical bugs, etc. that is blocking people and that need to be handled immediately
- `playground` - For when you are not creating an actual feature or something that will be integrated into the main branches. Branches with this name would be for playing around with ideas.
- `juggle` - this is for branches that are intermediates between others, think if you have some horrible merge conflicts or something and want to create a to help make your life a bunch easier, that is what this branch is for. For instance, you could make the `juggle_zingusMerge` branch that could be made. This one will be used rarely, but it is nice to have in our back pocket. For instance lets say two devs are working on branches that need to be merged together before everything is placed into `gmml-test`, that is when this naming will be used.

Some examples of good branch names are:

- `feature_basicGithubActions`
- `bugfix_addKDNOparams`
- `playground_llvmTooling`

### Pre-Commit Hooks

We run various pre-commit hooks to help ensure `gmml`'s commit history remains as clean and readable as possible.

#### Hooks for C-Files

All code must follow the format described in the `.clang-format` file, and the pre-commit hook will ensure the commited format is correct. The precommit hook will ensure all files you want to commit are correctly formatted. Any files that are not correctly formatted will be listed in the terminal you tried to commit from, if you are using something like `gitflow` or `gitkraken` check the logs. Many code editors, IDEs or text editors, have the ability to apply a specific format on save of the file, so save yourself headaches and set that up.

Now, how do you format a specific file?

```bash
user@host:.../gmml$ clang-tidy-15 -i path/to/bad/file.cpp 
```

What if you did a bunch of files and want to be lazy? This can miss a couple bits that need to be changed so run it a couple times, it also will use all your cores but hey it is pretty quick.

```bash
user@host:.../gmml$ find . -not -path "./cmakeBuild/*" -type f -iname "*.cpp" -o -iname "*.hpp" -o -iname "*.h" -o -iname "*.cc" | xargs -P $(nproc --all --ignore=2)  -I % sh -c 'clang-format-15 -i %'
```

#### Hooks for Shell Scripts

In order to commit any shell scripts, the files must adhear to both our formatting and linting commands. Do not worry, linting shell scripts is very quick. Our checks are defined directly below:

* Formatting a shell script (NOTE: the formatting will be immediately applied to the script in question):
```bash
user@host:.../gmml$ shfmt -i 4 -ci -fn -w path/to/bad/script.sh
```

* Linting a shell script (NOTE: you will have to edit your script so it no longer has the issues the linter displays):
```bash
user@host:.../gmml$ shellcheck --enable=require-variable-braces,quote-safe-variables,add-default-case path/to/bad/script.sh
```

---
## Depreciated Instructions

The GLYCAM Molecular Modeling Library (GMML) was designed to be used as a library accessed by GEMS (GLYCAM Extensible Modeling Script), but can be used as a standalone library.

More information about GEMS can be found here:

Website:  [http://glycam.org/gems](http://glycam.org/gems "GLYCAM GEMS")  
Github:  [https://github.com/GLYCAM-Web/gems](https://github.com/GLYCAM-Web/gems "GLYCAM GEMS - Github")

To get started, follow the [Download and Install](http://glycam.org/docs/gems/download-and-install/ "Download and Install") instructions. These instructions will walk you through the steps to obtain and configure the software, and also test the installation.
  
To compile and use the programs that are based on gmml (e.g. the carbohydrate or glycoprotein builders) go to their subfolders (e.g. internalPrograms/GlycoproteinBuilder/) and follow the compilation instructions in the readme there.
