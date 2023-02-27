//////////////////////////////////////////
// Instructions for making the ontology //
//////////////////////////////////////////

// You should be in the GRPC container of a running Development Environment
// If not, do that first.  

// If GMML wasn't compiled (IE you used bin/start.sh instead of Start-All-DevEnv.bash):

	cd $GEMSHOME/gmml/
	./make.sh

// Compile a modified detect_sugars with settings for building the database for Virtuoso
// See gmml/tests/tests/detect_sugars_RunArchive.cc for more

	cd tests/
  g++ -std=c++17 -I $GEMSHOME/gmml/ -L$GEMSHOME/gmml/bin/ -Wl,-rpath,$GEMSHOME/gmml/bin/ \
		tests/detect_sugars_RunArchive.cc -lgmml -pthread -o archiveRun_detect_sugars

// Run the files from the PDB

	cd runDetectSugars/scripts/

	nohup ./controlCPU.sh /programs/repos/PDB/all 4 &

// 4 'CPUs' used because the detect_sugars_runArchive.cc file uses 4 CPUs to build the structure
// and they don't build the structure the entire time (using 1 core for everything else)

// If you want to run a subset of PDBs, make a list of PDB IDs in a file called UpdatedPDBs.txt
// in the scripts directory (One PDB per line, case doesn't matter)
// and uncomment the "if grep ...UpdatedPDBs.txt; then" and the  "fi" statement atthe end of run_PDB
// in ontology.sh.

//Then, run as above with nohup
