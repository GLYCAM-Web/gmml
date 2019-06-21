// This is a C++ version of the GEMS test, detect_sugars.
// By: Davis Templeton
#include <cstdlib>
#include <iostream>
#include <vector>
#include "gmml.hpp"

// This looks horrible, but it works. Someone make it
// look better.
const std::string GEMSHOME_ERROR = "\n"
"Must set GEMSHOME environment variable.\n\n"
"    BASH:   export GEMSHOME=/path/to/gems\n"
"    SH:     setenv GEMSHOME /path/to/gems\n";

const std::string USAGE = "\n"
"Usage:\n\n"
"    detect_sugars PDB_file.pdb\n\n"
"The output goes to standard out (your terminal window, ususally).\n"
"So, alternately:\n\n"
"    detect_sugars PDB_file.pdb > output_file_name\n";

int main(int argc, char* argv[]) {
    // First get the GEMSHOME environment variable
    char* gemshome_env_var = std::getenv("GEMSHOME");
    std::string GEMSHOME(gemshome_env_var);

    // Check if the environment variable exists.
    if(GEMSHOME == "") {
        std::cout << GEMSHOME_ERROR << std::endl;
        return EXIT_FAILURE;
    }

    // Check to make sure we have enough command line arguments.
    if(argc < 2) {
        std::cout << USAGE << std::endl;
        return EXIT_FAILURE;
    }

    // Get the Amino Lib file from GMML.
    std::vector<std::string> aminolibs;
    aminolibs.push_back(GEMSHOME + "/gmml/dat/CurrentParams/leaprc.ff12SB_2014-04-24/amino12.lib");

    // Get command line argument, which should be an PDB file.
    std::string pdb_file(argv[1]);
    // Initialize an Assembly from the PDB file.
    MolecularModeling::Assembly assembly(pdb_file, gmml::PDB);
    // Build by Distance
    assembly.BuildStructureByDistance(10);
    // Find the Sugars.
    assembly.ExtractSugars(aminolibs, false, true);
    //Note that to have individual ontology (.ttl) files or to have CCD lookup, you must provide 
    //a bool (true) for individual ontologies, and the path to the CCD which right now is just in my home directory
    
    // YAY! We made it!
    return EXIT_SUCCESS;
}
