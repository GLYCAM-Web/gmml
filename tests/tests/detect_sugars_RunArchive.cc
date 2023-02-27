// This is a C++ version of the GEMS test, detect_sugars.
// By: Davis Templeton
// Edited by: Dave Montgomery

#include <cstdlib>
#include <iostream>
#include <vector>
#include "../../includes/gmml.hpp"

// This looks horrible, but it works. Someone make it
// look better.
const std::string GEMSHOME_ERROR = "\n"
"Must set GEMSHOME environment variable.\n\n"
"    BASH:   export GEMSHOME=/path/to/gems\n"
"    SH:     setenv GEMSHOME /path/to/gems\n";

const std::string CCDHOME_ERROR = "\n"
"Must set CCDHOME environment variable. \n\n"
"    BASH:   export CCDHOME=/path/to/CCD/ligand-dict\n"
"    SH:     setenv CCDHOME /path/to/CCD/ligand-dict\n";

const std::string USAGE = "\n"
"Usage:\n\n"
"    detect_sugars PDB_file.pdb\n\n"
"The output goes to standard out (your terminal window, ususally).\n"
"So, alternately:\n\n"
"    detect_sugars PDB_file.pdb > output_file_name\n";

int main(int argc, char* argv[]) {
  // First get the GEMSHOME & CCDHOME environment variables
  char* gemshome_env_var = std::getenv("GEMSHOME");
  char* ccdhome_env_var = std::getenv("CCDHOME");

  // Check if the environment variables exist.
  if(!gemshome_env_var)
  {
    std::cout << GEMSHOME_ERROR << std::endl;
    return EXIT_FAILURE;
  }
  if(!ccdhome_env_var)
  {
    std::cout << CCDHOME_ERROR << std::endl;
    return EXIT_FAILURE;
  }

  // Check to make sure we have enough command line arguments.
  if(argc < 2)
  {
    std::cout << USAGE << std::endl;
    return EXIT_FAILURE;
  }
  std::string GEMSHOME(gemshome_env_var);
  std::string CCDHOME(ccdhome_env_var);

  // Get the Amino Lib file from GMML.
  std::vector<std::string> aminolibs;
  aminolibs.push_back(GEMSHOME + "/gmml/dat/CurrentParams/leaprc.ff12SB_2014-04-24/amino12.lib");

  // Get command line argument, which should be an PDB file.
  std::string pdb_file(argv[1]);

  // Initialize an Assembly from the PDB file and remove hydrogens until the bug 
  // (involving hydrogens forming multiple bonds.. I think) is fixed
  MolecularModeling::Assembly assembly(pdb_file, gmml::PDB);
  assembly.RemoveAllHydrogenAtoms();

  // Create a copy of the assembly to check for sugars without protein or DNA/RNA to speed up building
  // by distance.  
  
  MolecularModeling::Assembly assembly2(pdb_file, gmml::PDB);
  assembly2.RemoveAllProteinResidues();
  assembly2.RemoveAllNucleicAcidResidues();
  assembly2.RemoveAllHydrogenAtoms();
  assembly2.BuildStructureByDistance(3);
  std::vector< Glycan::Oligosaccharide* > oligos = assembly2.ExtractSugars(aminolibs, false, false, false);

  // If there are sugars in the assembly, build the assembly with protein and DNA/RNA and run through sugar detection
  if(oligos.size() > 0)
  {
    // assembly.RemoveAllNucleicAcidResidues();
    assembly.BuildStructureByDistance(3);

    // Find the Sugars.
    bool glyprobity_report = false;
    bool populate_ontology = true;
    bool individualOntologies = true;
    assembly.ExtractSugars(aminolibs, glyprobity_report, populate_ontology, individualOntologies, CCDHOME);
  }

  // YAY! We made it!
  return EXIT_SUCCESS;
}
