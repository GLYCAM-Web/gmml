// This is a test for the PDB API for detect_sugars.
// By: Davis Templeton
// Modified by: Dave Montgomery
#include <cstdlib>
#include <iostream>
#include <vector>
#include "../../includes/gmml.hpp"

std::string TextFileToString(std::string filename)
{
   std::string line; 
   std::stringstream dosString;
   std::ifstream inFile;
   inFile.open (filename.c_str());
       while(getline(inFile, line))
        {           
           dosString<<line<<"\r\n";
        }
    return dosString.str();
}

const std::string GEMSHOME_ERROR = "\n"
"Must set GEMSHOME environment variable.\n\n"
"    BASH:   export GEMSHOME=/path/to/gems\n"
"    SH:     setenv GEMSHOME /path/to/gems\n";

const std::string CCDHOME_ERROR = "\n"
"Must set CCDHOME environment variable.\n\n"
"    BASH:   export CCDHOME=/path/to/CCD/ligand-dict-v3\n"
"    SH:     setenv CCDHOME /path/to/CCD/ligand-dict-v3\n";

const std::string USAGE = "\n" // @suppress("Type cannot be resolved")
"Usage:\n\n"
"    detect_sugars PDB_file.pdb\n\n"
"The output goes to standard out (your terminal window, ususally).\n"
"So, alternately:\n\n"
"    detect_sugars PDB_file.pdb > output_file_name\n";

int main(int argc, char* argv[]) {
    // First get the GEMSHOME & CCD_HOME environment variables
    char* gemshome_env_var = std::getenv("GEMSHOME");
    std::string GEMSHOME(gemshome_env_var);
    
    char* ccdhome_env_var = std::getenv("CCDHOME");
    std::string CCDHOME(ccdhome_env_var);

    // Check if the environment variables exists.
    if(GEMSHOME == "") 
    {
        std::cout << GEMSHOME_ERROR << std::endl;
        return EXIT_FAILURE;
    }
    if(CCDHOME == "") 
    {
        std::cout << CCDHOME_ERROR << std::endl;
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

    // Get command line argument, which should be a text file with only atom cards.
    std::string pdb_file(argv[1]);
    
    // Convert file into stringstream
    std::string atomString = TextFileToString(pdb_file);
    std::stringstream atomStream;
    atomStream << atomString;
    
    // Initialize an Assembly from the stringstream.
    MolecularModeling::Assembly assembly(atomStream);
    // Build by Distance
    assembly.BuildStructureByDistance(10);
    // Find the Sugars.
    std::vector<MolecularModeling::Assembly::gmml_api_output> output = assembly.PDBExtractSugars(aminolibs, CCDHOME);
    for(std::vector<MolecularModeling::Assembly::gmml_api_output>::iterator it = output.begin(); it != output.end(); it++)
    {
      MolecularModeling::Assembly::gmml_api_output thisOutput = (*it);
      std::cout << "\nLinear Descriptors\n";
      for(int i = 0; i < thisOutput.linear_descriptors.size(); i++)
      {
        std::cout << thisOutput.linear_descriptors[i].first << "\n" << thisOutput.linear_descriptors[i].second << "\n";
      }
      std::cout << "\nIndices\n";
      for(int i = 0; i < thisOutput.indices.size(); i++)
      {
        std::cout << thisOutput.indices[i].first << ": " << thisOutput.indices[i].second << "\n";
      }
      std::cout << "\nResidue Links\n";
      for(int i = 0; i < thisOutput.residue_links.size(); i++)
      {
        for(int j = 0; j < thisOutput.residue_links[i].size(); j++)
        {
          if(j == thisOutput.residue_links[i].size()-1)
          {
            std::cout << thisOutput.residue_links[i][j] << "\n";
          }
          else
          {
            std::cout << thisOutput.residue_links[i][j] << ", ";
          }
        }
      }
      std::cout << "\nErrors and Warnings\n";
      for(int i = 0; i < thisOutput.error_warning_messages.size(); i++)
      {
        std::cout << thisOutput.error_warning_messages[i] << "\n";
      }
    }  
    // YAY! We made it!
    return EXIT_SUCCESS;
}
