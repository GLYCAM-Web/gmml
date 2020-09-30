#include "../../includes/gmml.hpp"
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <map>

int main(int argc, char* argv[])
{

    MolecularModeling::Assembly assemblyA (argv[1], gmml::InputFileType::PDB);
    assemblyA.BuildStructureByDistance();
    PdbFileSpace::PdbFile *outputPdbFile = assemblyA.BuildPdbFileStructureFromAssembly();
    outputPdbFile->Write("011.output.original.pdb");
    outputPdbFile = assemblyA.BuildPdbFileStructureFromAssembly(-1,0,-1,false); // Change from default of setting useResidueNumbers true
    outputPdbFile->Write("011.output.newNumbers.pdb");
}
