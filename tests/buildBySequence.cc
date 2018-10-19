#include "../includes/gmml.hpp"
#include "../includes/MolecularModeling/assembly.hpp"
#include "../includes/ParameterSet/PrepFileSpace/prepfile.hpp"
#include "../includes/ParameterSet/PrepFileSpace/prepfileresidue.hpp"
#include "../includes/ParameterSet/PrepFileSpace/prepfileprocessingexception.hpp"
#include "../includes/ParameterSet/OffFileSpace/offfile.hpp"
#include "../includes/ParameterSet/OffFileSpace/offfileresidue.hpp"
#include "../includes/ParameterSet/OffFileSpace/offfileprocessingexception.hpp"
#include <iostream>
#include <string>

int main()
{
//prep
    std::string prep = "../dat/prep/GLYCAM_06j-1.prep";
    PrepFileSpace::PrepFile* prepA = new PrepFileSpace::PrepFile(prep);
    std::string condensed_sequence = "DManpb1-6[DManpa1-4]DManpa1-OH";
    MolecularModeling::Assembly assemblyA = MolecularModeling::Assembly();
    assemblyA.SetName("CONDENSEDSEQUENCE");
    assemblyA.BuildAssemblyFromCondensedSequence (condensed_sequence, prepA);
    PdbFileSpace::PdbFile *outputPdbFile = assemblyA.BuildPdbFileStructureFromAssembly();
    outputPdbFile->Write("buildBySequence.pdb");
    std::cout << "Done writing pdb." << std::endl;
}
//prep file

