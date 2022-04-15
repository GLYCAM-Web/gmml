#include "includes/gmml.hpp"
#include "includes/MolecularModeling/assembly.hpp"
#include "includes/ParameterSet/PrepFileSpace/prepfile.hpp"
#include "includes/ParameterSet/PrepFileSpace/prepfileresidue.hpp"
#include "includes/ParameterSet/PrepFileSpace/prepfileprocessingexception.hpp"
#include "includes/ParameterSet/OffFileSpace/offfile.hpp"
#include "includes/ParameterSet/OffFileSpace/offfileresidue.hpp"
#include "includes/ParameterSet/OffFileSpace/offfileprocessingexception.hpp"
#include <iostream>
#include <string>

int main()
{
    std::string prep = "../dat/prep/GLYCAM_06j-1.prep";
    PrepFileSpace::PrepFile* prepA = new PrepFileSpace::PrepFile(prep);
    std::string condensed_sequence = "DManp[2S,3Me]a1-6DManpa1-6[DGlcpNAcb1-4][DNeup5Aca2-6DGalpb1-4DGlcpNAc[3S]b1-2DManpa1-3]DManpb1-4DGlcpNAc[6Me]b1-4DGlcpNAcb1-OH";
    MolecularModeling::Assembly assemblyA = MolecularModeling::Assembly();
    assemblyA.SetName("CONDENSEDSEQUENCE");
    assemblyA.BuildAssemblyFromCondensedSequence (condensed_sequence, prepA);
    PdbFileSpace::PdbFile *outputPdbFile = assemblyA.BuildPdbFileStructureFromAssembly();
    outputPdbFile->Write("buildBySequence.pdb");
}
