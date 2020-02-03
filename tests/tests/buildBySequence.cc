#include "../../gmml/includes/gmml.hpp"
#include "../../gmml/includes/MolecularModeling/assembly.hpp"
#include "../../gmml/includes/ParameterSet/PrepFileSpace/prepfile.hpp"
#include "../../gmml/includes/ParameterSet/PrepFileSpace/prepfileresidue.hpp"
#include "../../gmml/includes/ParameterSet/PrepFileSpace/prepfileprocessingexception.hpp"
#include "../../gmml/includes/ParameterSet/OffFileSpace/offfile.hpp"
#include "../../gmml/includes/ParameterSet/OffFileSpace/offfileresidue.hpp"
#include "../../gmml/includes/ParameterSet/OffFileSpace/offfileprocessingexception.hpp"
#include <iostream>
#include <string>

int main()
{
    std::string prep = "../dat/prep/GLYCAM_06j-1.prep";
    PrepFileSpace::PrepFile* prepA = new PrepFileSpace::PrepFile(prep);
    std::string condensed_sequence = "DManp[2S,3Me]a1-6DManpa1-6[DGlcpNAcb1-4][DNeu5Aca2-6DGalpb1-4DGlcpNAc[3S]b1-2DManpa1-3]DManpb1-4DGlcpNAc[6Me]b1-4DGlcpNAcb1-OH";
    MolecularModeling::Assembly assemblyA = MolecularModeling::Assembly();
    assemblyA.SetName("CONDENSEDSEQUENCE");
    assemblyA.BuildAssemblyFromCondensedSequence (condensed_sequence, prepA);
    PdbFileSpace::PdbFile *outputPdbFile = assemblyA.BuildPdbFileStructureFromAssembly();
    outputPdbFile->Write("buildBySequence.pdb");
    //std::cout << "Done writing pdb." << std::endl;
    
    // Testing Yao's reordering thing
//    std::string sequence = "DManp[2S,3Me]a1-6DManpa1-6[DGlcpNAcb1-4][DNeu5Aca2-6DGalpb1-4DGlcpNAc[3S]b1-2DManpa1-3]DManpb1-4DGlcpNAc[6Me]b1-4DGlcpNAcb1-OH";
//    std::string sequence = "DManpa1-3[DGalpb1-4DGalpb1-4DGalpb1-4DGalpb1-4]LRhapa1-OH";
//    CondensedSequenceSpace::CondensedSequence cond_seq(sequence);
//    std::cout << "Output longest is: " << cond_seq.BuildLabeledCondensedSequence(CondensedSequenceSpace::CondensedSequence::Reordering_Approach::LONGEST_CHAIN, false) << std::endl;
//    std::cout << "Output sequence is: " << cond_seq.BuildLabeledCondensedSequence(CondensedSequenceSpace::CondensedSequence::Reordering_Approach::LOWEST_INDEX, false) << std::endl;

}
