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
    //std::string condensed_sequence = "DManp[2S,3Me]a1-6DManpa1-6[DGlcpNAcb1-4][DNeup5Aca2-6DGalpb1-4DGlcpNAc[3S]b1-2DManpa1-3]DManpb1-4DGlcpNAc[6Me]b1-4DGlcpNAcb1-OH";
    std::string condensed_sequence = "DNeup5Aca2-3DGalpb1-OH";
    //std::string condensed_sequence = "DNeup5Aca2-8DNeup5Aca2-3DGlcpb1-OH";
    //std::string condensed_sequence = "DNeup5Aca2-8DNeup5Aca2-8DNeup5Aca2-8DNeup5Aca2-8DNeup5Aca2-OH";
    //std::string condensed_sequence = "DGlcpNAca1-NLN";
    //std::string condensed_sequence = "DNeup5Aca2-7[DNeup5Aca2-8]DNeup5Aca2-OH";
    CondensedSequence::carbohydrateBuilder carbBuilder(condensed_sequence, prep);
    carbBuilder.Print();
    std::string fileOutputDirectory = "unspecified"; // If it's "unspecified" or you don't pass a directory string it will write to the current working directory. 
    bool likelyShapesOnly = true; // You can just pass "true" to the function. Not sure I like this. Two functions probably more readable.
    std::cout << "Number of likely shapes for this sequence is " << carbBuilder.GetNumberOfShapes(likelyShapesOnly) << "\n";
    // Default is to calculate all possible.
    std::cout << "Number of possible shapes for this sequence is " << carbBuilder.GetNumberOfShapes() << "\n";

    for(auto &linkageInfo : carbBuilder.GenerateUserOptionsDataStruct())
    {
        std::cout   << "Name: " << linkageInfo.linkageName_ << ", LinkageIndex: " << linkageInfo.indexOrderedLabel_ << ", Res1: "
                << linkageInfo.firstResidueNumber_ << ", Res2: " << linkageInfo.secondResidueNumber_ << "\n";
        for(auto &dihedralInfo : linkageInfo.likelyRotamers_)
        {
            std::cout << "    LikelyRotamers: " << dihedralInfo.dihedralName_ ;
            for(auto &rotamer : dihedralInfo.rotamers_)
            {
                std::cout << ", " << rotamer;
                // use when checking multiple rotamers, but only works for one linkage, and you must uncomment out the //fileName += "_" + rotamerInfo.selectedRotamer in GenerateSpecific3DStructure
//                CondensedSequence::SingleRotamerInfoVector rotamerInfoVector = { {linkageInfo.indexOrderedLabel_, linkageInfo.indexOrderedLabel_, dihedralInfo.dihedralName_, rotamer, ""} };
//                carbBuilder.GenerateSpecific3DStructure(rotamerInfoVector, "outputs/");
            }
            std::cout << "\n";
        } 
    }
    carbBuilder.GenerateSingle3DStructureDefaultFiles();
//    carbBuilder.GenerateUpToNRotamers(32);
    std::cout << "Fin\n";
}
//prep file

