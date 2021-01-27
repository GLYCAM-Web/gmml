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
    //std::string condensed_sequence = "DManp[2S,3Me]a1-6DManpa1-6[DGlcpNAcb1-4][DNeu5Aca2-6DGalpb1-4DGlcpNAc[3S]b1-2DManpa1-3]DManpb1-4DGlcpNAc[6Me]b1-4DGlcpNAcb1-OH";
    //std::string condensed_sequence = "DNeup5Aca2-8DNeup5Aca2-3DGlcpb1-OH";
    std::string condensed_sequence = "DNeup5Aca2-8DNeup5Aca2-8DNeup5Aca2-8DNeup5Aca2-8DNeup5Aca2-OH";
    //std::string condensed_sequence = "DGlcpNAca1-NLN";

    // std::string condensed_sequence = "DNeup5Aca2-7[DNeup5Aca2-8]DNeup5Aca2-OH";
    CondensedSequenceSpace::carbohydrateBuilder carbBuilder(condensed_sequence, prep);
    carbBuilder.Print();
    //std::cout << carbBuilder.GenerateUserOptionsJSON() << std::endl;
 //    CondensedSequenceSpace::SingleRotamerInfoVector rotamerInfoVector =
	// { 
	// 	{ "7" , "", "Phi", "t", ""},
	// 	{ "7" , "", "Omg", "gt", ""},
 //        { "9" , "", "Omg", "gg", ""},
	// 	{ "10" , "", "Omg" , "gt", ""},
	// };
    std::string fileOutputDirectory = "unspecified"; // If it's "unspecified" or you don't pass a directory string it will write to the current working directory. 
	//carbBuilder.GenerateSpecific3DStructure(rotamerInfoVector, fileOutputDirectory);
    bool likelyShapesOnly = true; // You can just pass "true" to the function. Not sure I like this. Two functions probably more readable.
    std::cout << "Number of likely shapes for this sequence is " << carbBuilder.GetNumberOfShapes(likelyShapesOnly) << "\n";
    // Default is to calculate all possible.
    std::cout << "Number of possible shapes for this sequence is " << carbBuilder.GetNumberOfShapes() << "\n";
    for(auto &linkageInfo : carbBuilder.GenerateUserOptionsDataStruct())
    {
        std::cout   << "Options: Name: " << linkageInfo.linkageName_ << ", Index: " << linkageInfo.indexOrderedLabel_ << ", Res1:" 
                << linkageInfo.firstResidueNumber_ << ", Res2:" << linkageInfo.secondResidueNumber_ << "\n";
        for(auto &dihedralInfo : linkageInfo.likelyRotamers_)
        {
            std::cout << "Likely: " << dihedralInfo.dihedralName_ ;
            for(auto &rotamer : dihedralInfo.rotamers_)
            {
                std::cout << ", " << rotamer;
            }
            std::cout << "\n";
        } 
    }
    carbBuilder.GenerateSingle3DStructureDefaultFiles();
}
//prep file

