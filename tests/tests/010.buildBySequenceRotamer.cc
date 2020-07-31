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
    std::string condensed_sequence = "DManp[2S,3Me]a1-6DManpa1-6[DGlcpNAcb1-4][DNeu5Aca2-6DGalpb1-4DGlcpNAc[3S]b1-2DManpa1-3]DManpb1-4DGlcpNAc[6Me]b1-4DGlcpNAcb1-OH";
    CondensedSequenceSpace::carbohydrateBuilder carbBuilder(condensed_sequence, prep);
    //std::cout << carbBuilder.GenerateUserOptionsJSON() << std::endl;
    CondensedSequenceSpace::singleRotamerInfoVector rotamerInfoVector =
	{ 
		{ "10" , "", "phi", "-g", ""},
		{ "10" , "", "omg", "gt", ""},
        { "11" , "", "omg", "gg", ""},
		{ "6" , "", "omg" , "gg", ""},
	};
    std::string fileOutputDirectory = "unspecified"; // If it's "unspecified" or you don't pass a directory string it will write to the current working directory. 
	carbBuilder.GenerateRotamer(rotamerInfoVector, fileOutputDirectory);
}
//prep file

