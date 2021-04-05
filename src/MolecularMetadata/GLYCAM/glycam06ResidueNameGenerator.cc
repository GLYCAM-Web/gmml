#include "glycam06ResidueNameGenerator.hpp"
#include "glycam06LinkageCodes.hpp"
#include "glycam06residuecodes.hpp"
#include <iostream> // for cout, can remove after debug
#include <locale> // for isLower()
#include <sstream> // for string stream


namespace gmml
{
namespace MolecularMetadata
{
namespace GLYCAM
{
std::string Glycam06ResidueNameGenerator(std::string linkages, char isomer, std::string inputResName, char ringType, std::string residueModifier, char configuration)
{
/* 	Example inputs: 
	linkages: "2,3" , "1" , "Terminal" , "4,7"
	isomer: Always D or L
	inputResName: "Glc" , "Neu" , "Ido"
	ringType: Always f or p
	residueModifier: "NAc" , "5Ac" , "A" 
	configuration: Always a or b

	Example output:
	UYB, 0GA etc. See glycam naming / nomenclature.
*/


	// Link code e.g. 0, 1, 2, W, Z etc
	std::cout << "\nInputs:\nlinkages: " << linkages << "\nisomer: " << isomer << "\ninputResName: " << inputResName << "\nringType: " << ringType << "\nresidueModifier: " << residueModifier << "\nconfiguration: " << configuration << std::endl;
	std::string linkCode = "";
	if(!linkages.empty())
	{
		Glycam06LinkageCodesLookupContainer linkageCodeLookup;
		linkCode = linkageCodeLookup.GetCodeForLinkages(linkages);
		if (linkCode.empty())
		{
			auto message = "No linkage code found in GMML metadata for linkage: " + linkages;
			throw message;
		}
	}

	// Configuration Code i.e. A/B/U/D
	std::string configurationCode(1, configuration); // Convert to string with string constuctor. 1 copy.
	if ((ringType == 'f') && (configuration == 'a'))
	{
		configurationCode = "D";
	}	
	else if ((ringType == 'f') && (configuration == 'b'))
	{
		configurationCode = "U";
	}
	else if ((ringType == 'p') && (configuration == 'a'))
	{
		configurationCode = "A";
	}
	else if ((ringType == 'p') && (configuration == 'b'))
	{
		configurationCode = "B";
	}

	// Residue Code e.g. G, U, A, KN, ROH, NLN
	Glycam06ResidueNamesToCodesLookupContainer ResidueCodeLookup;
	auto residueCode = ResidueCodeLookup.GetCodeForResidue(inputResName + residueModifier);
	if (residueCode.empty())
	{
		residueCode = ResidueCodeLookup.GetCodeForResidue(inputResName + ringType + residueModifier + configuration);
	}
	if (residueCode.empty())
	{
		residueCode = ResidueCodeLookup.GetCodeForResidue(isomer + inputResName + ringType + residueModifier + configuration);
	}
	if (residueCode.empty())
	{
		residueCode = ResidueCodeLookup.GetCodeForResidue(isomer + inputResName + ringType + residueModifier);
	}
	if (residueCode.empty())
	{
		auto message = "No residue code found in GMML metadata for residue: " + isomer + inputResName + ringType + residueModifier + configuration;
		throw message;
	}

	if (residueCode.size() > 1)
	{
		configurationCode = ""; // Make it empty, it will implied by residueCode
	}

	// D vs L sugars. Residue code will be lowercase for L sugars
	if ((isomer == 'L') && (residueCode.size() == 1))
	{
		residueCode = std::tolower(residueCode.at(0));
	}

	// ConfigurationCode may be empty.
	std::cout << "Returning: " << (linkCode + residueCode + configurationCode) << std::endl;
	return (linkCode + residueCode + configurationCode);
}
} // close namespace
} // close namespace
} // close namespace