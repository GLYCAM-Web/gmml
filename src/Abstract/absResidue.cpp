#include "includes/Abstract/absResidue.hpp"
#include "includes/CodeUtils/biology.hpp"
#include "includes/CodeUtils/templatedSelections.hpp"
using Abstract::absResidue;

Abstract::ResidueType absResidue::determineType(const std::string &residueName) const
{
	if ( codeUtils::isElementPresent(biology::proteinResidueNames.begin(), biology::proteinResidueNames.end(), residueName) )
	{
		return ResidueType::Protein;
	}
	// Do we want to figure out solvent, aglycone etc here too? Maybe.
	return ResidueType::Undefined;
}
