#include "includes/Abstract/absResidue.hpp"
#include "includes/common.hpp"
#include "includes/CodeUtils/templatedSelections.hpp"
using Abstract::absResidue;

absResidue::Type absResidue::determineType(const std::string &residueName) const
{
	if ( codeUtils::isThingPresentInContainer(gmml::proteinResidueNames.begin(), gmml::proteinResidueNames.end(), residueName) )
	{
		return absResidue::Type::Protein;
	}
	// Do we want to figure out solvent, aglycone etc here too? Maybe.
	return absResidue::Type::Undefined;
}
