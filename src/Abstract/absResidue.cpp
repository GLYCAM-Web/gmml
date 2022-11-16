#include "includes/Abstract/absResidue.hpp"
#include "includes/CodeUtils/biology.hpp"
#include "includes/CodeUtils/templatedSelections.hpp"
using Abstract::absResidue;

Abstract::ResidueType absResidue::determineType(const std::string &residueName) const
{
    if ( std::find(biology::proteinResidueNames.begin(), biology::proteinResidueNames.end(), residueName)  != biology::proteinResidueNames.end() )
    {
        return ResidueType::Protein;
    }
    // ToDo we want to figure out solvent, aglycone etc here too?.
    return ResidueType::Undefined;
}
