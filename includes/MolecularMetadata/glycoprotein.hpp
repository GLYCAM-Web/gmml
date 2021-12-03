#ifndef GMML_INCLUDES_MOLECULARMETADATA_GLYCOPROTEIN_HPP
#define GMML_INCLUDES_MOLECULARMETADATA_GLYCOPROTEIN_HPP

#include <string>
#include <vector>
#include "includes/MolecularModeling/residue.hpp"

namespace glycoproteinMetadata
{
std::string LookupCodeForAminoAcidName(const std::string queryName);
std::string LookupLinkTypeForAminoAcidName(const std::string queryName);
std::vector<std::string> GetTagsForSequence(const std::string firstRes = "", const std::string midRes = "", const std::string lastRes = "");
std::string GetSequenceContextAndDetermineTags(MolecularModeling::Residue* residue, std::vector<std::string> &tags);
}

#endif
