#ifndef GMML_INCLUDES_MOLECULARMETADATA_GLYCOPROTEIN_HPP
#define GMML_INCLUDES_MOLECULARMETADATA_GLYCOPROTEIN_HPP

#include <string>
#include <vector>
#include "includes/CentralDataStructure/residue.hpp"

namespace glycoproteinMetadata
{
    std::string LookupCodeForAminoAcidName(const std::string queryName);
    std::string LookupLinkTypeForAminoAcidName(const std::string queryName);
    std::vector<std::string> GetTagsForSequence(const std::string firstRes = "", const std::string midRes = "",
                                                const std::string lastRes = "");
    std::string GetSequenceContextAndDetermineTags(cds::Residue* residue, std::vector<std::string>& tags);
    std::string ConvertGlycosylatedResidueName(const std::string queryname);
    std::string GetGlycositeConnectionAtomName(const std::string queryname);
} // namespace glycoproteinMetadata

#endif
