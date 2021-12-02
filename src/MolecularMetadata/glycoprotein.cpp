#include "includes/CodeUtils/logging.hpp"
#include "includes/CodeUtils/strings.hpp"
#include "includes/MolecularMetadata/glycoprotein.hpp"
#include "includes/MolecularModeling/Selections/selections.hpp"

std::string glycoproteinMetadata::LookupCodeForAminoAcidName(const std::string queryName)
{
    std::unordered_map<std::string, std::string> aminoAcidNameToCodeMap ({
                { "ALA", "A" },
                { "ARG", "R" },
                { "ASN", "N" },
                { "ASP", "D" },
                { "CYS", "C" },
                { "GLN", "Q" },
                { "GLU", "E" },
                { "GLY", "G" },
                { "HIS", "H" },
                { "ILE", "I" },
                { "LEU", "L" },
                { "LYS", "K" },
                { "MET", "M" },
                { "PHE", "F" },
                { "PRO", "P" },
                { "SER", "S" },
                { "THR", "T" },
                { "TRP", "W" },
                { "TYR", "Y" },
                { "VAL", "V" },
                { "NLN", "N" },
                { "OLT", "T" },
                { "OLS", "S" },
                { "OLY", "Y" },
                { "CLW", "W" }
        });
    return codeUtils::FindStringInStringMap(queryName, aminoAcidNameToCodeMap);
}

std::string glycoproteinMetadata::LookupLinkTypeForAminoAcidName(const std::string queryName)
{
    std::unordered_map<std::string, std::string> residueLinkMap ({
        { "ASN", "nLink" },
        { "THR", "oLink" },
        { "SER", "oLink" },
        { "TYR", "oLink" },
        { "TRP", "cLink" },
        { "NLN", "nLink" },
        { "OLT", "oLink" },
        { "OLS", "oLink" },
        { "OLY", "oLink" },
        { "CLW", "cLink" } });
    return codeUtils::FindStringInStringMap(queryName, residueLinkMap);
}

std::vector<std::string> glycoproteinMetadata::GetTagsForSequence(const std::string firstRes, const std::string midRes, const std::string lastRes)
{
    std::vector<std::string> foundTags;
    // N-links
    if ( (firstRes == "N") && (midRes != "P") )
    {
        if ( (lastRes == "T") || (lastRes == "S") )
        {
            foundTags.push_back("nLinkLikely");
        }
        if (lastRes == "C")
        {
            foundTags.push_back("nLinkCysSequon");
        }
    }
    // O-links // The midRes == "" is funky, but if I'm coming in here with more than just firstRes, I already have this tag from the initial call.
    if ( ( (firstRes == "S") || (firstRes == "T") ) && (midRes == "") )
    {
        foundTags.push_back("oLinkLikely");
    }
    if ( (firstRes == "Y") && (midRes == "") )
    {
        foundTags.push_back("oLinkTyr");
    }
    return foundTags;
}

std::string glycoproteinMetadata::GetSequenceContextAndDetermineTags(MolecularModeling::Residue* residue, std::vector<std::string> &tags)
{
    // Tags based on residue name only:
    std::string conjugationResidueCode = glycoproteinMetadata::LookupCodeForAminoAcidName(residue->GetName());
    std::vector<std::string> nameBasedTags = glycoproteinMetadata::GetTagsForSequence(conjugationResidueCode);
    tags.insert(tags.end(), nameBasedTags.begin(), nameBasedTags.end());
    // Now look for context around residue via atom connectivity.
    std::string precedingContext = "";
    MolecularModeling::Residue* firstPrecedingNeighbor = selection::FindNeighborResidueConnectedViaSpecificAtom(residue, "N");
    MolecularModeling::Residue* secondPrecedingNeighbor = nullptr;
    if (firstPrecedingNeighbor)
    {
        precedingContext = glycoproteinMetadata::LookupCodeForAminoAcidName(firstPrecedingNeighbor->GetName());
        secondPrecedingNeighbor = selection::FindNeighborResidueConnectedViaSpecificAtom(firstPrecedingNeighbor, "N");
        if (secondPrecedingNeighbor)
        {
            precedingContext = glycoproteinMetadata::LookupCodeForAminoAcidName(secondPrecedingNeighbor->GetName()) + precedingContext;
        }
    }
    std::string followingContext = "";
    MolecularModeling::Residue* firstFollowingNeighbor = selection::FindNeighborResidueConnectedViaSpecificAtom(residue, "C");
    MolecularModeling::Residue* secondFollowingNeighbor = nullptr;
    if (firstFollowingNeighbor)
    {
        std::string midResCode = glycoproteinMetadata::LookupCodeForAminoAcidName(firstFollowingNeighbor->GetName());
        followingContext = midResCode;
        secondFollowingNeighbor = selection::FindNeighborResidueConnectedViaSpecificAtom(firstFollowingNeighbor, "C");
        if (secondFollowingNeighbor)
        {
            std::string lastResCode = glycoproteinMetadata::LookupCodeForAminoAcidName(secondFollowingNeighbor->GetName());
            followingContext += lastResCode;
            std::vector<std::string> contextTags = glycoproteinMetadata::GetTagsForSequence(conjugationResidueCode, midResCode, lastResCode);
            tags.insert(tags.end(), contextTags.begin(), contextTags.end());
        }
    }
    std::string context = precedingContext + "_" + conjugationResidueCode + "_" + followingContext;
    gmml::log(__LINE__, __FILE__, gmml::INF, "Context: " + context);
    return context;
}
