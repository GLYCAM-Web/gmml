#include <iostream>  // for cout, can remove after debug
#include <algorithm> // sort

#include "includes/MolecularMetadata/GLYCAM/glycam06LinkageCodes.hpp"
#include "includes/CodeUtils/strings.hpp"

using gmml::MolecularMetadata::GLYCAM::Glycam06LinkageCodesLookupContainer;

//////////////////////////////////////////////////////////
//                      QUERY FUNCTIONS                 //
//////////////////////////////////////////////////////////

std::string Glycam06LinkageCodesLookupContainer::GetCodeForLinkages(std::string queryLinkage)
{ // e.g. input is 2,3, or 3,2 and output is "Z".
    for (const auto& entry : Glycam06LinkageCodesLookup_)
    {
        std::vector<std::string> linkages      = codeUtils::split(entry.linkage_, ','); // These are already sorted.
        std::vector<std::string> queryLinkages = codeUtils::split(queryLinkage, ',');
        std::sort(queryLinkages.begin(), queryLinkages.end()); // These might not be sorted
        if (std::equal(linkages.begin(), linkages.end(), queryLinkages.begin(), queryLinkages.end()))
        {
            return entry.glycamLinkageCode_;
        }
    }
    return "";
}

//////////////////////////////////////////////////////////
//                    INITIALIZER                       //
//////////////////////////////////////////////////////////

Glycam06LinkageCodesLookupContainer::Glycam06LinkageCodesLookupContainer()
{
    Glycam06LinkageCodesLookup_ = {
  // linkage_    , glycamLinkageCode_
        {"Terminal", "0"},
        {       "1", "1"},
        {       "2", "2"},
        {       "3", "3"},
        {       "4", "4"},
        {       "5", "5"},
        {       "6", "6"},
        {       "7", "7"},
        {       "8", "8"},
        {       "9", "9"},
        {     "2,3", "Z"},
        {     "2,4", "Y"},
        {     "2,6", "X"},
        {     "3,4", "W"},
        {     "3,6", "V"},
        {     "4,6", "U"},
        {   "2,3,4", "T"},
        {   "2,3,6", "S"},
        {   "2,4,6", "R"},
        {   "3,4,6", "Q"},
        { "2,3,4,6", "P"},
        {     "4,7", "K"},
        {     "4,8", "J"},
        {     "4,9", "I"},
        {     "7,8", "H"},
        {     "7,9", "G"},
        {     "8,9", "F"},
        {   "4,7,8", "E"},
        {   "4,7,9", "D"},
        {   "4,8,9", "C"},
        {   "7,8,9", "B"},
        { "4,7,8,9", "A"},
    };
}
