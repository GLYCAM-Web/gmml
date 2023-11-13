#include <algorithm> // sort
#include <iostream>
#include "includes/MolecularMetadata/GLYCAM/glycam06LinkageCodes.hpp"
#include "includes/CodeUtils/strings.hpp"

std::string GlycamMetadata::GetGlycam06ResidueLinkageCode(const std::string query)
{ // If something like 3,2,6 or 6,4,2 comes in, we need to sort it first.
    std::vector<std::string> queryLinkages = codeUtils::split(query, ',');
    std::sort(queryLinkages.begin(), queryLinkages.end());
    std::string sortedQuery = "";
    std::string deliminator = "";
    for (auto& element : queryLinkages)
    {
        sortedQuery += (deliminator + element);
        deliminator = ",";
    }
    static const std::unordered_map<std::string, std::string> Glycam06LinkageCodeLookup = {
  // linkage    , glycamLinkageCode
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
    return codeUtils::FindStringInStringMap(sortedQuery, Glycam06LinkageCodeLookup);
}
