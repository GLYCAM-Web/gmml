#include "includes/CodeUtils/strings.hpp"

std::string codeUtils::FindStringInStringMap(const std::string s, const std::unordered_map<std::string, std::string> &sMap)
{
    std::unordered_map<std::string, std::string>::const_iterator found = sMap.find(s);
    if (found != sMap.end())
    {
        return found->second;
    }
    return "";
}

std::string RemoveWhiteSpace(std::string s)
{
    s.erase(s.find_last_not_of(" ") + 1);
    s.erase(0, s.find_first_not_of(" "));
    return s;
}

//bool codeUtils::isStringInVector(const std::string s, const std::vector<std::string> &v)
//{
//  if (std::find(v.begin(), v.end(), s) != v.end())
//      return true;
//  return false;
//}
