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

//bool codeUtils::isStringInVector(const std::string s, const std::vector<std::string> &v)
//{
//  if (std::find(v.begin(), v.end(), s) != v.end())
//      return true;
//  return false;
//}
