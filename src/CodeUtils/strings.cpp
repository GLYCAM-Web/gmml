#include <sstream>
#include <iomanip> // std::setw
#include <ctype.h> // isdigit
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

std::string codeUtils::RemoveWhiteSpace(std::string s)
{
    s.erase(s.find_last_not_of(" ") + 1);
    s.erase(0, s.find_first_not_of(" "));
    return s;
}

void codeUtils::ExpandLine(std::string &line, int length)
{
    int l = line.length();
    if(l >= length)
        return;
    else
    {
        int space = length - l;
        std::stringstream ss;
        ss << line << std::setw(space) << " ";
        line = ss.str();
        return;
    }
    return;
}

int codeUtils::GetSizeOfIntInString(const std::string str)
{
    int size = 0;
    for(const char& c : str)
    {
        if (isdigit(c))
        {
            ++size;
        }
        else
        {
            return size;
        }
    }
    return size;
}


//bool codeUtils::isStringInVector(const std::string s, const std::vector<std::string> &v)
//{
//  if (std::find(v.begin(), v.end(), s) != v.end())
//      return true;
//  return false;
//}

bool gmml::startsWith(std::string bigString, std::string smallString)
{
	return (bigString.compare(0, smallString.length(), smallString) == 0);
}