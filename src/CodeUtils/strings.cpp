#include "includes/CodeUtils/strings.hpp"
#include <sstream>
#include <iomanip> // std::setw
#include <ctype.h> // isdigit

bool codeUtils::startsWith(std::string bigString, std::string smallString)
{
	return (bigString.compare(0, smallString.length(), smallString) == 0);
}

std::string codeUtils::FindStringInStringMap(const std::string s, const std::unordered_map<std::string, std::string> &sMap)
{
    std::unordered_map<std::string, std::string>::const_iterator found = sMap.find(s);
    if (found != sMap.end())
    {
        return found->second;
    }
    return "";
}

std::string codeUtils::RemoveWhiteSpace(std::string str)
{
    codeUtils::RemoveSpaces(str);
//    s.erase(s.find_last_not_of(" ") + 1);
//    s.erase(0, s.find_first_not_of(" "));
    return str;
}

void codeUtils::RemoveSpaces(std::string& str)
{
    str.erase(std::remove_if(str.begin(), str.end(), ::isspace), str.end());
}

void codeUtils::RemoveQuotes(std::string& str)
{
    str.erase(std::remove(str.begin(), str.end(), '\"'), str.end());
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
