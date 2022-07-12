#ifndef GMML_INCLUDES_CODEUTILS_STRINGS_HPP
#define GMML_INCLUDES_CODEUTILS_STRINGS_HPP

#include <string>
#include <unordered_map>
#include <vector>
#include <sstream>

namespace codeUtils
{
// Functions
bool startsWith(std::string bigString, std::string smallString);
std::string FindStringInStringMap(const std::string s, const std::unordered_map<std::string, std::string> &sMap);
std::string RemoveWhiteSpace(std::string s);
void ExpandLine(std::string &line, int length);
int GetSizeOfIntInString(const std::string str);

template<typename Out>
inline void split(const std::string &s, char delim, Out result)
{
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim))
    {
        *(result++) = item;
    }
}

inline std::vector<std::string> split(const std::string &s, char delim)
{
    std::vector<std::string> elems;
    split(s, delim, std::back_inserter(elems));
    return elems;
}

}
#endif //GMML_INCLUDES_CODEUTILS_STRINGS_HPP
