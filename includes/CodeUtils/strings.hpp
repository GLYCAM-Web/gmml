#ifndef GMML_INCLUDES_CODEUTILS_STRINGS_HPP
#define GMML_INCLUDES_CODEUTILS_STRINGS_HPP
#include <string>
//#include <vector>
#include <unordered_map>

namespace codeUtils
{
// Functions
std::string FindStringInStringMap(const std::string s, const std::unordered_map<std::string, std::string> &sMap);
std::string RemoveWhiteSpace(std::string s);
void ExpandLine(std::string &line, int length);
int GetSizeOfIntInString(const std::string str);
}
#endif //GMML_INCLUDES_CODEUTILS_STRINGS_HPP
