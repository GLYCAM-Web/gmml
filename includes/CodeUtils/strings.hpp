#ifndef GMML_INCLUDES_CODEUTILS_STRINGS_HPP
#define GMML_INCLUDES_CODEUTILS_STRINGS_HPP
#include <string>
//#include <vector>
#include <unordered_map>

namespace codeUtils
{
     std::string FindStringInStringMap(const std::string s, const std::unordered_map<std::string, std::string> &sMap);
     //bool isStringInVector(const std::string s, const std::vector<std::string> &v)
}
#endif //GMML_INCLUDES_CODEUTILS_STRINGS_HPP
