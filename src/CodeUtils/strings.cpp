#include "includes/CodeUtils/string.hpp"


bool gmml::startsWith(std::string bigString, std::string smallString)
{
	return (bigString.compare(0, smallString.length(), smallString) == 0);
}

//std::vector<std::string> split(const std::string &s, char delim)
//{
//    std::vector<std::string> elems;
//    split(s, delim, std::back_inserter(elems));
//    return elems;
//}
