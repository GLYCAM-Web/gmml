#ifndef GMML_INCLUDES_CODEUTILS_STRING_HPP
#define GMML_INCLUDES_CODEUTILS_STRING_HPP

#include <vector>
#include <sstream>

#include "includes/CodeUtils/logging.hpp"

namespace gmml
{

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
#endif //GMML_INCLUDES_CODEUTILS_FILES_HPP
