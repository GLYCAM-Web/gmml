#ifndef UTILS_HPP
#define UTILS_HPP

#include <string>
#include <stdexcept>
#include <sstream>
#include <algorithm>
#include <vector>
#include "boost/tokenizer.hpp"
#include "boost/foreach.hpp"

namespace gmml
{
    /// Removes spaces on both sides of the string.
    inline std::string& Trim(std::string& str)
    {
        str.erase(str.find_last_not_of(" ") + 1);
        str.erase(0, str.find_first_not_of(" "));
        return str;
    }

    /// Removes quotation marks from the begining and the end of the given string
    inline void RemoveQuotes(std::string& str)
    {
        str.erase(std::remove(str.begin(), str.end(), '\"'), str.end());
    }

    /// Removes all spaces existing in the given string
    inline void RemoveSpaces(std::string& str)
    {
        str.erase(std::remove_if(str.begin(), str.end(), ::isspace), str.end());
    }

    /// String convertor
    template<typename T>
    inline T ConvertString(const std::string& str) {
        T val;
        std::stringstream ss(str);
        if (ss >> val)
            return val;

        throw std::invalid_argument("ConvertString: invalid conversion of string " + str);
    }

    inline std::vector<std::string> Split(std::string& line, std::string& delim)
    {
        boost::char_separator<char> separator(delim);
        boost::tokenizer< boost::char_separator<char> > tokens(line, separator);
        std::vector<string> vectorTokens = std::vector<std::string>();
        vectorTokens.assign(tokens.begin(), tokens.end());
        return vectorTokens;
    }
}

#endif // UTILS_HPP
