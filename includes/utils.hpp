#ifndef UTILS_HPP
#define UTILS_HPP

#include <string>
#include <stdexcept>
#include <sstream>
#include <algorithm>

namespace gmml
{
    // Remove spaces on both sides of the string.
    inline std::string& Trim(std::string& str)
    {
        str.erase(str.find_last_not_of(" ") + 1);
        str.erase(0, str.find_first_not_of(" "));
        return str;
    }

    inline void RemoveQuotes(std::string& str)
    {
        str.erase(std::remove(str.begin(), str.end(), '\"'), str.end());
    }

    inline void RemoveSpaces(std::string& str)
    {
        str.erase(std::remove_if(str.begin(), str.end(), ::isspace), str.end());
    }

    // String convertor
    template<typename T>
    inline T ConvertString(const std::string& str) {
        T val;
        std::stringstream ss(str);
        if (ss >> val)
            return val;

        throw std::invalid_argument("ConvertString: invalid conversion of string " + str);
    }
}

#endif // UTILS_HPP
