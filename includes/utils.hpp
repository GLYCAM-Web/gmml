#ifndef UTILS_HPP
#define UTILS_HPP

#include <string>
#include <stdexcept>
#include <sstream>

namespace gmml
{
    // Remove spaces on both sides of the string.
    inline std::string& trim(std::string& str)
    {
        str.erase(str.find_last_not_of(" ") + 1);
        str.erase(0, str.find_first_not_of(" "));
        return str;
    }

    // String convertor
    template<typename T>
    inline T convert_string(const std::string& str) {
        T val;
        std::stringstream ss(str);
        if (ss >> val)
            return val;

        throw std::invalid_argument("convert_string: invalid conversion of string " + str);
    }
}

#endif // UTILS_HPP
