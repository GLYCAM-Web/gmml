#ifndef UTILS_HPP
#define UTILS_HPP

#include <string>
#include <stdexcept>
#include <sstream>
#include <algorithm>
#include <vector>
#include <iomanip>
#include "boost/tokenizer.hpp"
#include "boost/foreach.hpp"
#include "common.hpp"

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

    inline std::vector<std::string> Split(std::string line, std::string delim)
    {
        boost::char_separator<char> separator(delim.c_str());
        boost::tokenizer< boost::char_separator<char> > tokens(line, separator);
        std::vector<std::string> vectorTokens = std::vector<std::string>();
        vectorTokens.assign(tokens.begin(), tokens.end());
        return vectorTokens;
    }

    /// String convertor
    template<typename T>
    inline T ConvertString(const std::string& str) {
        T val;
        std::stringstream ss(str);
        ss >> val;
//        if (ss >> val)
        return val;

        throw std::invalid_argument("ConvertString: invalid conversion of string " + str);
    }

    /// Expand a given line to a desired length by adding space at the end of the original one
    inline std::string& ExpandLine(std::string& line, int length)
    {
        if(line.length() >= length)
            return line;
        else
        {
            int space = length - line.length();
            std::stringstream ss;
            ss << line << std::setw(space) << " ";
            line = ss.str();
            return line;
        }
    }

    inline InputFileType ConvertString2AssemblyInputFileType(std::string type)
    {
        if(type.compare("PDB") == 0)
            return PDB;
        if(type.compare("LIB") == 0)
            return LIB;
        if(type.compare("PREP") == 0)
            return PREP;
        if(type.compare("TOP") == 0)
            return TOP;
        if(type.compare("TOP_CRD") == 0)
            return TOP_CRD;
        if(type.compare("MULTIPLE") == 0)
            return MULTIPLE;
    }
    inline std::string ConvertAssemblyInputFileType2String(InputFileType type)
    {
        switch(type)
        {
            case PDB:
                return "PDB";
            case LIB:
                return "LIB";
            case PREP:
                return "PREP";
            case TOP:
                return "TOP";
            case TOP_CRD:
                return "TOP_CRD";
            case MULTIPLE:
                return "MULTIPLE";
        }
    }
}

#endif // UTILS_HPP
