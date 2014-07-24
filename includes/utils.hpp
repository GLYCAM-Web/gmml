#ifndef UTILS_HPP
#define UTILS_HPP

#include <string>
#include <stdexcept>
#include <sstream>
#include <algorithm>
#include <vector>
#include <iomanip>
#include <math.h>
#include "boost/tokenizer.hpp"
#include "boost/foreach.hpp"
#include "common.hpp"
#include "Geometry/coordinate.hpp"

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
        if((int)line.length() >= length)
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
        return UNKNOWN;
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
            case UNKNOWN:
                return "UNKNOWN";
        }
    }

    inline BuildingStructureOption ConvertString2AssemblyBuildingStructureOption(std::string option)
    {
        if(option.compare("Distance") == 0)
            return DISTANCE;
        if(option.compare("Original") == 0)
            return ORIGINAL;
        if(option.compare("Database") == 0)
            return DATABASE;
    }

    inline std::string ConvertAssemblyBuildingStructureOption2String(BuildingStructureOption option)
    {
        switch(option)
        {
            case DISTANCE:
                return "Distance";
            case ORIGINAL:
                return "Original";
            case DATABASE:
                return "Database";
        }
    }

    inline double ConvertDegree2Radian(double degree)
    {
        return degree/PI_DEGREE*gmml::PI_RADIAN;
    }

    inline Geometry::Coordinate* ConvertInternalCoordinate2CartesianCoordinate(std::vector<Geometry::Coordinate*> coordinate_list, double distance, double angle, double torsion)
    {
        if(coordinate_list.size() == 0)
            return new Geometry::Coordinate();
        if(coordinate_list.size() == 1)
            return new Geometry::Coordinate(coordinate_list.at(0)->GetX() + distance, 0.0, 0.0);
        if(coordinate_list.size() == 2)
            return new Geometry::Coordinate(coordinate_list.at(1)->GetX() - cos(gmml::ConvertDegree2Radian(angle) * distance),
                                            sin(gmml::ConvertDegree2Radian(angle)) * distance, 0.0);
        else
        {
            torsion = gmml::PI_DEGREE - torsion;

            Geometry::Coordinate great_grandparent_vector = Geometry::Coordinate(coordinate_list.at(0)->GetX(), coordinate_list.at(0)->GetY(), coordinate_list.at(0)->GetZ());
            Geometry::Coordinate grandparent_vector = Geometry::Coordinate(coordinate_list.at(1)->GetX(), coordinate_list.at(1)->GetY(), coordinate_list.at(1)->GetZ());
            Geometry::Coordinate parent_vector = Geometry::Coordinate(coordinate_list.at(2)->GetX(), coordinate_list.at(2)->GetY(), coordinate_list.at(2)->GetZ());

            Geometry::Coordinate v1 = Geometry::Coordinate(grandparent_vector);
            Geometry::Coordinate v2 = Geometry::Coordinate(parent_vector);

            v1.operator-(great_grandparent_vector);
            v2.operator-(grandparent_vector);

            v1.Normalize();
            v2.Normalize();

            if(abs(v1.GetX() + v2.GetX()) < gmml::EPSILON &&
                    abs(v1.GetY() + v2.GetY()) < gmml::EPSILON &&
                    abs(v1.GetZ() + v2.GetZ()) < gmml::EPSILON)
            {
                great_grandparent_vector.Translate(1.0, -1.0, 2);
                v1 = Geometry::Coordinate(grandparent_vector);
                v1.operator-(great_grandparent_vector);
            }
            Geometry::Coordinate r = Geometry::Coordinate(v1);
            r.CrossProduct(v2);

            r.Normalize();

            Geometry::Coordinate p = Geometry::Coordinate(r);
            p.CrossProduct(v2);

            std::vector<double> v = std::vector<double>();

            v.push_back(distance * sin(gmml::ConvertDegree2Radian(angle)) * cos(gmml::ConvertDegree2Radian(torsion)));
            v.push_back(distance * sin(gmml::ConvertDegree2Radian(angle)) * sin(gmml::ConvertDegree2Radian(torsion)));
            v.push_back(distance * cos(gmml::ConvertDegree2Radian(angle)));
            v.push_back(1.0);

            return new Geometry::Coordinate(p.GetX() * v.at(0) + r.GetX() * v.at(1) + v2.GetX() * v.at(2) + parent_vector.GetX() * v.at(3),
                                            p.GetY() * v.at(0) + r.GetY() * v.at(1) + v2.GetY() * v.at(2) + parent_vector.GetY() * v.at(3),
                                            p.GetZ() * v.at(0) + r.GetZ() * v.at(1) + v2.GetZ() * v.at(2) + parent_vector.GetZ() * v.at(3));

        }
    }

}

#endif // UTILS_HPP
