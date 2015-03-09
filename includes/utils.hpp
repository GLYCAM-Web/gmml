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
    /*! \fn
      * Removes spaces on both sides of the string.
      * @param str String with spaces either in the beginning or at the end
      * @return Given string without spaces appeared in the beginning or at the end in the original one
      */
    inline std::string& Trim(std::string& str)
    {
        str.erase(str.find_last_not_of(" ") + 1);
        str.erase(0, str.find_first_not_of(" "));
        return str;
    }

    /*! \fn
      * Removes quotation marks from the begining and the end of the given string
      * @param @param str String with qouates
      * @return Given string without qouates appeared in the original one
      */
    inline void RemoveQuotes(std::string& str)
    {
        str.erase(std::remove(str.begin(), str.end(), '\"'), str.end());
    }

    /*! \fn
      * Removes all spaces existing in the given string
      * @param str String with space characters in it (anywhere in the string)
      * @return String without any space in between
      */
    inline void RemoveSpaces(std::string& str)
    {
        str.erase(std::remove_if(str.begin(), str.end(), ::isspace), str.end());
    }

    /*! \fn
      * Split a line (series of characters) with the given delimiter
      * @param line Input string in order to be split by the given delimiter
      * @param delim Series of delimiters (each delimiter character is followed by the next one in a single string variable) in order to split the given string
      * @return Vector of elements of the given string that have been split by the given delimiter(s)
      */
    inline std::vector<std::string> Split(std::string line, std::string delim)
    {
        boost::char_separator<char> separator(delim.c_str());
        boost::tokenizer< boost::char_separator<char> > tokens(line, separator);
        std::vector<std::string> vectorTokens = std::vector<std::string>();
        vectorTokens.assign(tokens.begin(), tokens.end());
        return vectorTokens;
    }

    /*! \fn
      * Convert string to the defined type
      * @param str String that has to be converted to the defined type
      * @return The given type version of the given string
      */
    template<typename T>
    inline T ConvertString(const std::string& str) {
        T val;
        std::stringstream ss(str);
        ss >> val;
//        if (ss >> val)
        return val;

        throw std::invalid_argument("ConvertString: invalid conversion of string " + str);
    }

    /*! \fn
      * Convert T (the given type) to string
      * @param T the given type has to be converted to string
      * @return The string version of the given type
      */
    template<typename T>
    std::string ConvertT(const T& given_type) {
        std::stringstream ss;
        ss << given_type;
        return ss.str();
    }

    /*! \fn
      * Expand a given line to a desired length by adding space at the end of the original one
      * @param line A line that have to be in a defined length
      * @param length Thel fixed length that the line has to be
      * @return An expanded line into the given length
      */
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

    /*! \fn
      * Convert string version of input file type to the corresponding enum value
      * @param type String indicates input file type
      * @return A value selected from InputFileType enumerator correspondence to the given string
      */
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

    /*! \fn
      * Convert a value of InputFileType enumerator to the string version of it
      * @param type A value of InputFileType has to be converted to string
      * @return String format of the given InputFileType enumerator
      */
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
            default:
                return "UNKNOWN";
        }
    }

    /*! \fn
      * Convert string version of assembly building structure option to the corresponding enum value
      * @param type String indicates assembly building structure option
      * @return A value selected from AssemblyBuildingStructureOption enumerator correspondence to the given string
      */
    inline BuildingStructureOption ConvertString2AssemblyBuildingStructureOption(std::string option)
    {
        if(option.compare("Distance") == 0)
            return DISTANCE;
        if(option.compare("Original") == 0)
            return ORIGINAL;
        if(option.compare("Database") == 0)
            return DATABASE;
        else
            return DISTANCE;
    }

    /*! \fn
      * Convert a value of AssemblyBuildingStrucutreOption enumerator to the string version of it
      * @param type A value of AssemblyBuildingStrucutreOption has to be converted to string
      * @return String format of the given value of AssemblyBuildingStrucutreOption enumerator
      */
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
            default:
                return "Distance";
        }
    }

    /*! \fn
      * Convert degree to radian
      * @param degree Magnitude of an angle in degree
      * @return Magnitude of the given angle in radian
      */
    inline double ConvertDegree2Radian(double degree)
    {
        return degree/PI_DEGREE*gmml::PI_RADIAN;
    }

    /*! \fn
      * Convert radian to degree
      * @param radian Magnitude of an angle in radian
      * @return Magnitude of the given angle in degree
      */
    inline double ConvertRadian2Degree(double radian)
    {
        return radian*PI_DEGREE/PI_RADIAN;
    }

    /*! \fn
      * Convert internal coordinate to the corresponding cartesian coordinate
      * @param coordinate_list List of at most three internal coordinates in order to calculate the cartesian coordinate of the given internal coordinate (distance, angle, torsion)
      * @param distance X value of the internal coordinate
      * @param angle Y value of the interanl coordinate
      * @param torsion Z value of the internal coordinate
      * @return Cartesian coordinate of the internal coordinate (distance, angle, torsion)
      */
    inline Geometry::Coordinate* ConvertInternalCoordinate2CartesianCoordinate(std::vector<Geometry::Coordinate*> coordinate_list, double distance, double angle, double torsion)
    {
        if(coordinate_list.size() == 0)
        {
            Geometry::Coordinate* coordinate = new Geometry::Coordinate();
            coordinate->Print(std::cout);
            std::cout << std::endl;
            return coordinate;
        }
        if(coordinate_list.size() == 1)
        {
            coordinate_list.at(0)->Print(std::cout);
            std::cout << std::endl;
            Geometry::Coordinate* coordinate = new Geometry::Coordinate(coordinate_list.at(0)->GetX() + distance, 0.0, 0.0);
            coordinate->Print(std::cout);
            std::cout << std::endl;
            return coordinate;
        }
        if(coordinate_list.size() == 2)
        {
            coordinate_list.at(0)->Print(std::cout);
            std::cout << std::endl;
            coordinate_list.at(1)->Print(std::cout);
            std::cout << std::endl;
            Geometry::Coordinate* coordinate = new Geometry::Coordinate(coordinate_list.at(1)->GetX() - cos(gmml::ConvertDegree2Radian(angle) * distance),
                                                                        sin(gmml::ConvertDegree2Radian(angle)) * distance, 0.0);
            coordinate->Print(std::cout);
            std::cout << std::endl;
            return coordinate;
        }
        else
        {
            coordinate_list.at(0)->Print(std::cout);
            std::cout << std::endl;
            coordinate_list.at(1)->Print(std::cout);
            std::cout << std::endl;
            coordinate_list.at(2)->Print(std::cout);
            std::cout << std::endl;
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

            Geometry::Coordinate* coordinate = new Geometry::Coordinate(p.GetX() * v.at(0) + r.GetX() * v.at(1) + v2.GetX() * v.at(2) + parent_vector.GetX() * v.at(3),
                                                                        p.GetY() * v.at(0) + r.GetY() * v.at(1) + v2.GetY() * v.at(2) + parent_vector.GetY() * v.at(3),
                                                                        p.GetZ() * v.at(0) + r.GetZ() * v.at(1) + v2.GetZ() * v.at(2) + parent_vector.GetZ() * v.at(3));
            coordinate->Print(std::cout);
            std::cout << std::endl;
            return coordinate;
        }
    }

    inline Geometry::Coordinate* ConvertCartesianCoordinate2InternalCoordinate(Geometry::Coordinate* coordinate, std::vector<Geometry::Coordinate*> coordinate_list)
    {
        if(coordinate_list.size() == 0)
            return new Geometry::Coordinate();
        if(coordinate_list.size() == 1)
        {
            coordinate_list.at(0)->Print(std::cout);
            std::cout << std::endl;
            coordinate->Print(std::cout);
            std::cout << std::endl;
            Geometry::Coordinate parent_vector = Geometry::Coordinate(*coordinate_list.at(0));
            double distance = coordinate->Distance(parent_vector);
            return new Geometry::Coordinate(distance, 0.0, 0.0);
        }
        if(coordinate_list.size() == 2)
        {
            Geometry::Coordinate grandparent_vector = Geometry::Coordinate(*coordinate_list.at(0));
            Geometry::Coordinate parent_vector = Geometry::Coordinate(*coordinate_list.at(1));
            double distance = coordinate->Distance(parent_vector);

            Geometry::Coordinate dist_current_parent_vector = Geometry::Coordinate(*coordinate);
            Geometry::Coordinate dist_grandparent_parent_vector = Geometry::Coordinate(grandparent_vector);
            dist_current_parent_vector.operator -(parent_vector);
            dist_grandparent_parent_vector.operator -(parent_vector);
            double dist_current_parent = dist_current_parent_vector.length();
            double dist_grandparent_parent = dist_grandparent_parent_vector.length();
            double dist_current_parent_dot_grandparent_parent = dist_current_parent_vector.DotProduct(dist_grandparent_parent_vector);
            double angle = ConvertRadian2Degree(acos(dist_current_parent_dot_grandparent_parent/(dist_current_parent * dist_grandparent_parent)));

            return new Geometry::Coordinate(distance, angle, 0.0);
        }
        else
        {
            Geometry::Coordinate greatgrandparent_vector = Geometry::Coordinate(*coordinate_list.at(0));
            Geometry::Coordinate grandparent_vector = Geometry::Coordinate(*coordinate_list.at(1));
            Geometry::Coordinate parent_vector = Geometry::Coordinate(*coordinate_list.at(2));
            double distance = coordinate->Distance(parent_vector);

            Geometry::Coordinate dist_current_parent_vector = Geometry::Coordinate(*coordinate);
            Geometry::Coordinate dist_grandparent_parent_vector = Geometry::Coordinate(grandparent_vector);
            dist_current_parent_vector.operator -(parent_vector);
            dist_grandparent_parent_vector.operator -(parent_vector);
            double dist_current_parent = dist_current_parent_vector.length();
            double dist_grandparent_parent = dist_grandparent_parent_vector.length();
            double dist_current_parent_dot_grandparent_parent = dist_current_parent_vector.DotProduct(dist_grandparent_parent_vector);
            double angle = ConvertRadian2Degree(acos(dist_current_parent_dot_grandparent_parent/(dist_current_parent * dist_grandparent_parent)));

            Geometry::Coordinate dist_parent_current_vector = Geometry::Coordinate(parent_vector);
            dist_parent_current_vector.operator -(*coordinate);
            Geometry::Coordinate dist_grandparent_parent_vector_1 = Geometry::Coordinate(grandparent_vector);
            dist_grandparent_parent_vector_1.operator -(parent_vector);
            Geometry::Coordinate dist_greatgrandparent_grandparent_vector = Geometry::Coordinate(greatgrandparent_vector);
            dist_greatgrandparent_grandparent_vector.operator -(grandparent_vector);
            Geometry::Coordinate dist_grandparent_parent_cross_dist_greatgrandparent_grandparent_vector =
                    Geometry::Coordinate(dist_grandparent_parent_vector);
            dist_grandparent_parent_cross_dist_greatgrandparent_grandparent_vector.CrossProduct(dist_greatgrandparent_grandparent_vector);
            Geometry::Coordinate dist_parent_current_cross_dist_grandparent_parent_vector = Geometry::Coordinate(dist_parent_current_vector);
            dist_parent_current_cross_dist_grandparent_parent_vector.CrossProduct(dist_grandparent_parent_vector_1);
            Geometry::Coordinate dist_parent_current_multiply_dist_grandparent_parent_vector = Geometry::Coordinate(dist_parent_current_vector);
            dist_parent_current_multiply_dist_grandparent_parent_vector.operator *(dist_grandparent_parent_vector.length());

            double torsion = ConvertRadian2Degree(
                        atan2(dist_parent_current_multiply_dist_grandparent_parent_vector.DotProduct(
                                  dist_grandparent_parent_cross_dist_greatgrandparent_grandparent_vector),
                              dist_parent_current_cross_dist_grandparent_parent_vector.DotProduct(
                                  dist_grandparent_parent_cross_dist_greatgrandparent_grandparent_vector)));

            return new Geometry::Coordinate(distance, angle, torsion);
        }
    }
}

#endif // UTILS_HPP
