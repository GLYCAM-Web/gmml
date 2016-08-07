#ifndef UTILS_HPP
#define UTILS_HPP

#include <ctime>
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
#include "GeometryTopology/coordinate.hpp"

#include <fstream>

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
    inline std::string ExpandLine(std::string line, int length)
    {
        int l = line.length();
        if(l >= length)
            return line;
        else
        {
            int space = length - l;
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
        if(type.compare("PDBQT") == 0)
            return PDBQT;
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
            case PDBQT:
                return "PDBQT";
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

    inline std::string ConvertTopologicalType2String(TopologicalType type)
    {
        switch(type)
        {
            case kTopTypeE:
                return "E";
            case kTopTypeS:
                return "S";
            case kTopTypeB:
                return "B";
            case kTopType3:
                return "3";
            case kTopType4:
                return "4";
            case kTopTypeM:
                return "M";
            default:
                return "E";
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
    inline GeometryTopology::Coordinate* ConvertInternalCoordinate2CartesianCoordinate(std::vector<GeometryTopology::Coordinate*> coordinate_list, double distance, double angle, double torsion)
    {
        if(coordinate_list.size() == 0)
        {
            GeometryTopology::Coordinate* coordinate = new GeometryTopology::Coordinate();
            return coordinate;
        }
        if(coordinate_list.size() == 1)
        {
            GeometryTopology::Coordinate* coordinate = new GeometryTopology::Coordinate(coordinate_list.at(0)->GetX() + distance, 0.0, 0.0);
            return coordinate;
        }
        if(coordinate_list.size() == 2)
        {
            GeometryTopology::Coordinate* coordinate = new GeometryTopology::Coordinate(coordinate_list.at(1)->GetX() - cos(gmml::ConvertDegree2Radian(angle) * distance),
                                                                        sin(gmml::ConvertDegree2Radian(angle)) * distance, 0.0);
            return coordinate;
        }
        else
        {
            torsion = gmml::PI_DEGREE - torsion;

            GeometryTopology::Coordinate great_grandparent_vector = GeometryTopology::Coordinate(coordinate_list.at(0)->GetX(), coordinate_list.at(0)->GetY(), coordinate_list.at(0)->GetZ());
            GeometryTopology::Coordinate grandparent_vector = GeometryTopology::Coordinate(coordinate_list.at(1)->GetX(), coordinate_list.at(1)->GetY(), coordinate_list.at(1)->GetZ());
            GeometryTopology::Coordinate parent_vector = GeometryTopology::Coordinate(coordinate_list.at(2)->GetX(), coordinate_list.at(2)->GetY(), coordinate_list.at(2)->GetZ());

            GeometryTopology::Coordinate v1 = GeometryTopology::Coordinate(great_grandparent_vector);
            GeometryTopology::Coordinate v2 = GeometryTopology::Coordinate(grandparent_vector);

            v1.operator-(grandparent_vector);
            v2.operator-(parent_vector);

            v1.Normalize();
            v2.Normalize();

            if(abs(v1.GetX() + v2.GetX()) < gmml::EPSILON &&
                    abs(v1.GetY() + v2.GetY()) < gmml::EPSILON &&
                    abs(v1.GetZ() + v2.GetZ()) < gmml::EPSILON)
            {
                great_grandparent_vector.Translate(10.0, -1.0, 3);
                v1 = GeometryTopology::Coordinate(grandparent_vector);
                v1.operator-(great_grandparent_vector);
            }

            GeometryTopology::Coordinate v1_cross_v2 = GeometryTopology::Coordinate(v1);
            v1_cross_v2.CrossProduct(v2);
            v1_cross_v2.Normalize();
            double matrix[3][4];
            matrix[0][1] = v1_cross_v2.GetX();
            matrix[1][1] = v1_cross_v2.GetY();
            matrix[2][1] = v1_cross_v2.GetZ();

            matrix[0][2] = v2.GetX();
            matrix[1][2] = v2.GetY();
            matrix[2][2] = v2.GetZ();

            GeometryTopology::Coordinate v1_cross_v2_cross_v2 = GeometryTopology::Coordinate(v1_cross_v2);
            v1_cross_v2_cross_v2.CrossProduct(v2);

            matrix[0][0] = v1_cross_v2_cross_v2.GetX();
            matrix[1][0] = v1_cross_v2_cross_v2.GetY();
            matrix[2][0] = v1_cross_v2_cross_v2.GetZ();

            matrix[0][3] = parent_vector.GetX();
            matrix[1][3] = parent_vector.GetY();
            matrix[2][3] = parent_vector.GetZ();

            GeometryTopology::Coordinate v = GeometryTopology::Coordinate(distance * sin(ConvertDegree2Radian(angle)) * cos(ConvertDegree2Radian(torsion)),
                                                                          distance * sin(ConvertDegree2Radian(angle)) * sin(ConvertDegree2Radian(torsion)),
                                                                          distance * cos(ConvertDegree2Radian(angle)));

            GeometryTopology::Coordinate* coordinate = new GeometryTopology::Coordinate(matrix[0][0] * v.GetX() + matrix[0][1] * v.GetY() + matrix[0][2] * v.GetZ() + matrix[0][3],
                                                                                        matrix[1][0] * v.GetX() + matrix[1][1] * v.GetY() + matrix[1][2] * v.GetZ() + matrix[1][3],
                                                                                        matrix[2][0] * v.GetX() + matrix[2][1] * v.GetY() + matrix[2][2] * v.GetZ() + matrix[2][3]);
            return coordinate;
        }
    }

    inline GeometryTopology::Coordinate* ConvertCartesianCoordinate2InternalCoordinate(GeometryTopology::Coordinate* coordinate, std::vector<GeometryTopology::Coordinate*> coordinate_list)
    {
        if(coordinate_list.size() == 0)
            return new GeometryTopology::Coordinate();
        if(coordinate_list.size() == 1)
        {
            GeometryTopology::Coordinate parent_vector = GeometryTopology::Coordinate(*coordinate_list.at(0));
            double distance = coordinate->Distance(parent_vector);
            return new GeometryTopology::Coordinate(distance, 0.0, 0.0);
        }
        if(coordinate_list.size() == 2)
        {
            GeometryTopology::Coordinate grandparent_vector = GeometryTopology::Coordinate(*coordinate_list.at(0));
            GeometryTopology::Coordinate parent_vector = GeometryTopology::Coordinate(*coordinate_list.at(1));
            double distance = coordinate->Distance(parent_vector);

            GeometryTopology::Coordinate dist_current_parent_vector = GeometryTopology::Coordinate(*coordinate);
            GeometryTopology::Coordinate dist_grandparent_parent_vector = GeometryTopology::Coordinate(grandparent_vector);
            dist_current_parent_vector.operator -(parent_vector);
            dist_grandparent_parent_vector.operator -(parent_vector);
            double dist_current_parent = dist_current_parent_vector.length();
            double dist_grandparent_parent = dist_grandparent_parent_vector.length();
            double dist_current_parent_dot_grandparent_parent = dist_current_parent_vector.DotProduct(dist_grandparent_parent_vector);
            double angle = ConvertRadian2Degree(acos(dist_current_parent_dot_grandparent_parent/(dist_current_parent * dist_grandparent_parent)));

            return new GeometryTopology::Coordinate(distance, angle, 0.0);
        }
        else
        {
            GeometryTopology::Coordinate greatgrandparent_vector = GeometryTopology::Coordinate(*coordinate_list.at(0));
            GeometryTopology::Coordinate grandparent_vector = GeometryTopology::Coordinate(*coordinate_list.at(1));
            GeometryTopology::Coordinate parent_vector = GeometryTopology::Coordinate(*coordinate_list.at(2));
            double distance = coordinate->Distance(parent_vector);

            GeometryTopology::Coordinate dist_current_parent_vector = GeometryTopology::Coordinate(*coordinate);
            GeometryTopology::Coordinate dist_grandparent_parent_vector = GeometryTopology::Coordinate(grandparent_vector);
            dist_current_parent_vector.operator -(parent_vector);
            dist_grandparent_parent_vector.operator -(parent_vector);
            double dist_current_parent = dist_current_parent_vector.length();
            double dist_grandparent_parent = dist_grandparent_parent_vector.length();
            double dist_current_parent_dot_grandparent_parent = dist_current_parent_vector.DotProduct(dist_grandparent_parent_vector);
            double angle = ConvertRadian2Degree(acos(dist_current_parent_dot_grandparent_parent/(dist_current_parent * dist_grandparent_parent)));

            GeometryTopology::Coordinate dist_parent_current_vector = GeometryTopology::Coordinate(parent_vector);
            dist_parent_current_vector.operator -(*coordinate);
            GeometryTopology::Coordinate dist_grandparent_parent_vector_1 = GeometryTopology::Coordinate(grandparent_vector);
            dist_grandparent_parent_vector_1.operator -(parent_vector);
            GeometryTopology::Coordinate dist_greatgrandparent_grandparent_vector = GeometryTopology::Coordinate(greatgrandparent_vector);
            dist_greatgrandparent_grandparent_vector.operator -(grandparent_vector);
            GeometryTopology::Coordinate dist_grandparent_parent_cross_dist_greatgrandparent_grandparent_vector =
                    GeometryTopology::Coordinate(dist_grandparent_parent_vector);
            dist_grandparent_parent_cross_dist_greatgrandparent_grandparent_vector.CrossProduct(dist_greatgrandparent_grandparent_vector);
            GeometryTopology::Coordinate dist_parent_current_cross_dist_grandparent_parent_vector = GeometryTopology::Coordinate(dist_parent_current_vector);
            dist_parent_current_cross_dist_grandparent_parent_vector.CrossProduct(dist_grandparent_parent_vector_1);
            GeometryTopology::Coordinate dist_parent_current_multiply_dist_grandparent_parent_vector = GeometryTopology::Coordinate(dist_parent_current_vector);
            dist_parent_current_multiply_dist_grandparent_parent_vector.operator *(dist_grandparent_parent_vector.length());

            double torsion = ConvertRadian2Degree(
                        atan2(dist_parent_current_multiply_dist_grandparent_parent_vector.DotProduct(
                                  dist_grandparent_parent_cross_dist_greatgrandparent_grandparent_vector),
                              dist_parent_current_cross_dist_grandparent_parent_vector.DotProduct(
                                  dist_grandparent_parent_cross_dist_greatgrandparent_grandparent_vector)));

            return new GeometryTopology::Coordinate(distance, angle, torsion);
        }
    }

    /*! \fn
      * A function in order to look up the stereochemistry name of the sugar structure based on the given string version of the chemical code structure
      * @param code The string chemical code structure
      * @return SUGARNAMELOOKUP The matched row of the lookup table with the given code
      */
    inline Glycan::SugarName SugarStereoChemistryNameLookup(std::string code)
    {
        for(int i = 0; i < SUGARNAMELOOKUPSIZE; i++)
        {
            if(code.compare(SUGARNAMELOOKUP[i].chemical_code_string_) == 0)
                return SUGARNAMELOOKUP[i];
        }
        return SUGARNAMELOOKUP[0];
    }

    /*! \fn
      * A function in order to search the stereochemistry lookup table and identify the closest match for the sugar structure based on the given string version of the chemical code structure
      * @param code The string chemical code structure
      * @return SUGARNAMELOOKUP The closest row of the lookup table which matches th emost with the given code
      */
    inline Glycan::SugarName ClosestMatchSugarStereoChemistryNameLookup(std::string code, std::vector<Glycan::SugarName>& closest_matches)
    {
        std::string vocab[] = {"2", "3", "4", "a", "+1", "+2", "+3", "-1"};
        int vocab_size = (sizeof(vocab)/sizeof(vocab[0]));
        std::string pos = "^_";
        std::string stat = "d";
        int min_diff_score = -1000;
        std::map<int, std::vector<Glycan::SugarName> > score_map = std::map<int, std::vector<Glycan::SugarName> >();
        for(int i = 0; i < SUGARNAMELOOKUPSIZE; i++)
        {
            int diff_score = 0;
            std::string target_code = SUGARNAMELOOKUP[i].chemical_code_string_;
            if((target_code.find("P") != std::string::npos && code.find("P") != std::string::npos) ||
                    (target_code.find("F") != std::string::npos && code.find("F") != std::string::npos))
            {
                for(int j = 0; j < vocab_size; j++)
                {
                    int target_index = target_code.find(vocab[j]);
                    int code_index = code.find(vocab[j]);
                    if(target_index != std::string::npos)
                    {
                        //target code and the key code both contain the vocab[j]
                        if(code_index != std::string::npos)
                        {
                            if(pos.find(code[code_index-1]) != std::string::npos || pos.find(target_code[target_index-1]) != std::string::npos)
                                if(code[code_index-1] != target_code[target_index-1])
                                    diff_score--;
                            if(code_index != code.size()-1 && target_index != target_code.size()-1)
                                if(stat.find(code[code_index+1]) != std::string::npos || stat.find(target_code[target_index+1]) != std::string::npos)
                                    if(code[code_index+1] != target_code[target_index+1])
                                        diff_score--;
                        }
                        //target code contains the vocab[j] but the key code does not
                        else
                            diff_score--;
                    }
                    else
                    {
                        //the key code contains the vocab[j] but the target code does not
                        if(code_index != std::string::npos)
                            diff_score--;
                        //target code and the key code both do not contain the vocab[j]
                        else
                        {
                        }
                    }
                }
                if(diff_score > min_diff_score)
                    min_diff_score = diff_score;
                score_map[diff_score].push_back(SUGARNAMELOOKUP[i]);
            }
        }
        closest_matches = score_map[min_diff_score];
        ///SELECTING ONE MATCH FROM CLOSEST MATCHES
        ///RULE1: choose D over L isomer in closest matches
        ///RULE2: if the input chemical code has 'd's and the matches don't have corresponding 'd's. for D sugars choose the match that has '^' instead of 'd' for that index and '_' for L sugars
        Glycan::SugarName selected_match;
        selected_match.chemical_code_string_ = "";
        std::vector<size_t> positions; // holds all the positions that 'd' occurs within the chemical code stirng
        if(code.find("d") != std::string::npos)
        {
            size_t pos = code.find("d");
            while(pos != std::string::npos)
            {
                positions.push_back(pos);
                if(pos+1 < code.size())
                    pos = code.find("d",pos+1);
                else
                    break;
            }
        }
        int max_sugar_name_score = 0;
        std::string to_be_check_orientation = "";
        for(std::vector<Glycan::SugarName>::iterator it = closest_matches.begin(); it != closest_matches.end(); it++)
        {
            Glycan::SugarName name = (*it);
            int sugar_name_score = 0;
            to_be_check_orientation = "^";
            if(name.isomer_.compare("D") == 0)
                sugar_name_score = 100;
            else
                to_be_check_orientation = "_";

            if(positions.size() == 0)
            {
                selected_match = name;
                break;
            }
            else
            {
                std::string current = "";
                std::string non_deoxy_position = "";
                std::string deoxy_position = "";
                int count = 0;
                int corresponding_d_counts = 0;
                for(unsigned int i = 0; i < positions.size(); i++)
                {
                    if(positions.at(i) < name.chemical_code_string_.size())
                    {
                        deoxy_position = name.chemical_code_string_.at(positions.at(i));
                        if(deoxy_position.compare("d") == 0)
                            corresponding_d_counts = corresponding_d_counts + 2;
                        else
                        {
                            current = name.chemical_code_string_.at(positions.at(i));
                            if(current.compare("1") == 0) ///check for -1 and +1. in this case orientation is 2 positions back from the current position
                                non_deoxy_position = name.chemical_code_string_.at(positions.at(i) - 2);
                            else
                                non_deoxy_position = name.chemical_code_string_.at(positions.at(i) - 1);
                            if(non_deoxy_position.compare(to_be_check_orientation) == 0)
                                count++;
                        }
                    }
                }
                sugar_name_score += corresponding_d_counts + count;
                if(sugar_name_score > max_sugar_name_score)
                {
                    max_sugar_name_score = sugar_name_score;
                    selected_match = name;
                }
            }
        }
        if(selected_match.chemical_code_string_.compare("") == 0)
            selected_match = closest_matches.at(0);
        return selected_match;
    }

    /*! \fn
      * A function in order to look up the complex name of the sugar structure based on the given string version of the chemical code structure
      * @param code The string complex chemical code structure
      * @return COMPLEXSUGARNAMELOOKUP The matched row of the lookup table with the given code
      */
    inline Glycan::SugarName ComplexSugarNameLookup(std::string code)
    {
        for(int i = 0; i < COMPLEXSUGARNAMELOOKUPSIZE; i++)
        {
            if(code.compare(COMPLEXSUGARNAMELOOKUP[i].chemical_code_string_) == 0)
                return COMPLEXSUGARNAMELOOKUP[i];
        }
        return COMPLEXSUGARNAMELOOKUP[0];
    }    

    /*! \fn
      * A function in order to initializing the common terminal residue map
      * @return COMMON_TERMINAL_REDSIDUES The mapping between the terminal residue names and the terminal residue names
      */
    inline ResidueNameMap InitializeCommonTerminalResidueMap()
    {
        ResidueNameMap COMMON_TERMINAL_REDSIDUES = ResidueNameMap();
        COMMON_TERMINAL_REDSIDUES["ROH"] = "ROH";
        COMMON_TERMINAL_REDSIDUES["TBT"] = "TBT";
        COMMON_TERMINAL_REDSIDUES["OME"] = "OME";
        return COMMON_TERMINAL_REDSIDUES;
    }

    inline gmml::ResidueCodeName ResidueNameCodeLookup(std::string residue_name)
    {
        for(int i = 0; i < RESIDUENAMECODELOOKUPSIZE; i++)
        {
            if(residue_name.compare(RESIDUENAMECODELOOKUP[i].name_) == 0)
                return RESIDUENAMECODELOOKUP[i];
        }
        return RESIDUENAMECODELOOKUP[0];
    }

    inline gmml::ResidueCodeName ResidueCodeNameLookup(std::string residue_code)
    {
        for(int i = 0; i < RESIDUENAMECODELOOKUPSIZE; i++)
        {
            if(residue_code.compare(RESIDUENAMECODELOOKUP[i].code_) == 0)
                return RESIDUENAMECODELOOKUP[i];
        }
        return RESIDUENAMECODELOOKUP[0];
    }

    inline gmml::ResidueCodeName ResidueNameIndexLookup(std::string residue_name)
    {
        for(int i = 0; i < RESIDUENAMECODELOOKUPSIZE; i++)
        {
            if(residue_name.compare(RESIDUENAMECODELOOKUP[i].name_) == 0)
                return RESIDUENAMECODELOOKUP[i];
        }
        return RESIDUENAMECODELOOKUP[0];
    }

    inline gmml::AmberGlycamMap AmberGlycamLookup(std::string amber_residue_name)
    {
        for(int i = 0; i < AMBERGLYCAMLOOKUPSIZE; i++)
        {
            if(amber_residue_name.compare(AMBERGLYCAMLOOKUP[i].amber_name_) == 0)
                return AMBERGLYCAMLOOKUP[i];
        }
        return AMBERGLYCAMLOOKUP[0];
    }

    inline gmml::AmberGlycamMap GlycamAmberLookup(std::string glycam_residue_name)
    {
        for(int i = 0; i < AMBERGLYCAMLOOKUPSIZE; i++)
        {
            if(glycam_residue_name.compare(AMBERGLYCAMLOOKUP[i].glycam_name_) == 0)
                return AMBERGLYCAMLOOKUP[i];
        }
        return AMBERGLYCAMLOOKUP[0];
    }

    inline gmml::AtomTypesInfo AtomTypesLookup(std::string atom_type)
    {
        for(int i = 0; i < ATOMTYPESINFOLOOKUPSIZE; i++)
        {
            if(atom_type.compare(ATOMTYPESINFOLOOKUP[i].atom_type_) == 0)
                return ATOMTYPESINFOLOOKUP[i];
        }
        return ATOMTYPESINFOLOOKUP[0];
    }

    /*! \fn
      * A function in order to write the information/warning/error messages produced by the program into a log file
      * @param line The line number producing the message
      * @param file_path The file path of the file which the message has been produced within in
      * @param level The type of the produced message INF/WAR/ERR
      * @param msg The message content that has been produced
      * @param out_file_name The name of the output log file
      */
    inline void log(int line, std::string file_path, LogLevel level, std::string msg, std::string out_file_name = "log.log")
    {
//        std::ofstream file;
//        file.open(out_file_name.c_str(), std::ios_base::app);

//        time_t t = time(0);
//        std::string time_str = std::asctime(std::localtime(&t));
//        file << time_str.substr(0, time_str.size() - 1) << " >>> " << file_path << ":" << line << " >>>";
//        switch(level)
//        {
//            case INF:
//                file << " [INFO]: ";
//                break;
//            case ERR:
//                file << " [ERROR]: ";
//                break;
//            case WAR:
//                file << " [WARNING]: ";
//                break;
//        }
//        file << msg << std::endl;

//        file.close();
    }

    inline double** GenerateRotationMatrix(GeometryTopology::Coordinate* direction, GeometryTopology::Coordinate* parent, double angle)
    {
        double** rotation_matrix = new double*[3];
        for(int i = 0; i < 3; i++)
            rotation_matrix[i] = new double[4];

        direction->Normalize();
        double u = direction->GetX();
        double v = direction->GetY();
        double w = direction->GetZ();

        double a = parent->GetX();
        double b = parent->GetY();
        double c = parent->GetZ();

        double u2 = u*u;
        double v2 = v*v;
        double w2 = w*w;
        double cos_rotation_angle = cos(angle);
        double sin_rotation_angle = sin(angle);

        rotation_matrix[0][3] = a * (v2 + w2) - u * (b * v + c * w) + (u * (b * v + c * w) - a * (v2 + w2)) * cos_rotation_angle + (b * w - c * v) * sin_rotation_angle;
        rotation_matrix[1][3] = b * (u2 + w2) - v * (a * u + c * w) + (v * (a * u + c * w) - b * (u2 + w2)) * cos_rotation_angle + (c * u - a * w) * sin_rotation_angle;
        rotation_matrix[2][3] = c * (u2 + v2) - w * (a * u + b * v) + (w * (a * u + b * v) - c * (u2 + v2)) * cos_rotation_angle + (a * v - b * u) * sin_rotation_angle;

        rotation_matrix[0][0] = u2 + (v2 + w2) * cos_rotation_angle;
        rotation_matrix[0][1] = u * v * (1 - cos_rotation_angle) - w * sin_rotation_angle;
        rotation_matrix[0][2] = u * w * (1 - cos_rotation_angle) + v * sin_rotation_angle;

        rotation_matrix[1][0] = u * v * (1 - cos_rotation_angle) + w * sin_rotation_angle;
        rotation_matrix[1][1] = v2 + (u2 + w2) * cos_rotation_angle;
        rotation_matrix[1][2] = v * w * (1 - cos_rotation_angle) - u * sin_rotation_angle;

        rotation_matrix[2][0] = u * w * (1 - cos_rotation_angle) - v * sin_rotation_angle;
        rotation_matrix[2][1] = v * w * (1 - cos_rotation_angle) + u * sin_rotation_angle;
        rotation_matrix[2][2] = w2 + (u2 + v2) * cos_rotation_angle;

        return rotation_matrix;
    }
    /*! \fn
      * A function in order to replace all occurrences of a sub-string with another sub-string in a string
      * @param str The string that is going to be manipulated
      * @param search The sub-string that is going to be searched and later replaced in the string
      * @param replace The replacement sub-string
      */
    inline void FindReplaceString(std::string &str, std::string search, std::string replace)
    {
        for(size_t pos = 0; ; pos += replace.length())
        {
            pos = str.find( search, pos);
            if(pos == std::string::npos) break;

            str.erase(pos, search.length());
            str.insert(pos, replace);
        }
    }

    inline std::string ConvertVectorString2String(std::vector<std::string> vector_string)
    {
        std::string result = "";
        for(unsigned int i = 0; i < vector_string.size(); i++)
        {
            result += vector_string.at(i);
        }
        return result;
    }
}


#endif // UTILS_HPP
