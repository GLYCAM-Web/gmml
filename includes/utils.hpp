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
#include "MolecularModeling/atom.hpp"

#include <fstream>

#include <iostream>

namespace gmml
{
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
      * Removes duplicate spaces inside of the string.
      * @param str String with duplicate spaces
      * @return Given string without duplicate spaces within the original one
      */
    inline std::string& TrimSpaces(std::string& str)
    {
          std::string s;
          bool first = true;
          bool space = false;
          std::string::iterator iter;
          for(iter = str.begin(); iter != str.end(); ++iter){
              if(*iter == ' '){
                  if(first == false){
                      space = true;
                  }
              }else{
                  if(*iter != ',' && *iter != '.'){
                      if(space){
                          s.push_back(' ');
                      }
                  }
                  s.push_back(*iter);
                  space = false;
                  first = false;
              }
          }
          str = s;
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
        // if(ss >> val)
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
    inline std::string execcommand(const char* cmd) {
        char buffer[128];
        std::string result = "";
        FILE* pipe = popen(cmd, "r");
        if (!pipe) throw std::runtime_error("popen() failed!");
        try {
            while (!feof(pipe)) {
                if (fgets(buffer, 128, pipe) != NULL)
                    result += buffer;
            }
        } catch (...) {
            pclose(pipe);
            throw;
        }
        pclose(pipe);
        return result;
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
      * A function in order to look up the stereochemistry name of the sugar structure based on the given string version of the chemical code structure
      * @param code The string chemical code structure
      * @return SUGARNAMELOOKUP The matched row of the lookup table with the given code
      */
    inline Glycan::SugarName SugarStereoChemistryNameLookup(std::string code)
    {
        for(int i = 0; i < SUGARNAMELOOKUPSIZE; i++)
        {
            if(code.compare(SUGARNAMELOOKUP[i].chemical_code_string_) == 0){
                return SUGARNAMELOOKUP[i];
            }
        }
        return SUGARNAMELOOKUP[0];
    }
    
    inline Glycan::SugarName ResidueSugarNameLookup(std::string residue)
    {
        for(int i = 0; i < SUGARNAMELOOKUPSIZE; i++)
        {
            if(residue.compare(SUGARNAMELOOKUP[i].pdb_code_) == 0){
                return SUGARNAMELOOKUP[i];
            }
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
        std::string vocab[] = {"2", "3", "4", "a", "+1", "+2", "+3", "-1", "NAc"};
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
                    size_t target_index = target_code.find(vocab[j]);
                    size_t code_index = code.find(vocab[j]);
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
    
    inline Glycan::SugarName ResidueComplexSugarNameLookup(std::string residue)
    {
        for(int i = 0; i < COMPLEXSUGARNAMELOOKUPSIZE; i++)
        {
            if(residue.compare(COMPLEXSUGARNAMELOOKUP[i].pdb_code_) == 0)
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

    inline gmml::AminoacidGlycamMap AminoacidGlycamLookup(std::string aminoacid_residue_name)
    {
        for(int i = 0; i < AMINOACIDGLYCAMLOOKUPSIZE; i++)
        {
            if(aminoacid_residue_name.compare(AMINOACIDGLYCAMLOOKUP[i].aminoacid_name_) == 0)
                return AMINOACIDGLYCAMLOOKUP[i];
        }
        return AMINOACIDGLYCAMLOOKUP[0];
    }

    inline gmml::AminoacidGlycamMap GlycamAminoacidLookup(std::string glycam_residue_name)
    {
        for(int i = 0; i < AMINOACIDGLYCAMLOOKUPSIZE; i++)
        {
            if(glycam_residue_name.compare(AMINOACIDGLYCAMLOOKUP[i].glycam_name_) == 0)
                return AMINOACIDGLYCAMLOOKUP[i];
        }
        return AMINOACIDGLYCAMLOOKUP[0];
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
    inline void log(int line, std::string file_path, LogLevel level, std::string msg, std::string out_file_name = "")
    {
      std::ofstream file;
      if(out_file_name == "")
      {
        std::string GEMSHOME_ERROR = "\nMust set GEMSHOME environment variable.\n\n    BASH:   export GEMSHOME=/path/to/gems\n    SH:     setenv GEMSHOME /path/to/gems\n";
        char* gemshome_env_var = std::getenv("GEMSHOME");
        // Check if the environment variables exist.
        if(!gemshome_env_var) 
        {
          std::cout << GEMSHOME_ERROR << std::endl;
        }
        std::string GEMSHOME(gemshome_env_var);
        out_file_name = GEMSHOME + "/log.log";
      }
      file.open(out_file_name.c_str(), std::ios_base::app);

       time_t t = time(0);
       std::string time_str = std::asctime(std::localtime(&t));
       file << time_str.substr(0, time_str.size() - 1) << " >>> " << file_path << ":" << line << " >>>";
       switch(level)
       {
           case INF:
               file << " [INFO]: ";
               break;
           case ERR:
               file << " [ERROR]: ";
               break;
           case WAR:
               file << " [WARNING]: ";
               break;
       }
       file << msg << std::endl;

       file.close();
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

    inline void AddTriple(std::string s, std::string p, std::string o, std::stringstream& stream)
    {
        // Find and replace escape characters to and '\'
        stream << s << " " << p << " " << o << "." << std::endl;
    }

    inline void AddLiteral(std::string s, std::string p, std::string o, std::stringstream& stream)
    {
        // Find and replace
        stream << s << " " << p << " \"\"\"" << o << "\"\"\"." << std::endl;
    }

    inline void AddDecimal(std::string s, std::string p, float o, std::stringstream& stream)
    {
        stream << s << " " << p << " \"" << o << "\"^^xsd:decimal." << std::endl;
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
