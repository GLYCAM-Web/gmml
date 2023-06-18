#ifndef SUGARNAME_HPP
#define SUGARNAME_HPP

#include <string>
#include <sstream>

namespace Glycan
{
    struct SugarName
    {
        std::string chemical_code_string_;                /*!< The string version of the chemical code structure >*/
        std::string monosaccharide_stereochemistry_name_; /*!< The stereochemistry name of the monosaccharide >*/
        std::string
            monosaccharide_stereochemistry_short_name_; /*!< The condensed stereochemisty name of the monosaccharide >*/
        std::string isomer_;                            /*!< The isomer(D/L) info of the monosaccharide >*/
        std::string name_;                              /*!< The base name of the monosaccharide >*/
        std::string ring_type_; /*!< The ring type of the monosaccharide which can be P/F which represents pyranose and
                                   furanose accordingly >*/
        std::string configuration_;             /*!< The configuration(a/b, alpha/beta) info of the monosaccharide >*/
        std::string monosaccharide_name_;       /*!< The name of the monosaccharide >*/
        std::string monosaccharide_short_name_; /*!< The condensed name of the monosaccharide >*/
        std::string pdb_code_;                  /*!< The condensed name of the monosaccharide >*/

        /** A function to convert the sugar name structure into a string
         * @return ss The string version of the sugar name structure
         */
        std::string toString()
        {
            std::stringstream ss;

            ss << "Chemical code: ";
            ss << chemical_code_string_;
            ss << "\n";

            ss << "Stereochemistry name: ";
            ss << monosaccharide_stereochemistry_name_;
            ss << "\n";

            ss << "Stereochemistry short name: ";
            ss << monosaccharide_stereochemistry_short_name_;
            ss << "\n";

            ss << "Isomer: ";
            ss << isomer_;
            ss << "\n";

            ss << "Base name: ";
            ss << name_;
            ss << "\n";

            ss << "Ring type: ";
            ss << ring_type_;
            ss << "\n";

            ss << "Configuration: ";
            ss << configuration_;
            ss << "\n";

            ss << "Name: ";
            ss << monosaccharide_name_;
            ss << "\n";

            ss << "Short name: ";
            ss << monosaccharide_short_name_;
            ss << "\n";

            ss << "PDB code: ";
            ss << pdb_code_;
            ss << "\n";

            return ss.str();
        }
    };
} // namespace Glycan

#endif // SUGARNAME_HPP
