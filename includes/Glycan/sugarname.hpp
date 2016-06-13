#ifndef SUGARNAME_HPP
#define SUGARNAME_HPP

#include <string>

namespace Glycan
{
    struct SugarName {
            std::string chemical_code_string_;                          /*!< The string version of the chemical code structure >*/
            std::string monosaccharide_stereochemistry_name_;           /*!< The stereochemistry name of the monosacchride >*/
            std::string monosaccharide_stereochemistry_short_name_;     /*!< The condensed stereochemisty name of the monosacchride >*/
            std::string isomer_;                                        /*!< The isomer(D/L) info of the monosacchride >*/
            std::string name_;                                          /*!< The base name of the monosacchride >*/
            std::string ring_type_;                                     /*!< The ring type of the monosacchride which can be P/F which represents pyranose and furanose accordingly >*/
            std::string configuration_;                                 /*!< The configuration(a/b, alpha/beta) info of the monosacchride >*/
            std::string monosaccharide_name_;                           /*!< The name of the monosacchride >*/
            std::string monosaccharide_short_name_;                     /*!< The condensed name of the monosacchride >*/
            std::string pdb_code_;                     /*!< The condensed name of the monosacchride >*/
    } ;
}

#endif // SUGARNAME_HPP
