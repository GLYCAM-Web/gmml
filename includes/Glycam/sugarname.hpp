#ifndef SUGARNAME_HPP
#define SUGARNAME_HPP

#include <string>

namespace Glycam
{
    struct SugarName {
            std::string chemical_code_string_;
            std::string monosaccharide_stereochemistry_name_;
            std::string monosaccharide_stereochemistry_short_name_;
            std::string isomer_;
            std::string name_;
            std::string ring_type_;
            std::string configuration_;
            std::string monosaccharide_name_;
            std::string monosaccharide_short_name_;
    } ;
}

#endif // SUGARNAME_HPP
