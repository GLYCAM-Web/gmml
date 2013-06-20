#include <iomanip>
#include <iostream>

#include "../../../includes/ParameterSet/ParameterFileSpace/parameterfiledihedral.hpp"
#include "../../../includes/ParameterSet/ParameterFileSpace/parameterfiledihedralterm.hpp"
#include "../../../includes/common.hpp"

using namespace gmml;
using namespace ParameterFileSpace;

/////////////////////////////// CONSTRUCTOR ////////////////////////////////
ParameterFileDihedral::ParameterFileDihedral()
    : types_(), terms_(), scee_(kNotSet), scnb_(kNotSet), is_generic_(false), is_improper_(false) {}

ParameterFileDihedral::ParameterFileDihedral(const std::vector<std::string> &types, const ParameterFileDihedralTerm& term,
                                             double scee, double scnb, bool is_generic, bool is_improper)
    : types_(types), terms_(), scee_(scee), scnb_(scnb), is_generic_(is_generic), is_improper_(is_improper)
{
    terms_.push_back(term);
}

////////////////////////// DISPLAY FUNCTION ////////////////////////////////
void ParameterFileDihedral::Print(std::ostream& out)
{
    for(std::vector<std::string>::iterator it = types_.begin(); it != types_.end(); it++)
    {
        if(it != types_.end() - 1)
        {
            out << std::setw(4) << (*it)
                << std::setw(2) << "-";
        }
        else
            out << std::setw(4) << (*it);
    }
    if(is_generic_ == true)
        out << std::setw(10) << "YES";
    else
        out << std::setw(10) << "NO";
    if(is_improper_ == true)
        out << std::setw(10) << "YES";
    else
        out << std::setw(10) << "NO";

    if(scee_ == kNotSet)
        out << std::setw(6) << "--";
    else
        out << std::setw(6) << scee_;

    if(scnb_ == kNotSet)
        out << std::setw(6) << "--";
    else
        out << std::setw(6) << scnb_;

    for(std::vector<ParameterFileDihedralTerm>::iterator it = terms_.begin(); it != terms_.end(); it++)
    {
        if(it == terms_.begin())
            it->Print(out);
        else
        {
            out << std::setw(42) << "";
            it -> Print(out);
        }
    }
}
