#include <iomanip>
#include <iostream>

#include "../../../includes/ParameterSet/ParameterFileSpace/parameterfileatom.hpp"
#include "../../../includes/common.hpp"

using namespace ParameterFileSpace;
using namespace gmml;

/////////////////////////////// CONSTRUCTOR ////////////////////////////////
ParameterFileAtom::ParameterFileAtom()
    : type_(""), mass_(kNotSet), polarizability_(kNotSet), dscr_(""), radius_(kNotSet), well_depth_(kNotSet), mod4_dscr_("") {}

ParameterFileAtom::ParameterFileAtom(const std::string &type, double mass, double polarizability, const std::string &dscr)
    : type_(type), mass_(mass), polarizability_(polarizability),dscr_(dscr), radius_(kNotSet), well_depth_(kNotSet) {}

ParameterFileAtom::ParameterFileAtom(const std::string &type, double mass, double polarizability, double radius,
                                     double well_depth, const std::string &dscr, const std::string &mod4_dscr)
    : type_(type), mass_(mass), polarizability_(polarizability), dscr_(dscr), radius_(radius), well_depth_(well_depth), mod4_dscr_(mod4_dscr) {}

ParameterFileAtom::ParameterFileAtom(const std::string &type, double mass, double polarizability, double radius,
                                     double well_depth, const std::vector<std::string> &equivalent_list, const std::string &dscr, const std::string &mod4_dscr, bool is_hydrophilic)
    : type_(type), mass_(mass), polarizability_(polarizability), dscr_(dscr), radius_(radius), well_depth_(well_depth), mod4_dscr_(mod4_dscr),
      is_hydrophilic_(is_hydrophilic), equivalent_list_(equivalent_list) {}

////////////////////////// DISPLAY FUNCTION ////////////////////////////////
void ParameterFileAtom::Print(std::ostream& out)
{
    out << std::setw(6) << type_;
    if(mass_ == kNotSet)
        out << std::setw(6) << "--";
    else
        out << std::setw(6) << mass_;

    if(polarizability_ == kNotSet)
        out << std::setw(16) << "--";
    else
        out << std::setw(16) << polarizability_;

    if(radius_ == kNotSet)
        out << std::setw(8) << "--";
    else
        out << std::setw(8) << radius_;

    if(well_depth_ == kNotSet)
        out << std::setw(12) << "--";
    else
        out << std::setw(12) << well_depth_;

    if(is_hydrophilic_ == true)
        out << std::setw(13) << "YES";
    else
        out << std::setw(13) << "NO";

    out << std::setw(60) << dscr_
        << std::setw(60) << mod4_dscr_
        << std::endl;

    if(!equivalent_list_.empty())
    {
        out << std::setw(20) << "" << "Equivalent symbols: ";
        for(std::vector<std::string>::iterator it = equivalent_list_.begin(); it != equivalent_list_.end(); it++)
        {
            out << (*it) << "   ";
        }
        out << std::endl;
    }
}
