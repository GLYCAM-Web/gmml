#include <iomanip>
#include <iostream>

#include "../../../includes/ParameterSet/ParameterFileSpace/parameterfilebond.hpp"
#include "../../../includes/common.hpp"

using namespace ParameterFileSpace;
using namespace gmml;

/////////////////////////////// CONSTRUCTOR ////////////////////////////////
ParameterFileBond::ParameterFileBond()
    : types_(), force_constant_(kNotSet), length_(kNotSet), dscr_("") {}

ParameterFileBond::ParameterFileBond(const std::vector<std::string> &types, double force_constant, double length, const std::string &dscr)
    : types_(types), force_constant_(force_constant), length_(length), dscr_(dscr) {}

ParameterFileBond::ParameterFileBond(const std::vector<std::string> &types, double force_constant, double length,
                                     const std::vector<double> &hbond_coefficients, const std::string &dscr)
    : types_(types), force_constant_(force_constant), length_(length), dscr_(dscr), hbond_coefficients_(hbond_coefficients) {}

////////////////////////// DISPLAY FUNCTION ////////////////////////////////
void ParameterFileBond::Print(std::ostream& out)
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

    if(force_constant_ == kNotSet)
        out << std::setw(15) << "--";
    else
        out << std::setw(15) << force_constant_;

    if(length_ == kNotSet)
        out << std::setw(10) << "--";
    else
        out << std::setw(10) << length_;

    out << std::setw(60) << dscr_
        << std::endl;
}
