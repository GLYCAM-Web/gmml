#include <iomanip>
#include <iostream>

#include "../../../includes/common.hpp"
#include "../../../includes/ParameterSet/ParameterFileSpace/parameterfileangle.hpp"
using namespace gmml;
using namespace ParameterFileSpace;

/////////////////////////////// CONSTRUCTOR ////////////////////////////////
ParameterFileAngle::ParameterFileAngle()
    : types_(), force_constant_(kNotSet), angle_(kNotSet), dscr_("") {}

ParameterFileAngle::ParameterFileAngle(const std::vector<std::string> &types, double force_constant, double angle, const std::string &dscr)
    : types_(types), force_constant_(force_constant), angle_(angle), dscr_(dscr) {}

////////////////////////// DISPLAY FUNCTION ////////////////////////////////
void ParameterFileAngle::Print(std::ostream& out)
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
    if (force_constant_ == kNotSet)
        out << std::setw(15) << "--";
    else
        out << std::setw(15) << force_constant_;

    if(angle_ == kNotSet)
        out << std::setw(10) << "--";
    else
        out << std::setw(10) << angle_;

    out << std::setw(60) << dscr_
        << std::endl;
}

