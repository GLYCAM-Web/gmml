#include <iomanip>
#include <iostream>

#include "../../../includes/ParameterSet/ParameterFileSpace/parameterfilebond.hpp"
#include "../../../includes/common.hpp"

using ParameterFileSpace::ParameterFileBond;

///////////////////////////////////////////////////////////
//                       Constructor                    //
//////////////////////////////////////////////////////////
ParameterFileBond::ParameterFileBond()
    : types_(), force_constant_(gmml::dNotSet), length_(gmml::dNotSet), dscr_("") {}

ParameterFileBond::ParameterFileBond(const std::vector<std::string> &types, double force_constant, double length, const std::string &dscr)
    : types_(types), force_constant_(force_constant), length_(length), dscr_(dscr) {}

ParameterFileBond::ParameterFileBond(const std::vector<std::string> &types, double force_constant, double length,
                                     const std::vector<double> &hbond_coefficients, const std::string &dscr)
    : types_(types), force_constant_(force_constant), length_(length), dscr_(dscr), hbond_coefficients_(hbond_coefficients) {}

//////////////////////////////////////////////////////////
//                           ACCESSOR                   //
//////////////////////////////////////////////////////////
/// Return the list of atom types involved in the current bond object
std::vector<std::string> ParameterFileBond::GetTypes()
{
    return types_;
}

/// Return force constant
double ParameterFileBond::GetForceConstant()
{
    return force_constant_;
}

/// Return angle
double ParameterFileBond::GetLength()
{
    return length_;
}

/// Return dscr
std::string ParameterFileBond::GetDscr()
{
    return dscr_;
}

/// Return hbonded_coefficients
std::vector<double> ParameterFileBond::GetHbondCoefficients()
{
    return hbond_coefficients_;
}

//////////////////////////////////////////////////////////
//                           MUTATOR                    //
//////////////////////////////////////////////////////////
/// Set the list of atom types involved in the current bond object
void ParameterFileBond::SetTypes(std::vector<std::string> types)
{
    types_.clear();
    for(std::vector<std::string>::iterator it = types.begin(); it != types.end(); it++)
    {
        types_.push_back(*it);
    }
}

/// Set force constant
void ParameterFileBond::SetForceConstant(double force_constant)
{
    force_constant_ = force_constant;
}

/// Set angle
void ParameterFileBond::SetLength(double length){
    length_ = length;
}

/// Set dscr
void ParameterFileBond::SetDscr(const std::string dscr){
    dscr_ = dscr;
}

void ParameterFileBond::SetHbondCoefficients(std::vector<double> hbond_coefficients){
    hbond_coefficients_.clear();
    for(std::vector<double>::iterator it = hbond_coefficients.begin(); it != hbond_coefficients.end(); it++)
    {
        hbond_coefficients_.push_back(*it);
    }
}

//////////////////////////////////////////////////////////
//                     DISPLAY FUNCTIONS                //
//////////////////////////////////////////////////////////
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

    if(force_constant_ == gmml::dNotSet)
        out << std::setw(15) << "--";
    else
        out << std::setw(15) << force_constant_;

    if(length_ == gmml::dNotSet)
        out << std::setw(10) << "--";
    else
        out << std::setw(10) << length_;

    out << std::setw(60) << dscr_
        << std::endl;
}
