#include <iomanip>
#include <iostream>

#include "../../../includes/common.hpp"
#include "../../../includes/ParameterSet/ParameterFileSpace/parameterfileangle.hpp"

using ParameterFileSpace::ParameterFileAngle;

//////////////////////////////////////////////////////////
//                       Constructor                    //
//////////////////////////////////////////////////////////
ParameterFileAngle::ParameterFileAngle()
    : types_(), force_constant_(gmml::dNotSet), angle_(gmml::dNotSet), dscr_("") {}

ParameterFileAngle::ParameterFileAngle(const std::vector<std::string> &types, double force_constant, double angle, const std::string &dscr)
    : types_(types), force_constant_(force_constant), angle_(angle), dscr_(dscr) {}

//////////////////////////////////////////////////////////
//                           ACCESSOR                   //
//////////////////////////////////////////////////////////
/// Return the list of atom types involved in the current angle object
std::vector<std::string> ParameterFileAngle::GetTypes()
{
    return types_;
}

/// Return force constant
double ParameterFileAngle::GetForceConstant()
{
    return force_constant_;
}

/// Return angle
double ParameterFileAngle::GetAngle()
{
    return angle_;
}

/// Return dscr
std::string ParameterFileAngle::GetDscr()
{
    return dscr_;
}

//////////////////////////////////////////////////////////
//                           MUTATOR                    //
//////////////////////////////////////////////////////////
/// Set the list of atom types involved in the current angle object
void ParameterFileAngle::SetTypes(std::vector<std::string> types)
{
    types_.clear();
    for(std::vector<std::string>::iterator it = types.begin(); it != types.end(); it++)
    {
        types_.push_back(*it);
    }
}

/// Set force constant
void ParameterFileAngle::SetForceConstant(double force_constant)
{
    force_constant_ = force_constant;
}

/// Set angle
void ParameterFileAngle::SetAngle(double angle){
    angle_ = angle;
}

/// Set dscr
void ParameterFileAngle::SetDscr(const std::string dscr){
    dscr_ = dscr;
}

//////////////////////////////////////////////////////////
//                     DISPLAY FUNCTIONS                //
//////////////////////////////////////////////////////////
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
    if (force_constant_ == gmml::dNotSet)
        out << std::setw(15) << "--";
    else
        out << std::setw(15) << force_constant_;

    if(angle_ == gmml::dNotSet)
        out << std::setw(10) << "--";
    else
        out << std::setw(10) << angle_;

    out << std::setw(60) << dscr_
        << std::endl;
}
