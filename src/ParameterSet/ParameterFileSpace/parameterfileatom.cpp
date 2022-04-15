#include <iomanip>
#include <iostream>

#include "../../../includes/ParameterSet/ParameterFileSpace/parameterfileatom.hpp"
#include "../../../includes/common.hpp"

using ParameterFileSpace::ParameterFileAtom;

//////////////////////////////////////////////////////////
//                       Constructor                    //
//////////////////////////////////////////////////////////
ParameterFileAtom::ParameterFileAtom()
    : type_(""), mass_(gmml::dNotSet), polarizability_(gmml::dNotSet), dscr_(""), radius_(gmml::dNotSet), well_depth_(gmml::dNotSet), mod4_dscr_("") {}

ParameterFileAtom::ParameterFileAtom(const std::string &type, double mass, double polarizability, const std::string &dscr)
    : type_(type), mass_(mass), polarizability_(polarizability),dscr_(dscr), radius_(gmml::dNotSet), well_depth_(gmml::dNotSet), is_hydrophilic_(false) {}

ParameterFileAtom::ParameterFileAtom(const std::string &type, double mass, double polarizability, double radius,
                                     double well_depth, const std::string &dscr, const std::string &mod4_dscr)
    : type_(type), mass_(mass), polarizability_(polarizability), dscr_(dscr), radius_(radius), well_depth_(well_depth), mod4_dscr_(mod4_dscr) {}

ParameterFileAtom::ParameterFileAtom(const std::string &type, double mass, double polarizability, double radius,
                                     double well_depth, const std::vector<std::string> &equivalent_list, const std::string &dscr, const std::string &mod4_dscr, bool is_hydrophilic)
    : type_(type), mass_(mass), polarizability_(polarizability), dscr_(dscr), radius_(radius), well_depth_(well_depth), mod4_dscr_(mod4_dscr),
      is_hydrophilic_(is_hydrophilic), equivalent_list_(equivalent_list) {}

//////////////////////////////////////////////////////////
//                           ACCESSOR                   //
//////////////////////////////////////////////////////////
/// Retuns type
std::string ParameterFileAtom::GetType()
{
    return type_;
}

/// Returns mass
double ParameterFileAtom::GetMass()
{
    return mass_;
}


/// Returns polarizability
double ParameterFileAtom::GetPolarizability()
{
    return polarizability_;
}

/// Return dscr
std::string ParameterFileAtom::GetDscr()
{
    return dscr_;
}

/// Returns Radius
double ParameterFileAtom::GetRadius()
{
    return radius_;
}

/// Retunrs well debpth
double ParameterFileAtom::GetWellDepth()
{
    return well_depth_;
}

/// Returns mod4 dscr
std::string ParameterFileAtom::GetMod4Dscr()
{
    return mod4_dscr_;
}

/// Returns is_hydrophilic
bool ParameterFileAtom::GetIsHydrophilic()
{
    return is_hydrophilic_;
}

/// Returns equivalent_list
std::vector<std::string> ParameterFileAtom::GetEquivalentList()
{
    return equivalent_list_;
}

//////////////////////////////////////////////////////////
//                           MUTATOR                    //
//////////////////////////////////////////////////////////

void ParameterFileAtom::SetType( std::string type){
    type_ = type;
}

void ParameterFileAtom::SetMass(double mass){
    mass_ = mass;
}

void ParameterFileAtom::SetPolarizability(double polarizability){
    polarizability_ = polarizability;
}

void ParameterFileAtom::SetDscr( std::string dscr){
    dscr_ = dscr;
}

void ParameterFileAtom::SetRadius(double radius){
    radius_ = radius;
}

void ParameterFileAtom::SetWellDepth(double well_depth){
    well_depth_ = well_depth;
}

void ParameterFileAtom::SetMod4Dscr( std::string mod4_dscr)
{
    mod4_dscr_ = mod4_dscr;
}

void ParameterFileAtom::SetIsHydrophilic(bool is_hydrophilic){
    is_hydrophilic_ = is_hydrophilic;
}

void ParameterFileAtom::SetEquivalentList( std::vector<std::string> equivalent_list){
    equivalent_list_.clear();
    for(std::vector<std::string>::iterator it = equivalent_list.begin(); it != equivalent_list.end(); it++)
    {
        equivalent_list_.push_back(*it);
    }
}

//////////////////////////////////////////////////////////
//                     DISPLAY FUNCTIONS                //
//////////////////////////////////////////////////////////
void ParameterFileAtom::Print(std::ostream& out)
{
    out << std::setw(6) << type_;
    if(mass_ == gmml::dNotSet)
        out << std::setw(6) << "--";
    else
        out << std::setw(6) << mass_;

    if(polarizability_ == gmml::dNotSet)
        out << std::setw(16) << "--";
    else
        out << std::setw(16) << polarizability_;

    if(radius_ == gmml::dNotSet)
        out << std::setw(8) << "--";
    else
        out << std::setw(8) << radius_;

    if(well_depth_ == gmml::dNotSet)
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
