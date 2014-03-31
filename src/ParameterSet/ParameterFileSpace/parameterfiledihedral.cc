#include <iomanip>
#include <iostream>

#include "../../../includes/ParameterSet/ParameterFileSpace/parameterfiledihedral.hpp"
#include "../../../includes/ParameterSet/ParameterFileSpace/parameterfiledihedralterm.hpp"
#include "../../../includes/common.hpp"

using namespace gmml;
using namespace ParameterFileSpace;

//////////////////////////////////////////////////////////
//                       Constructor                    //
//////////////////////////////////////////////////////////
ParameterFileDihedral::ParameterFileDihedral()
    : types_(), terms_(), scee_(dNotSet), scnb_(dNotSet), is_generic_(false), is_improper_(false) {}

ParameterFileDihedral::ParameterFileDihedral(const std::vector<std::string> &types, const ParameterFileDihedralTerm& term,
                                             double scee, double scnb, bool is_generic, bool is_improper)
    : types_(types), terms_(), scee_(scee), scnb_(scnb), is_generic_(is_generic), is_improper_(is_improper)
{
    terms_.push_back(term);
}

///Delaram
//////////////////////////////////////////////////////////
//                           ACCESSOR                   //
//////////////////////////////////////////////////////////

std::vector<std::string> ParameterFileDihedral::GetTypes()
{
    return types_;
}

std::vector<ParameterFileDihedralTerm> ParameterFileDihedral::GetTerms()
{
    return terms_;
}

double ParameterFileDihedral::GetScee()
{
    return scee_;
}

double ParameterFileDihedral::GetScnb()
{
    return scnb_;
}

bool ParameterFileDihedral::GetIsGeneric()
{
    return is_generic_;
}

bool ParameterFileDihedral::GetIsImproper()
{
    return is_improper_;
}

//////////////////////////////////////////////////////////
//                           MUTATOR                    //
//////////////////////////////////////////////////////////

void ParameterFileDihedral::SetTypes(const std::vector<std::string> types){
    types_ = types;
}

void ParameterFileDihedral::SetTerms(const std::vector<ParameterFileDihedralTerm> terms){
    terms_ = terms;
}

void ParameterFileDihedral::AddTerm(ParameterFileDihedralTerm term){
    terms_.push_back(term);
}

void ParameterFileDihedral::SetScee(double scee){
    scee_ = scee;
}

void ParameterFileDihedral::SetScnb(double scnb){
    scnb_ = scnb;
}

void ParameterFileDihedral::SetIsGeneric(bool is_generic){
    is_generic_ = is_generic;
}

void ParameterFileDihedral::SetIsImproper(bool is_improper){
    is_improper_ = is_improper;
}

///Delaram

//////////////////////////////////////////////////////////
//                     DISPLAY FUNCTIONS                //
//////////////////////////////////////////////////////////
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

    if(scee_ == dNotSet)
        out << std::setw(6) << "--";
    else
        out << std::setw(6) << scee_;

    if(scnb_ == dNotSet)
        out << std::setw(6) << "--";
    else
        out << std::setw(6) << scnb_;

    for(std::vector<ParameterFileDihedralTerm>::iterator it = terms_.begin(); it != terms_.end(); it++)
    {
        if(it == terms_.begin())
            it->Print(out);
        else
        {
            out << std::setw(54) << "";
            it -> Print(out);
        }
    }
}
