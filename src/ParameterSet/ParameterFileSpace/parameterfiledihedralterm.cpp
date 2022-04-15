#include <iomanip>
#include <iostream>

#include "../../../includes/common.hpp"
#include "../../../includes/ParameterSet/ParameterFileSpace/parameterfiledihedralterm.hpp"

using ParameterFileSpace::ParameterFileDihedralTerm;

//////////////////////////////////////////////////////////
//                       Constructor                    //
//////////////////////////////////////////////////////////
ParameterFileDihedralTerm::ParameterFileDihedralTerm() :
    factor_(gmml::dNotSet), force_constant_(gmml::dNotSet), phase_(gmml::dNotSet), periodicity_(gmml::dNotSet), dscr_("") {}

ParameterFileDihedralTerm::ParameterFileDihedralTerm(int factor, double force_constant, double phase, double periodicity, const std::string& dscr) :
    factor_(factor), force_constant_(force_constant), phase_(phase), periodicity_(periodicity), dscr_(dscr) {}

//////////////////////////////////////////////////////////
//                           ACCESSOR                   //
//////////////////////////////////////////////////////////

int ParameterFileDihedralTerm::GetFactor(){
    return factor_;
}

double ParameterFileDihedralTerm::GetForceConstant(){
    return force_constant_;
}

double ParameterFileDihedralTerm::GetPhase(){
    return phase_;
}

double ParameterFileDihedralTerm::GetPeriodicity(){
    return periodicity_;
}

std::string ParameterFileDihedralTerm::GetDscr(){
    return dscr_;
}

//////////////////////////////////////////////////////////
//                           MUTATOR                    //
//////////////////////////////////////////////////////////

void ParameterFileDihedralTerm::SetFactor(int factor){
    factor_ = factor;
}

void ParameterFileDihedralTerm::SetForceConstant(double force_constant){
    force_constant_ = force_constant;
}

void ParameterFileDihedralTerm::SetPhase(double phase){
    phase_ = phase;
}

void ParameterFileDihedralTerm::SetPeriodicity(double periodicity){
    periodicity_ = periodicity;
}

void ParameterFileDihedralTerm::SetDscr(const std::string dscr){
    dscr_ = dscr;
}

//////////////////////////////////////////////////////////
//                     DISPLAY FUNCTIONS                //
//////////////////////////////////////////////////////////
void ParameterFileDihedralTerm::Print(std::ostream& out)
{
    if(factor_ == gmml::iNotSet)
        out << std::setw(10) << "--";
    else
        out << std::setw(10) << factor_;

    if(force_constant_ == gmml::dNotSet)
        out << std::setw(10) << "--";
    else
        out << std::setw(10) << force_constant_;

    if(phase_ == gmml::dNotSet)
        out << std::setw(10) << "--";
    else
        out << std::setw(10) << phase_;

    if(periodicity_ == gmml::dNotSet)
        out << std::setw(12) << "--";
    else
        out << std::setw(12) << periodicity_;

    out << std::setw(60) << dscr_
        << std::endl;
}
