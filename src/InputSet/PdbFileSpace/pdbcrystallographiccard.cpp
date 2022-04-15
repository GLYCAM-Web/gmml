#include "../../../includes/InputSet/PdbFileSpace/pdbcrystallographiccard.hpp"
#include "../../../includes/utils.hpp"
#include "../../../includes/common.hpp"

using PdbFileSpace::PdbCrystallographicCard;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbCrystallographicCard::PdbCrystallographicCard() {}

PdbCrystallographicCard::PdbCrystallographicCard(std::stringstream& stream_block)
{
    std::string line;
    getline(stream_block, line);
    std::string temp = line;
    while (!gmml::Trim(temp).empty())
    {
        record_name_ = line.substr(0,6);
        gmml::Trim(record_name_);
        if(line.substr(6, 9) != "         ")
            a_ = gmml::ConvertString<double>(line.substr(6,9));
        else
            a_ = gmml::dNotSet;
        if(line.substr(15, 9) != "         ")
            b_ = gmml::ConvertString<double>(line.substr(15,9));
        else
            b_ = gmml::dNotSet;
        if(line.substr(24, 9) != "         ")
            c_ = gmml::ConvertString<double>(line.substr(24,9));
        else
            c_ = gmml::dNotSet;
        if(line.substr(33, 7) != "       ")
            alpha_ = gmml::ConvertString<double>(line.substr(33,7));
        else
            alpha_ = gmml::dNotSet;
        if(line.substr(40, 7) != "       ")
            beta_ = gmml::ConvertString<double>(line.substr(40,7));
        else
            beta_ = gmml::dNotSet;
        if(line.substr(47, 7) != "       ")
            gamma_ = gmml::ConvertString<double>(line.substr(47,7));
        else
            gamma_ = gmml::dNotSet;
        space_group_ = line.substr(55,11);
        gmml::Trim(space_group_);
        if(line.substr(66, 4) != "    ")
            z_value_ = gmml::ConvertString<int>(line.substr(66,4));
        else
            z_value_ = gmml::iNotSet;

        getline(stream_block, line);
        temp = line;
    }
}


//////////////////////////////////////////////////////////
//                       ACCESSOR                       //
//////////////////////////////////////////////////////////

std::string PdbCrystallographicCard::GetRecordName(){
    return record_name_;
}

double PdbCrystallographicCard::GetA(){
    return a_;
}

double PdbCrystallographicCard::GetB(){
    return b_;
}

double PdbCrystallographicCard::GetC(){
    return c_;
}

double PdbCrystallographicCard::GetAlpha(){
    return alpha_;
}

double PdbCrystallographicCard::GetBeta(){
    return beta_;
}

double PdbCrystallographicCard::GetGamma(){
    return gamma_;
}

std::string PdbCrystallographicCard::GetSpaceGroup(){
    return space_group_;
}

int PdbCrystallographicCard::GetZValue(){
    return z_value_;
}

//////////////////////////////////////////////////////////
//                       MUTATOR                        //
//////////////////////////////////////////////////////////

void PdbCrystallographicCard::SetRecordName(const std::string record_name){
    record_name_ = record_name;
}

void PdbCrystallographicCard::SetA(double a){
    a_ = a;
}

void PdbCrystallographicCard::SetB(double b){
    b_ = b;
}

void PdbCrystallographicCard::SetC(double c){
    c_ = c;
}

void PdbCrystallographicCard::SetAlpha(double alpha){
    alpha_ = alpha;
}

void PdbCrystallographicCard::SetBeta(double beta){
    beta_ = beta;
}

void PdbCrystallographicCard::SetGamma(double gamma){
    gamma_ = gamma;
}

void PdbCrystallographicCard::SetSpaceGroup(const std::string space_group){
    space_group_ = space_group;
}

void PdbCrystallographicCard::SetZValue(int z_value){
    z_value_ = z_value;
}

//////////////////////////////////////////////////////////
//                       DISPLAY FUNCTION               //
//////////////////////////////////////////////////////////
void PdbCrystallographicCard::Print(std::ostream &out)
{
    out << "Record Name: " << record_name_
        << ", A: ";
    if(a_ != gmml::dNotSet)
        out << a_;
    else
        out << " ";
    out << ", B: ";
    if(b_ != gmml::dNotSet)
        out << b_;
    else
        out << " ";
    out << ", C: ";
    if(c_ != gmml::dNotSet)
        out << c_;
    else
        out << " ";
    out << ", Alpha: ";
    if(alpha_ != gmml::dNotSet)
        out << alpha_;
    else
        out << " ";
    out << ", Beta: ";
    if(beta_ != gmml::dNotSet)
        out << beta_;
    else
        out << " ";
    out << ", Gamma: ";
    if(gamma_ != gmml::dNotSet)
        out << gamma_;
    else
        out << " ";
    out << "Space Group: " << space_group_
        << ", Z Value: ";
    if(z_value_ != gmml::iNotSet)
        out << z_value_;
    else
        out << " ";
    out << std::endl << std::endl;
}
