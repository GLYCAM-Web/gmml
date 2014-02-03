#include "../../../includes/FileSet/PdbFileSpace/pdbcrystallographiccard.hpp"
#include "../../../includes/utils.hpp"

using namespace std;
using namespace PdbFileSpace;
using namespace gmml;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbCrystallographicCard::PdbCrystallographicCard() {}

PdbCrystallographicCard::PdbCrystallographicCard(stringstream& stream_block)
{
    string line;
    getline(stream_block, line);
    string temp = line;
    while (!Trim(temp).empty())
    {
        record_name_ = line.substr(0,6);
        Trim(record_name_);
        a_ = ConvertString<double>(line.substr(6,9));
        b_ = ConvertString<double>(line.substr(15,9));
        c_ = ConvertString<double>(line.substr(24,9));
        alpha_ = ConvertString<double>(line.substr(33,7));
        beta_ = ConvertString<double>(line.substr(40,7));
        gamma_ = ConvertString<double>(line.substr(47,7));
        space_group_ = line.substr(55,11);
        Trim(space_group_);
        z_value_ = ConvertString<int>(line.substr(66,4));

        getline(stream_block, line);
        temp = line;
    }
}


//////////////////////////////////////////////////////////
//                       ACCESSOR                       //
//////////////////////////////////////////////////////////

string PdbCrystallographicCard::GetRecordName(){
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

string PdbCrystallographicCard::GetSpaceGroup(){
    return space_group_;
}

int PdbCrystallographicCard::GetZValue(){
    return z_value_;
}

//////////////////////////////////////////////////////////
//                       MUTATOR                        //
//////////////////////////////////////////////////////////

void PdbCrystallographicCard::SetRecordName(const string record_name){
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

void PdbCrystallographicCard::SetSpaceGroup(const string space_group){
    space_group_ = space_group;
}

void PdbCrystallographicCard::SetZValue(int z_value){
    z_value_ = z_value;
}

//////////////////////////////////////////////////////////
//                       DISPLAY FUNCTION               //
//////////////////////////////////////////////////////////
void PdbCrystallographicCard::Print(ostream &out)
{
    out << "Record Name: " << record_name_ << ", A: " << a_ << ", B: " << b_ << ", C: " << c_ <<
           ", Alpha: " << alpha_ << ", Beta: " << beta_ << ", Gamma: " << gamma_ << "Space Group: " << space_group_ << ", Z Value: " << z_value_ << endl << endl;
}
