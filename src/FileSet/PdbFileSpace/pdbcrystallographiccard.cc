#include "../../../includes/FileSet/PdbFileSpace/pdbcrystallographiccard.hpp"

using namespace std;
using namespace PdbFileSpace;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbCrystallographicCard::PdbCrystallographicCard() {}


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



