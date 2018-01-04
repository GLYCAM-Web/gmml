#include "../../../includes/InputSet/CondensedSequenceSpace/condensedsequenceglycam06residue.hpp"
#include "../../../includes/InputSet/CondensedSequenceSpace/condensedsequenceprocessingexception.hpp"
#include "../../../includes/common.hpp"
#include "../../../includes/utils.hpp"

using namespace std;
using namespace gmml;
using namespace CondensedSequenceSpace;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
CondensedSequenceGlycam06Residue::CondensedSequenceGlycam06Residue(){}

CondensedSequenceGlycam06Residue::CondensedSequenceGlycam06Residue(string name) : name_(name), parent_id_(-1), is_derivative_(false) {} //has_parent_(false){}

CondensedSequenceGlycam06Residue::CondensedSequenceGlycam06Residue(string name, string anomeric_carbon, string parent_oxygen, bool is_derivative) :
    name_(name), anomeric_carbon_(anomeric_carbon), parent_oxygen_(parent_oxygen), parent_id_(0), is_derivative_(is_derivative) {}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
string CondensedSequenceGlycam06Residue::GetName()
{
    return name_;
}

string CondensedSequenceGlycam06Residue::GetAnomericCarbon()
{
    return anomeric_carbon_;
}

string CondensedSequenceGlycam06Residue::GetParentOxygen()
{
    return parent_oxygen_;
}

bool CondensedSequenceGlycam06Residue::GetIsDerivative()
{
    return is_derivative_;
}

int CondensedSequenceGlycam06Residue::GetParentId()
{
    return parent_id_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void CondensedSequenceGlycam06Residue::SetName(string name)
{
    name_ = name;
}

void CondensedSequenceGlycam06Residue::SetAnomericCarbon(string anomeric_carbon)
{
    anomeric_carbon_ = anomeric_carbon;
}

void CondensedSequenceGlycam06Residue::SetParentOxygen(string parent_oxygen)
{
    parent_oxygen_ = parent_oxygen;
}

void CondensedSequenceGlycam06Residue::SetIsDerivative(bool is_derivative)
{
    is_derivative_ = is_derivative;
}

void CondensedSequenceGlycam06Residue::SetParentId(int parent_id)
{
    parent_id_ = parent_id;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void CondensedSequenceGlycam06Residue::Print(ostream &out)
{}






