#include "../../../includes/InputSet/CondensedSequenceSpace/condensedsequenceglycam06residue.hpp"
#include "../../../includes/InputSet/CondensedSequenceSpace/condensedsequenceprocessingexception.hpp"
#include "../../../includes/common.hpp"
#include "../../../includes/utils.hpp"

using CondensedSequenceSpace::CondensedSequenceGlycam06Residue;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
CondensedSequenceGlycam06Residue::CondensedSequenceGlycam06Residue(){}

CondensedSequenceGlycam06Residue::CondensedSequenceGlycam06Residue(std::string name) : name_(name), is_derivative_(false), parent_id_(gmml::iNotSet) {} //has_parent_(false){},originally pare_id -1

CondensedSequenceGlycam06Residue::CondensedSequenceGlycam06Residue(std::string name, std::string anomeric_carbon, std::string parent_oxygen, bool is_derivative) :
    name_(name), anomeric_carbon_(anomeric_carbon), parent_oxygen_(parent_oxygen), is_derivative_(is_derivative), parent_id_(gmml::iNotSet) {} //originally parent_id = 0

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
std::string CondensedSequenceGlycam06Residue::GetName()
{
    return name_;
}

std::string CondensedSequenceGlycam06Residue::GetAnomericCarbon()
{
    return anomeric_carbon_;
}

std::string CondensedSequenceGlycam06Residue::GetParentOxygen()
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
int CondensedSequenceGlycam06Residue::GetBondId()  //Added by Yao 08/03/2018
{
    return bond_id_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void CondensedSequenceGlycam06Residue::SetName(std::string name)
{
    name_ = name;
}

void CondensedSequenceGlycam06Residue::SetAnomericCarbon(std::string anomeric_carbon)
{
    anomeric_carbon_ = anomeric_carbon;
}

void CondensedSequenceGlycam06Residue::SetParentOxygen(std::string parent_oxygen)
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
void CondensedSequenceGlycam06Residue::SetBondId(int bond_id)  //Added by Yao 08/03/2018
{
    bond_id_ = bond_id;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void CondensedSequenceGlycam06Residue::Print(std::ostream &out)
{
    out << "";
}
