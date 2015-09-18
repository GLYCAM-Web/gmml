#include "../../../includes/InputSet/CondensedSequenceSpace/condensedsequenceamberprepresidue.hpp"
#include "../../../includes/InputSet/CondensedSequenceSpace/condensedsequenceprocessingexception.hpp"
#include "../../../includes/common.hpp"
#include "../../../includes/utils.hpp"

using namespace std;
using namespace gmml;
using namespace CondensedSequenceSpace;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
CondensedSequenceAmberPrepResidue::CondensedSequenceAmberPrepResidue(){}

CondensedSequenceAmberPrepResidue::CondensedSequenceAmberPrepResidue(string name) : name_(name), has_parent_(false){}

CondensedSequenceAmberPrepResidue::CondensedSequenceAmberPrepResidue(string name, string anomeric_carbon, string parent_oxygen) :
    name_(name), anomeric_carbon_(anomeric_carbon), parent_oxygen_(parent_oxygen), has_parent_(true) {}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
string CondensedSequenceAmberPrepResidue::GetName()
{
    return name_;
}

string CondensedSequenceAmberPrepResidue::GetAnomericCarbon()
{
    return anomeric_carbon_;
}

string CondensedSequenceAmberPrepResidue::GetParentOxygen()
{
    return parent_oxygen_;
}

bool CondensedSequenceAmberPrepResidue::GetHasParent()
{
    return has_parent_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void CondensedSequenceAmberPrepResidue::SetName(string name)
{
    name_ = name;
}

void CondensedSequenceAmberPrepResidue::SetAnomericCarbon(string anomeric_carbon)
{
    anomeric_carbon_ = anomeric_carbon;
}

void CondensedSequenceAmberPrepResidue::SetParentOxygen(string parent_oxygen)
{
    parent_oxygen_ = parent_oxygen;
}

void CondensedSequenceAmberPrepResidue::SetHasParent(bool has_parent)
{
    has_parent_ = has_parent;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void CondensedSequenceAmberPrepResidue::Print(ostream &out)
{}






