#include "../../../includes/FileSet/PdbqtFileSpace/pdbqttorsionaldofcard.hpp"
#include "../../../includes/utils.hpp"
#include "../../../includes/common.hpp"

using namespace std;
using namespace gmml;
using namespace PdbqtFileSpace;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbqtTorsionalDoFCard::PdbqtTorsionalDoFCard(){}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
string PdbqtTorsionalDoFCard::GetRecordName()
{
    return record_name_;
}

int PdbqtTorsionalDoFCard::GetNumberofTorsionalDoF()
{
    return number_of_tosional_dof_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbqtTorsionalDoFCard::SetRecordName(const string record_name)
{
    record_name_ = record_name;
}

void PdbqtTorsionalDoFCard::SetNumberOfTorsionalDoF(int number_of_torsional_dof)
{
    number_of_tosional_dof_ = number_of_torsional_dof;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbqtTorsionalDoFCard::Print(ostream &out)
{
}


