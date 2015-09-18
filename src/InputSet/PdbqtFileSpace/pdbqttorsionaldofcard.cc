#include "../../../includes/InputSet/PdbqtFileSpace/pdbqttorsionaldofcard.hpp"
#include "../../../includes/utils.hpp"
#include "../../../includes/common.hpp"

using namespace std;
using namespace gmml;
using namespace PdbqtFileSpace;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbqtTorsionalDoFCard::PdbqtTorsionalDoFCard() : record_name_("TORSDOF"){}

PdbqtTorsionalDoFCard::PdbqtTorsionalDoFCard(string line)
{
    record_name_ = line.substr(0, 7);
    string temp = line.substr(7);
    temp = Trim(temp);
    number_of_tosional_dof_ = ConvertString<int>(temp);
}

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
    out << "Torsional DOF: ";
    if(number_of_tosional_dof_ != iNotSet)
        out << number_of_tosional_dof_;
    else
        out << "";
    out << endl;
}


