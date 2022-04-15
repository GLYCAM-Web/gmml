#include "../../../includes/InputSet/PdbqtFileSpace/pdbqttorsionaldofcard.hpp"
#include "../../../includes/utils.hpp"
#include "../../../includes/common.hpp"

using PdbqtFileSpace::PdbqtTorsionalDoFCard;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbqtTorsionalDoFCard::PdbqtTorsionalDoFCard() : record_name_("TORSDOF"){}

PdbqtTorsionalDoFCard::PdbqtTorsionalDoFCard(std::string line)
{
    record_name_ = line.substr(0, 7);
    std::string temp = line.substr(7);
    temp = gmml::Trim(temp);
    number_of_tosional_dof_ = gmml::ConvertString<int>(temp);
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
std::string PdbqtTorsionalDoFCard::GetRecordName()
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
void PdbqtTorsionalDoFCard::SetRecordName(const std::string record_name)
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
void PdbqtTorsionalDoFCard::Print(std::ostream &out)
{
    out << "Torsional DOF: ";
    if(number_of_tosional_dof_ != gmml::iNotSet)
        out << number_of_tosional_dof_;
    else
        out << "";
    out << std::endl;
}
