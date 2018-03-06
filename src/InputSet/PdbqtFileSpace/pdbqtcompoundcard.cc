#include "../../../includes/InputSet/PdbqtFileSpace/pdbqtcompoundcard.hpp"
#include "../../../includes/utils.hpp"
#include "../../../includes/common.hpp"

using PdbqtFileSpace::PdbqtCompoundCard;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////

PdbqtCompoundCard::PdbqtCompoundCard() : record_name_("COMPND"){}

PdbqtCompoundCard::PdbqtCompoundCard(std::string line)
{
    std::vector<std::string> tokens = gmml::Split(line, " ");
    record_name_ = tokens.at(0);
    value_ = tokens.at(1);
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
std::string PdbqtCompoundCard::GetRecordName()
{
    return record_name_;
}

std::string PdbqtCompoundCard::GetValue()
{
    return value_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbqtCompoundCard::SetRecordName(const std::string record_name)
{
    record_name_ = record_name;
}

void PdbqtCompoundCard::SetValue(const std::string value)
{
    value_ = value;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbqtCompoundCard::Print(std::ostream &out)
{
    out << record_name_ << "   " << value_ << std::endl;
}
