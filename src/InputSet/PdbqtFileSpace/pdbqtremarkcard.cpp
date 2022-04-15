#include "../../../includes/InputSet/PdbqtFileSpace/pdbqtremarkcard.hpp"
#include "../../../includes/utils.hpp"
#include "../../../includes/common.hpp"

using PdbqtFileSpace::PdbqtRemarkCard;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////

PdbqtRemarkCard::PdbqtRemarkCard() : record_name_("REMARK"){}

PdbqtRemarkCard::PdbqtRemarkCard(std::string line)
{
    record_name_ = line.substr(0,6);
    value_ = line.substr(6);
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
std::string PdbqtRemarkCard::GetRecordName()
{
    return record_name_;
}

std::string PdbqtRemarkCard::GetValue()
{
    return value_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbqtRemarkCard::SetRecordName(const std::string record_name)
{
    record_name_ = record_name;
}

void PdbqtRemarkCard::SetValue(const std::string value)
{
    value_ = value;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbqtRemarkCard::Print(std::ostream &out)
{
    out << "REMARK: " << value_ << std::endl;
}
