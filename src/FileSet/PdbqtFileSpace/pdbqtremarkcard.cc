#include "../../../includes/FileSet/PdbqtFileSpace/pdbqtremarkcard.hpp"
#include "../../../includes/utils.hpp"
#include "../../../includes/common.hpp"

using namespace std;
using namespace gmml;
using namespace PdbqtFileSpace;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////

PdbqtRemarkCard::PdbqtRemarkCard() : record_name_("REMARK"){}

PdbqtRemarkCard::PdbqtRemarkCard(string line)
{
    record_name_ = line.substr(0,6);
    value_ = line.substr(10);
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
string PdbqtRemarkCard::GetRecordName()
{
    return record_name_;
}

string PdbqtRemarkCard::GetValue()
{
    return value_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbqtRemarkCard::SetRecordName(const string record_name)
{
    record_name_ = record_name;
}

void PdbqtRemarkCard::SetValue(const string value)
{
    value_ = value;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbqtRemarkCard::Print(ostream &out)
{
    out << "REMARK: " << value_ << endl;
}

