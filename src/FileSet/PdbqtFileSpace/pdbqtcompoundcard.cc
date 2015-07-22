#include "../../../includes/FileSet/PdbqtFileSpace/pdbqtcompoundcard.hpp"
#include "../../../includes/utils.hpp"
#include "../../../includes/common.hpp"

using namespace std;
using namespace gmml;
using namespace PdbqtFileSpace;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////

PdbqtCompoundCard::PdbqtCompoundCard() : record_name_("COMPND"){}

PdbqtCompoundCard::PdbqtCompoundCard(string line)
{
    record_name_ = line.substr(0, 6);
    value_ = line.substr(10);
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
string PdbqtCompoundCard::GetRecordName()
{
    return record_name_;
}

string PdbqtCompoundCard::GetValue()
{
    return value_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbqtCompoundCard::SetRecordName(const string record_name)
{
    record_name_ = record_name;
}

void PdbqtCompoundCard::SetValue(const string value)
{
    value_ = value;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbqtCompoundCard::Print(ostream &out)
{
}


