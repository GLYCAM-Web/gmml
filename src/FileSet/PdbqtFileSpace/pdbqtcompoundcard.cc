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
    vector<string> tokens = Split(line, " ");
    record_name_ = tokens.at(0);
    value_ = tokens.at(1);
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
    out << record_name_ << "   " << value_ << endl;
}


