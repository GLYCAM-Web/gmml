// Author: Alireza Khatamian

#include "../../../includes/FileSet/PdbFileSpace/pdbheterogennamecard.hpp"

using namespace std;
using namespace PdbFileSpace;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbHeterogenNameCard::PdbHeterogenNameCard() : record_name_("HETNAM") {}
PdbHeterogenNameCard::PdbHeterogenNameCard(const string &record_name) : record_name_(record_name) {}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
string PdbHeterogenNameCard::GetRecordName()
{
    return record_name_;
}

PdbHeterogenNameCard::HeterogenNameMap PdbHeterogenNameCard::GetHeterogenNames()
{
    return heterogen_names_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbHeterogenNameCard::SetRecordName(const string record_name)
{
    record_name_ = record_name;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////

