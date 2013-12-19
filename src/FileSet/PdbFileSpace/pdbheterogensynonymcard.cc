// Author: Alireza Khatamian

#include "../../../includes/FileSet/PdbFileSpace/pdbheterogensynonymcard.hpp"

using namespace std;
using namespace PdbFileSpace;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbHeterogenSynonymCard::PdbHeterogenSynonymCard() : record_name_("HETSYN") {}
PdbHeterogenSynonymCard::PdbHeterogenSynonymCard(const string& record_name) : record_name_(record_name) {}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
string PdbHeterogenSynonymCard::GetRecordName()
{
    return record_name_;
}

PdbHeterogenSynonymCard::HeterogenSynonymMap PdbHeterogenSynonymCard::GetHeterogensSynonyms()
{
    return heterogens_synonyms_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbHeterogenSynonymCard::SetRecordName(const string record_name)
{
    record_name_ = record_name;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////


