// Author: Alireza Khatamian

#include "../../../includes/FileSet/PdbFileSpace/pdbheterogenname.hpp"

using namespace std;
using namespace PdbFileSpace;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbHeterogenName::PdbHeterogenName() : heterogen_identifier_(""), heterogen_name_("") {}
PdbHeterogenName::PdbHeterogenName(const string &heterogen_identifier, const string &heterogen_name) :
    heterogen_identifier_(heterogen_identifier), heterogen_name_(heterogen_name) {}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
string PdbHeterogenName::GetHeterogenIdentifier()
{
    return heterogen_identifier_;
}

string PdbHeterogenName::GetHeterogenName()
{
    return heterogen_name_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbHeterogenName::SetHeterogenIdentifier(const string heterogen_identifier)
{
    heterogen_identifier_ = heterogen_identifier;
}

void PdbHeterogenName::SetHeterogenName(const string heterogen_name)
{
    heterogen_name_ = heterogen_name;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////

