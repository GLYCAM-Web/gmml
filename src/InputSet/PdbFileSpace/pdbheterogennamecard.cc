// Author: Alireza Khatamian

#include "../../../includes/InputSet/PdbFileSpace/pdbheterogennamecard.hpp"
#include "../../../includes/utils.hpp"

using namespace std;
using namespace PdbFileSpace;
using namespace gmml;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbHeterogenNameCard::PdbHeterogenNameCard() : heterogen_identifier_(""), heterogen_name_("") {}
PdbHeterogenNameCard::PdbHeterogenNameCard(const string &heterogen_identifier, const string &heterogen_name) :
    heterogen_identifier_(heterogen_identifier), heterogen_name_(heterogen_name) {}

PdbHeterogenNameCard::PdbHeterogenNameCard(stringstream& stream_block)
{
    string line;
    bool is_heterogen_identifier_set = false;
    stringstream ss;
    getline(stream_block, line);
    string temp = line;
    while (!Trim(temp).empty())
    {
        if(!is_heterogen_identifier_set){
            heterogen_identifier_ = line.substr(11,3);
            Trim(heterogen_identifier_);
            is_heterogen_identifier_set = true;
        }

        ss << line.substr(15,55) << " ";

        getline(stream_block, line);
        temp = line;
    }
    heterogen_name_ = ss.str();
    heterogen_name_ = Trim(heterogen_name_);
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
string PdbHeterogenNameCard::GetHeterogenIdentifier()
{
    return heterogen_identifier_;
}

string PdbHeterogenNameCard::GetHeterogenName()
{
    return heterogen_name_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbHeterogenNameCard::SetHeterogenIdentifier(const string heterogen_identifier)
{
    heterogen_identifier_ = heterogen_identifier;
}

void PdbHeterogenNameCard::SetHeterogenName(const string heterogen_name)
{
    heterogen_name_ = heterogen_name;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbHeterogenNameCard::Print(ostream &out)
{
    out << "Heterogen ID: " << heterogen_identifier_ << ", Heterogen Name: " << heterogen_name_ << endl;
}
