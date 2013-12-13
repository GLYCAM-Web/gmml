#include "../../../includes/FileSet/PdbFileSpace/pdbmodeltypecard.hpp"

using namespace std;
using namespace PdbFileSpace;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbModelTypeCard::PdbModelTypeCard() : record_name_("MDLTYP") {}

PdbModelTypeCard::PdbModelTypeCard(const string &record_name, const vector<string> &comments) : record_name_(record_name), comments_(comments) {}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
string PdbModelTypeCard::GetRecordName()
{
    return record_name_;
}

vector<string> PdbModelTypeCard::GetComments()
{
    return comments_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbModelTypeCard::SetRecordName(const string record_name)
{
    record_name_ = record_name;
}

void PdbModelTypeCard::SetComments(const vector<string> comments)
{
    comments_.clear();
    for(vector<string>::const_iterator it = comments.begin(); it != comments.end(); it++)
    {
        comments_.push_back(*it);
    }
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
