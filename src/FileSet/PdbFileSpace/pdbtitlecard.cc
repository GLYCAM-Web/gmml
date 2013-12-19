#include "../../../includes/FileSet/PdbFileSpace/pdbtitlecard.hpp"

using namespace std;
using namespace PdbFileSpace;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbTitleCard::PdbTitleCard() : record_name_("TITLE"), title_(""){}

PdbTitleCard::PdbTitleCard(const string &record_name, const string &title)
{
    record_name_ = record_name;
    title_ = title;
}


//////////////////////////////////////////////////////////
//                       ACCESSOR                       //
//////////////////////////////////////////////////////////
string PdbTitleCard::GetRecordName()
{
    return record_name_;
}

string PdbTitleCard::GetTitle()
{
    return title_;
}

//////////////////////////////////////////////////////////
//                       MUTATOR                        //
//////////////////////////////////////////////////////////
void PdbTitleCard::SetRecordName(const string record_name)
{
    record_name_ = record_name;
}

void PdbTitleCard::SetTitle(const string title)
{
    title_ = title;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////