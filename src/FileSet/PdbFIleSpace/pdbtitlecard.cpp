#include "../../../includes/FileSet/PdbFileSpace/pdbtitlecard.h"

using namespace std;
using namespace PdbFileSpace;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbTitleCard::PdbTitleCard() : record_name_("TITLE"), title_(""){}

PdbTitleCard::PdbTitleCard(string record_name, string title)
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
void PdbTitleCard::SetRecordName(string record_name)
{
    record_name_ = record_name;
}

void PdbTitleCard::SetTitle(string title)
{
    title_ = title;
}

//////////////////////////////////////////////////////////
//                       DISPLAY FUNCTION               //
//////////////////////////////////////////////////////////
