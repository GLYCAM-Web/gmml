#include "../../../includes/FileSet/PdbFileSpace/pdbtitlecard.hpp"
#include "../../../includes/utils.hpp"

using namespace std;
using namespace PdbFileSpace;
using namespace gmml;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbTitleCard::PdbTitleCard() : record_name_("TITLE"), title_(""){}

PdbTitleCard::PdbTitleCard(const string &record_name, const string &title)
{
    record_name_ = record_name;
    title_ = title;
}

PdbTitleCard::PdbTitleCard(istringstream& stream_block)
{
    string line;
    getline(stream_block, line);
    line = Trim(line);
    while (!Trim(line).empty())
    {
        record_name_ = line.substr(0,6);
        classification_ = line.substr(10,40);
        deposition_date_ = line.substr(50, 9);
        identifier_code_ = line.substr(62,4);
        getline(stream_block, line);
    }
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
