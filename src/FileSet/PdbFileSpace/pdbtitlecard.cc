#include <iostream>

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

PdbTitleCard::PdbTitleCard(stringstream& stream_block)
{
    string line;
    bool is_record_name_set = false;
    stringstream ss;
    getline(stream_block, line);
    string temp = line;
    while (!Trim(temp).empty())
    {
        if(!is_record_name_set){
            record_name_ = line.substr(0,6);
            Trim(record_name_);
            is_record_name_set=true;
        }
        ss << line.substr(10,70);

        getline(stream_block, line);
        temp = line;
    }
    title_ = ss.str();
    Trim(title_);
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
void PdbTitleCard::Print(ostream &out)
{
    out << "Record Name: " << record_name_ << ", Title: " << title_ << endl << endl;
}
