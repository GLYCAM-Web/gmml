#include <iostream>

#include "../../../includes/InputSet/PdbFileSpace/pdbauthorsection.hpp"
#include "../../../includes/utils.hpp"

using namespace std;
using namespace PdbFileSpace;
using namespace gmml;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbAuthorSection::PdbAuthorSection() : record_name_("AUTHOR"), author_(""){}

PdbAuthorSection::PdbAuthorSection(const string &record_name, const string &author)
{
    record_name_ = record_name;
    author_ = author;
}

PdbAuthorSection::PdbAuthorSection(stringstream& stream_block)
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
    author_ = ss.str();
    Trim(author_);
}

//////////////////////////////////////////////////////////
//                       ACCESSOR                       //
//////////////////////////////////////////////////////////
string PdbAuthorSection::GetRecordName()
{
    return record_name_;
}

string PdbAuthorSection::GetAuthor()
{
    return author_;
}

//////////////////////////////////////////////////////////
//                       MUTATOR                        //
//////////////////////////////////////////////////////////
void PdbAuthorSection::SetRecordName(const string record_name)
{
    record_name_ = record_name;
}

void PdbAuthorSection::SetAuthor(const string author)
{
    author_ = author;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbAuthorSection::Print(ostream &out)
{
    out << "Record Name: " << record_name_ << ", Author(s): " << author_ << endl << endl;
}
