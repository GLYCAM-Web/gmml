#include <iostream>

#include "../../../includes/InputSet/PdbFileSpace/pdbkeywordssection.hpp"
#include "../../../includes/utils.hpp"

using namespace std;
using namespace PdbFileSpace;
using namespace gmml;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbKeywordsSection::PdbKeywordsSection() : record_name_("KEYWDS"), keywords_(""){}

PdbKeywordsSection::PdbKeywordsSection(const string &record_name, const string &keywords)
{
    record_name_ = record_name;
    keywords_ = keywords;
}

PdbKeywordsSection::PdbKeywordsSection(stringstream& stream_block)
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
    keywords_ = ss.str();
    Trim(keywords_);
}

//////////////////////////////////////////////////////////
//                       ACCESSOR                       //
//////////////////////////////////////////////////////////
string PdbKeywordsSection::GetRecordName()
{
    return record_name_;
}

string PdbKeywordsSection::GetKeywords()
{
    return keywords_;
}

//////////////////////////////////////////////////////////
//                       MUTATOR                        //
//////////////////////////////////////////////////////////
void PdbKeywordsSection::SetRecordName(const string record_name)
{
    record_name_ = record_name;
}

void PdbKeywordsSection::SetKeywords(const string keywords)
{
    keywords_ = keywords;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbKeywordsSection::Print(ostream &out)
{
    out << "Record Name: " << record_name_ << ", Keywords: " << keywords_ << endl << endl;
}
