#include <iostream>

#include "../../../includes/InputSet/PdbFileSpace/pdbsplitsection.hpp"
#include "../../../includes/utils.hpp"

using namespace std;
using namespace PdbFileSpace;
using namespace gmml;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbSplitSection::PdbSplitSection() : record_name_("SPLIT"), split_(""){}

PdbSplitSection::PdbSplitSection(const string &record_name, const string &split)
{
    record_name_ = record_name;
    split_ = split;
}

PdbSplitSection::PdbSplitSection(stringstream& stream_block)
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
    split_ = ss.str();
    Trim(split_);
}

//////////////////////////////////////////////////////////
//                       ACCESSOR                       //
//////////////////////////////////////////////////////////
string PdbSplitSection::GetRecordName()
{
    return record_name_;
}

string PdbSplitSection::GetSplit()
{
    return split_;
}

//////////////////////////////////////////////////////////
//                       MUTATOR                        //
//////////////////////////////////////////////////////////
void PdbSplitSection::SetRecordName(const string record_name)
{
    record_name_ = record_name;
}

void PdbSplitSection::SetSplit(const string split)
{
    split_ = split;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbSplitSection::Print(ostream &out)
{
    out << "Record Name: " << record_name_ << ", Split: " << split_ << endl << endl;
}
