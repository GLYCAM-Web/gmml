#include "../../../includes/FileSet/PdbFileSpace/pdbmodeltypecard.hpp"
#include "../../../includes/utils.hpp"

using namespace std;
using namespace PdbFileSpace;
using namespace gmml;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbModelTypeCard::PdbModelTypeCard() : record_name_("MDLTYP") {}

PdbModelTypeCard::PdbModelTypeCard(const string &record_name, const vector<string> &comments) : record_name_(record_name), comments_(comments) {}

PdbModelTypeCard::PdbModelTypeCard(stringstream& stream_block)
{
    string line;
    bool is_record_name_set = false;
    stringstream ss;
    string temp;
    getline(stream_block, line);
    line = Trim(line);
    while (!Trim(line).empty())
    {
        if(!is_record_name_set){
            record_name_ = line.substr(0,6);
            is_record_name_set=true;
        }
        ss << line.substr(10,70);

        getline(stream_block, line);
    }
    temp = ss.str();
    comments_ = Split(temp, ",");
}

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
