#include "../../../includes/InputSet/PdbFileSpace/pdbmodeltypesection.hpp"
#include "../../../includes/utils.hpp"

using namespace std;
using namespace PdbFileSpace;
using namespace gmml;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbModelTypeSection::PdbModelTypeSection() : record_name_("MDLTYP") {}

PdbModelTypeSection::PdbModelTypeSection(const string &record_name, const vector<string> &comments) : record_name_(record_name), comments_(comments) {}

PdbModelTypeSection::PdbModelTypeSection(stringstream& stream_block)
{
    string line;
    bool is_record_name_set = false;
    stringstream ss;
    string temp_comments;
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
    temp_comments = ss.str();
    comments_ = Split(temp_comments, ",");
    for(vector<string>::iterator it = comments_.begin(); it != comments_.end(); it++)
    {
        Trim(*it);
    }
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
string PdbModelTypeSection::GetRecordName()
{
    return record_name_;
}

vector<string> PdbModelTypeSection::GetComments()
{
    return comments_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbModelTypeSection::SetRecordName(const string record_name)
{
    record_name_ = record_name;
}

void PdbModelTypeSection::SetComments(const vector<string> comments)
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
void PdbModelTypeSection::Print(ostream &out)
{
    out << "Record Name: " << record_name_ << endl << "Comments: ";
    for(vector<string>::iterator it = comments_.begin(); it != comments_.end(); it++)
    {
        out << (*it) << ", ";
    }
    out << endl << endl;
}
