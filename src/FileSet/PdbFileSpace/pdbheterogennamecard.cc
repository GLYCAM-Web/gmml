// Author: Alireza Khatamian

#include "../../../includes/FileSet/PdbFileSpace/pdbheterogennamecard.hpp"
#include "../../../includes/FileSet/PdbFileSpace/pdbheterogenname.hpp"
#include "../../../includes/utils.hpp"

using namespace std;
using namespace PdbFileSpace;
using namespace gmml;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbHeterogenNameCard::PdbHeterogenNameCard() : record_name_("HETNAM") {}
PdbHeterogenNameCard::PdbHeterogenNameCard(const string &record_name) : record_name_(record_name) {}

PdbHeterogenNameCard::PdbHeterogenNameCard(stringstream &stream_block)
{
    string line;
    bool is_record_name_set = false;
    getline(stream_block, line);
    string temp = line;
    while (!Trim(temp).empty())
    {
        if(!is_record_name_set){
            record_name_ = line.substr(0,6);
            Trim(record_name_);
            is_record_name_set=true;
        }
        stringstream heterogen_name_block;
        heterogen_name_block << line << endl;
        string heterogen_id = line.substr(11,3);

        getline(stream_block, line);
        temp = line;
        while (!Trim(temp).empty() && line.substr(11,3) == heterogen_id){
            heterogen_name_block << line << endl;
            getline(stream_block, line);
            temp = line;
        }
        PdbHeterogenName* heterogen_name = new PdbHeterogenName(heterogen_name_block);
        heterogen_names_[heterogen_id] = heterogen_name;
    }
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
string PdbHeterogenNameCard::GetRecordName()
{
    return record_name_;
}

PdbHeterogenNameCard::HeterogenNameMap PdbHeterogenNameCard::GetHeterogenNames()
{
    return heterogen_names_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbHeterogenNameCard::SetRecordName(const string record_name)
{
    record_name_ = record_name;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbHeterogenNameCard::Print(ostream &out)
{
    out << "Record Name: " << record_name_ << endl <<
           "========== Heterogen Names ==========" << endl;
    for(PdbHeterogenNameCard::HeterogenNameMap::iterator it = heterogen_names_.begin(); it != heterogen_names_.end(); it++)
    {
        out << "Heterogen ID: " << (it)->first << endl;
        (it)->second->Print();
        out << endl;
    }
}
