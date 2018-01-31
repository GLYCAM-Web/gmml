// Author: Alireza Khatamian

#include "../../../includes/InputSet/PdbFileSpace/pdbheterogennamecard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbheterogennamesection.hpp"
#include "../../../includes/utils.hpp"

using namespace std;
using namespace PdbFileSpace;
using namespace gmml;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbHeterogenNameSection::PdbHeterogenNameSection() : record_name_("HETNAM") {}
PdbHeterogenNameSection::PdbHeterogenNameSection(const string &record_name) : record_name_(record_name) {}

PdbHeterogenNameSection::PdbHeterogenNameSection(stringstream &stream_block)
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
        PdbHeterogenNameCard* heterogen_name = new PdbHeterogenNameCard(heterogen_name_block);
        heterogen_id = Trim(heterogen_id);
        heterogen_name_cards_[heterogen_id] = heterogen_name;
    }
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
string PdbHeterogenNameSection::GetRecordName()
{
    return record_name_;
}

PdbHeterogenNameSection::HeterogenNameCardMap PdbHeterogenNameSection::GetHeterogenNameCards()
{
    return heterogen_name_cards_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbHeterogenNameSection::SetRecordName(const string record_name)
{
    record_name_ = record_name;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbHeterogenNameSection::Print(ostream &out)
{
    out << "Record Name: " << record_name_ << endl <<
           "========== Heterogen Names ==========" << endl;
    for(PdbHeterogenNameSection::HeterogenNameCardMap::iterator it = heterogen_name_cards_.begin(); it != heterogen_name_cards_.end(); it++)
    {
        out << "Heterogen ID: " << (it)->first << endl;
        (it)->second->Print();
        out << endl;
    }
}
