// Author: Alireza Khatamian

#include "../../../includes/InputSet/PdbFileSpace/pdbheterogencard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbheterogensection.hpp"
#include "../../../includes/utils.hpp"

using namespace std;
using namespace PdbFileSpace;
using namespace gmml;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbHeterogenSection::PdbHeterogenSection() : record_name_("HET") {}

PdbHeterogenSection::PdbHeterogenSection(const string &record_name) : record_name_(record_name) {}

PdbHeterogenSection::PdbHeterogenSection(stringstream &stream_block)
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
        stringstream ss;
        ss << line << endl;
        PdbHeterogenCard* heterogen = new PdbHeterogenCard(ss);
        stringstream key;
        key << heterogen->GetChainId() << "_" << heterogen->GetSequenceNumber() << "_" << heterogen->GetInsertionCode();
        heterogen_cards_[key.str()] = heterogen;
        getline(stream_block, line);
        temp = line;
    }
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
std::string PdbHeterogenSection::GetRecordName()
{
    return record_name_;
}

PdbHeterogenSection::HeterogenCardMap PdbHeterogenSection::GetHeterogenCards()
{
    return heterogen_cards_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbHeterogenSection::SetRecordName(const string record_name)
{
    record_name_ = record_name;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbHeterogenSection::Print(ostream &out)
{
    out << "Record Name: " << record_name_ << endl <<
           "=============== Heterogen ==============" << endl;
    for(PdbHeterogenSection::HeterogenCardMap::iterator it = heterogen_cards_.begin(); it != heterogen_cards_.end(); it++)
    {
        out << "Heterogen ID: " << (it)->first << endl;
        (it)->second->Print();
        out << endl;
    }
}
