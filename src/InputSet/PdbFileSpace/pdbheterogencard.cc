// Author: Alireza Khatamian

#include "../../../includes/InputSet/PdbFileSpace/pdbheterogencard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbheterogen.hpp"
#include "../../../includes/utils.hpp"

using namespace std;
using namespace PdbFileSpace;
using namespace gmml;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbHeterogenCard::PdbHeterogenCard() : record_name_("HET") {}

PdbHeterogenCard::PdbHeterogenCard(const string &record_name) : record_name_(record_name) {}

PdbHeterogenCard::PdbHeterogenCard(stringstream &stream_block)
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
        PdbHeterogen* heterogen = new PdbHeterogen(ss);
        stringstream key;
        key << heterogen->GetChainId() << "_" << heterogen->GetSequenceNumber() << "_" << heterogen->GetInsertionCode();
        heterogens_[key.str()] = heterogen;
        getline(stream_block, line);
        temp = line;
    }
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
std::string PdbHeterogenCard::GetRecordName()
{
    return record_name_;
}

PdbHeterogenCard::HeterogenMap PdbHeterogenCard::GetHeterogens()
{
    return heterogens_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbHeterogenCard::SetRecordName(const string record_name)
{
    record_name_ = record_name;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbHeterogenCard::Print(ostream &out)
{
    out << "Record Name: " << record_name_ << endl <<
           "=============== Heterogen ==============" << endl;
    for(PdbHeterogenCard::HeterogenMap::iterator it = heterogens_.begin(); it != heterogens_.end(); it++)
    {
        out << "Heterogen ID: " << (it)->first << endl;
        (it)->second->Print();
        out << endl;
    }
}
