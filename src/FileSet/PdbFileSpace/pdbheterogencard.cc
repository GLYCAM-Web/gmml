// Author: Alireza Khatamian

#include "../../../includes/FileSet/PdbFileSpace/pdbheterogencard.hpp"
#include "../../../includes/FileSet/PdbFileSpace/pdbheterogen.hpp"
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
    stringstream ss;
    bool is_record_name_set = false;
    getline(stream_block, line);
    line = Trim(line);
    while (!Trim(line).empty())
    {
        if(!is_record_name_set){
            record_name_ = line.substr(0,6);
            is_record_name_set=true;
        }

        ss << line;
        PdbHeterogen* heterogen = new PdbHeterogen(ss);
        heterogens_[line.substr(7,3)] = heterogen;
        getline(stream_block, line);
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

