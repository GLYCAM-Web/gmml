// Author: Alireza Khatamian

#include "../../../includes/InputSet/PdbFileSpace/pdbheterogencard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbheterogensection.hpp"
#include "../../../includes/utils.hpp"

using PdbFileSpace::PdbHeterogenSection;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbHeterogenSection::PdbHeterogenSection() : record_name_("HET") {}

PdbHeterogenSection::PdbHeterogenSection(const std::string &record_name) : record_name_(record_name) {}

PdbHeterogenSection::PdbHeterogenSection(std::stringstream &stream_block)
{
    std::string line;
    bool is_record_name_set = false;
    getline(stream_block, line);
    std::string temp = line;
    while (!gmml::Trim(temp).empty())
    {
        if(!is_record_name_set){
            record_name_ = line.substr(0,6);
            gmml::Trim(record_name_);
            is_record_name_set=true;
        }
        std::stringstream ss;
        ss << line << std::endl;
        PdbHeterogenCard* heterogen = new PdbHeterogenCard(ss);
        std::stringstream key;
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
void PdbHeterogenSection::SetRecordName(const std::string record_name)
{
    record_name_ = record_name;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbHeterogenSection::Print(std::ostream &out)
{
    out << "Record Name: " << record_name_ << std::endl <<
           "=============== Heterogen ==============" << std::endl;
    for(PdbHeterogenSection::HeterogenCardMap::iterator it = heterogen_cards_.begin(); it != heterogen_cards_.end(); it++)
    {
        out << "Heterogen ID: " << (it)->first << std::endl;
        (it)->second->Print();
        out << std::endl;
    }
}
