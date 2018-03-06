// Author: Alireza Khatamian

#include "../../../includes/InputSet/PdbFileSpace/pdbheterogennamecard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbheterogennamesection.hpp"
#include "../../../includes/utils.hpp"

using PdbFileSpace::PdbHeterogenNameSection;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbHeterogenNameSection::PdbHeterogenNameSection() : record_name_("HETNAM") {}
PdbHeterogenNameSection::PdbHeterogenNameSection(const std::string &record_name) : record_name_(record_name) {}

PdbHeterogenNameSection::PdbHeterogenNameSection(std::stringstream &stream_block)
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
        std::stringstream heterogen_name_block;
        heterogen_name_block << line << std::endl;
        std::string heterogen_id = line.substr(11,3);

        getline(stream_block, line);
        temp = line;
        while (!gmml::Trim(temp).empty() && line.substr(11,3) == heterogen_id){
            heterogen_name_block << line << std::endl;
            getline(stream_block, line);
            temp = line;
        }
        PdbHeterogenNameCard* heterogen_name = new PdbHeterogenNameCard(heterogen_name_block);
        heterogen_id = gmml::Trim(heterogen_id);
        heterogen_name_cards_[heterogen_id] = heterogen_name;
    }
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
std::string PdbHeterogenNameSection::GetRecordName()
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
void PdbHeterogenNameSection::SetRecordName(const std::string record_name)
{
    record_name_ = record_name;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbHeterogenNameSection::Print(std::ostream &out)
{
    out << "Record Name: " << record_name_ << std::endl <<
           "========== Heterogen Names ==========" << std::endl;
    for(PdbHeterogenNameSection::HeterogenNameCardMap::iterator it = heterogen_name_cards_.begin(); it != heterogen_name_cards_.end(); it++)
    {
        out << "Heterogen ID: " << (it)->first << std::endl;
        (it)->second->Print();
        out << std::endl;
    }
}
