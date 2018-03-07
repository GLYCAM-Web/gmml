// Author: Alireza Khatamian

#include "../../../includes/InputSet/PdbFileSpace/pdbheterogensynonymcard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbheterogensynonymsection.hpp"
#include "../../../includes/utils.hpp"

using PdbFileSpace::PdbHeterogenSynonymSection;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbHeterogenSynonymSection::PdbHeterogenSynonymSection() : record_name_("HETSYN") {}
PdbHeterogenSynonymSection::PdbHeterogenSynonymSection(const std::string& record_name) : record_name_(record_name) {}
PdbHeterogenSynonymSection::PdbHeterogenSynonymSection(std::stringstream &stream_block)
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
        std::stringstream heterogen_synonym_block;
        heterogen_synonym_block << line << std::endl;
        std::string heterogen_identifier = line.substr(11,3);

        getline(stream_block, line);
        temp = line;

        while (!gmml::Trim(temp).empty() && line.substr(11,3) == heterogen_identifier){
            heterogen_synonym_block << line << std::endl;
            getline(stream_block, line);
            temp = line;
        }
        PdbHeterogenSynonymCard* heterogen_synonym = new PdbHeterogenSynonymCard(heterogen_synonym_block);
        heterogen_identifier = gmml::Trim(heterogen_identifier);
        heterogens_synonym_cards_[heterogen_identifier] = heterogen_synonym;
    }
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
std::string PdbHeterogenSynonymSection::GetRecordName()
{
    return record_name_;
}

PdbHeterogenSynonymSection::HeterogenSynonymCardMap PdbHeterogenSynonymSection::GetHeterogensSynonymCards()
{
    return heterogens_synonym_cards_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbHeterogenSynonymSection::SetRecordName(const std::string record_name)
{
    record_name_ = record_name;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbHeterogenSynonymSection::Print(std::ostream &out)
{
    out << "Record Name: " << record_name_ << std::endl <<
           "============ Heterogen Synonyms ===========" << std::endl;
    for(PdbHeterogenSynonymSection::HeterogenSynonymCardMap::iterator it = heterogens_synonym_cards_.begin(); it != heterogens_synonym_cards_.end(); it++)
    {
        out << "Heterogen ID: " << (it)->first << std::endl;
        (it)->second->Print();
        out << std::endl;
    }
}
