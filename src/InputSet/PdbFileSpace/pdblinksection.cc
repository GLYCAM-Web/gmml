#include "../../../includes/InputSet/PdbFileSpace/pdblinkcard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdblinksection.hpp"
#include "../../../includes/utils.hpp"

using PdbFileSpace::PdbLinkSection;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbLinkSection::PdbLinkSection() {}

PdbLinkSection::PdbLinkSection(std::stringstream &stream_block)
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

        PdbLinkCard* link_card = new PdbLinkCard(line);
        AddResidueLinkCard(link_card);
        getline(stream_block, line);
        temp = line;
    }
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////

std::string PdbLinkSection::GetRecordName(){
    return record_name_;
}

PdbLinkSection::LinkCardVector PdbLinkSection::GetResidueLinkCards(){
    return residue_link_cards_;
}

//////////////////////////////////////////////////////////
//                       MUTATOR                        //
//////////////////////////////////////////////////////////

void PdbLinkSection::SetRecordName(std::string record_name){
    record_name_ = record_name;
}

void PdbLinkSection::SetResidueLinkCards(LinkCardVector residue_link_cards){
    residue_link_cards_.clear();
    for(LinkCardVector::iterator it = residue_link_cards.begin(); it != residue_link_cards.end(); it++)
    {
        residue_link_cards_.push_back(*it);
    }
}

void PdbLinkSection::AddResidueLinkCard(PdbLinkCard *residue_link_card)
{
    residue_link_cards_.push_back(residue_link_card);
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbLinkSection::Print(std::ostream &out)
{
    out << "Record Name: " << record_name_ << std::endl <<
           "================== Residue Links ================" << std::endl;
    for(PdbLinkSection::LinkCardVector::iterator it = residue_link_cards_.begin(); it != residue_link_cards_.end(); it++)
    {
        (*it)->Print(out);
        out << std::endl;
    }
    out << std::endl;
}
