// Author: Alireza Khatamian

#include "../../../includes/InputSet/PdbFileSpace/pdbresiduemodificationcard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbresiduemodificationsection.hpp"
#include "../../../includes/utils.hpp"

using PdbFileSpace::PdbResidueModificationSection;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbResidueModificationSection::PdbResidueModificationSection() : record_name_("MODRES") {}

PdbResidueModificationSection::PdbResidueModificationSection(const std::string &record_name) : record_name_(record_name) {}

PdbResidueModificationSection::PdbResidueModificationSection(std::stringstream &stream_block)
{
    std::string line;
    std::stringstream ss;
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

        ss << line << std::endl;
        PdbResidueModificationCard* residue_modification_cards = new PdbResidueModificationCard(ss);
        residue_modification_cards_[residue_modification_cards->GetIdCode()] = residue_modification_cards;
        getline(stream_block, line);
        temp = line;
    }
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
std::string PdbResidueModificationSection::GetRecordName()
{
    return record_name_;
}

PdbResidueModificationSection::ResidueModificationCardMap PdbResidueModificationSection::GetResidueModificationCards()
{
    return residue_modification_cards_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbResidueModificationSection::SetRecordName(const std::string record_name)
{
    record_name_ = record_name;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbResidueModificationSection::Print(std::ostream &out)
{
    out << "Record Name: " << record_name_ << std::endl <<
           "============ Residue Modification ===========" << std::endl;
    for(PdbResidueModificationSection::ResidueModificationCardMap::iterator it = residue_modification_cards_.begin(); it != residue_modification_cards_.end(); it++)
    {
        out << "ID Code: " << (it)->first << std::endl;
        (it)->second->Print();
        out << std::endl;
    }
}
