#include "../../../includes/InputSet/PdbFileSpace/pdbsitesection.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbsitecard.hpp"
#include "../../../includes/utils.hpp"

using PdbFileSpace::PdbSiteSection;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbSiteSection::PdbSiteSection() {}

PdbSiteSection::PdbSiteSection(std::stringstream &stream_block)
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
        std::stringstream site_block;
        site_block << line << std::endl;
        std::string site_id = line.substr(11,3);

        getline(stream_block, line);
        temp = line;

        while (!gmml::Trim(temp).empty() && line.substr(11,3) == site_id){
            site_block << line << std::endl;
            getline(stream_block, line);
            temp = line;
        }
        PdbSiteCard* site = new PdbSiteCard(site_block);
        site_id = gmml::Trim(site_id);
        residue_site_cards_[site->GetSiteId()] = site;
    }
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////

std::string PdbSiteSection::GetRecordName(){
    return record_name_;
}

PdbSiteSection::PdbSiteCardMap PdbSiteSection::GetResidueSiteCards(){
    return residue_site_cards_;
}

//////////////////////////////////////////////////////////
//                       MUTATOR                        //
//////////////////////////////////////////////////////////

void PdbSiteSection::SetRecordName(const std::string record_name){
    record_name_ = record_name;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbSiteSection::Print(std::ostream &out)
{
    out << "Record Name: " << record_name_ << std::endl <<
           "==================== Residue Site ==============" << std::endl;
    for(PdbSiteSection::PdbSiteCardMap::iterator it = residue_site_cards_.begin(); it != residue_site_cards_.end(); it++)
    {
        out << "Site ID: " << (it)->first << std::endl;
        (it)->second->Print();
        out << std::endl;
    }
}
