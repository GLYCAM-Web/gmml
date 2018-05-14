#include "../../../includes/InputSet/PdbFileSpace/pdbsitecard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbsiteresidue.hpp"
#include "../../../includes/utils.hpp"
#include "../../../includes/common.hpp"

using PdbFileSpace::PdbSiteCard;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbSiteCard::PdbSiteCard() {}
PdbSiteCard::PdbSiteCard(std::stringstream &site_block)
{
    std::string line;
    bool is_site_id_set = false;
    getline(site_block, line);
    std::string temp = line;
    int residue_counter = 0;
    while (!gmml::Trim(temp).empty())
    {
        if(!is_site_id_set){
            site_id_ = line.substr(11,3);
            gmml::Trim(site_id_);
            is_site_id_set=true;
        }
        if(line.substr(15, 2) == "  ")
            number_of_residues_ = gmml::iNotSet;
        else
            number_of_residues_ = gmml::ConvertString<int>(line.substr(15, 2));

        int residue_counter_per_line = 0;
        while(residue_counter < number_of_residues_ && residue_counter_per_line < 4)
        {
            int index = 18 + (residue_counter%4)*11;
            PdbSiteResidue *residue = new PdbSiteResidue(line.substr(index, 10));
            residues_.push_back(residue);
            residue_counter++;
            residue_counter_per_line++;
        }

        getline(site_block, line);
        temp = line;
    }
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////

std::string PdbSiteCard::GetSiteId(){
    return site_id_;
}

PdbSiteCard::SiteResidueVector PdbSiteCard::GetResidues(){
    return residues_;
}

int PdbSiteCard::GetNumberOfResidues(){
    return number_of_residues_;
}

//////////////////////////////////////////////////////////
//                       MUTATOR                        //
//////////////////////////////////////////////////////////

void PdbSiteCard::SetSiteId(const std::string site_id){
    site_id_ = site_id;
}

void PdbSiteCard::SetResidues(SiteResidueVector residues){
    residues_.clear();
    for(SiteResidueVector::iterator it = residues.begin(); it != residues.end(); it++)
    {
        residues_.push_back(*it);
    }
}

void PdbSiteCard::AddResidue(PdbSiteResidue *residue)
{
    residues_.push_back(residue);
}

void PdbSiteCard::SetNumberOfResidues(int number_of_residues){
    number_of_residues_ = number_of_residues;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbSiteCard::Print(std::ostream &out)
{
    out << "Site ID: " << site_id_
        << ", Number of Residues: ";
    if(number_of_residues_ != gmml::iNotSet)
        out << number_of_residues_;
    else
        out << " ";
    out << std::endl << "---------------- Residues ----------------" << std::endl;
    for(PdbSiteCard::SiteResidueVector::iterator it = residues_.begin(); it != residues_.end(); it++)
    {
        (*it)->Print(out);
        out << std::endl;
    }
    out << std::endl;
}
