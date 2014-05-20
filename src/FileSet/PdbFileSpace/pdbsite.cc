#include "../../../includes/FileSet/PdbFileSpace/pdbsite.hpp"
#include "../../../includes/FileSet/PdbFileSpace/pdbsiteresidue.hpp"
#include "../../../includes/utils.hpp"
#include "../../../includes/common.hpp"

using namespace std;
using namespace gmml;
using namespace PdbFileSpace;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbSite::PdbSite() {}
PdbSite::PdbSite(stringstream &site_block)
{
    string line;
    bool is_site_id_set = false;
    getline(site_block, line);
    string temp = line;
    int residue_counter = 0;
    while (!Trim(temp).empty())
    {
        if(!is_site_id_set){
            site_id_ = line.substr(11,3);
            Trim(site_id_);
            is_site_id_set=true;
        }
        if(line.substr(15, 2) == "  ")
            number_of_residues_ = iNotSet;
        else
            number_of_residues_ = ConvertString<int>(line.substr(15, 2));

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

string PdbSite::GetSiteId(){
    return site_id_;
}

PdbSite::SiteResidueVector PdbSite::GetResidues(){
    return residues_;
}

int PdbSite::GetNumberOfResidues(){
    return number_of_residues_;
}

//////////////////////////////////////////////////////////
//                       MUTATOR                        //
//////////////////////////////////////////////////////////

void PdbSite::SetSiteId(const string site_id){
    site_id_ = site_id;
}

void PdbSite::SetResidues(SiteResidueVector residues){
    residues_.clear();
    for(SiteResidueVector::iterator it = residues.begin(); it != residues.end(); it++)
    {
        residues_.push_back(*it);
    }
}

void PdbSite::AddResidue(PdbSiteResidue *residue)
{
    residues_.push_back(residue);
}

void PdbSite::SetNumberOfResidues(int number_of_residues){
    number_of_residues_ = number_of_residues;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbSite::Print(ostream &out)
{
    out << "Site ID: " << site_id_
        << ", Number of Residues: ";
    if(number_of_residues_ != iNotSet)
        out << number_of_residues_;
    else
        out << " ";
    out << endl << "---------------- Residues ----------------" << endl;
    for(PdbSite::SiteResidueVector::iterator it = residues_.begin(); it != residues_.end(); it++)
    {
        (*it)->Print(out);
        out << endl;
    }
    out << endl;
}
