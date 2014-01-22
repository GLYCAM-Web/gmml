#include "../../../includes/FileSet/PdbFileSpace/pdbsite.hpp"
#include "../../../includes/FileSet/PdbFileSpace/pdbsiteresidue.hpp"
#include "../../../includes/utils.hpp"

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
    line = Trim(line);
    int residue_counter = 0;
    while (!Trim(line).empty())
    {
        if(!is_site_id_set){
            site_id_ = line.substr(0,6);
            is_site_id_set=true;
        }
        number_of_residues_ = ConvertString<int>(line.substr(15, 2));

        while(residue_counter < number_of_residues_ || residue_counter % 4 != 0)
        {
            PdbSiteResidue *residue = new PdbSiteResidue(line.substr(18, 10));
            residues_.push_back(residue);
            residue_counter++;
        }

        getline(site_block, line);
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

void PdbSite::SetResidues(const SiteResidueVector residues){
    residues_ = residues;
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




