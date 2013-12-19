#include "../../../includes/FileSet/PdbFileSpace/pdbsite.hpp"
#include "../../../includes/FileSet/PdbFileSpace/pdbsiteresidue.hpp"

using namespace std;
using namespace PdbFileSpace;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbSite::PdbSite() {}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////

int PdbSite::GetSiteId(){
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

void PdbSite::SetSiteId(int site_id){
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




