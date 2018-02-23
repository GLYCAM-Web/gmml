#include "../../../includes/InputSet/PdbFileSpace/pdbsitesection.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbsitecard.hpp"
#include "../../../includes/utils.hpp"

using namespace std;
using namespace gmml;
using namespace PdbFileSpace;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbSiteSection::PdbSiteSection() {}

PdbSiteSection::PdbSiteSection(stringstream &stream_block)
{
    string line;
    bool is_record_name_set = false;
    getline(stream_block, line);
    string temp = line;
    while (!Trim(temp).empty())
    {
        if(!is_record_name_set){
            record_name_ = line.substr(0,6);
            Trim(record_name_);
            is_record_name_set=true;
        }
        stringstream site_block;
        site_block << line << endl;
        string site_id = line.substr(11,3);

        getline(stream_block, line);
        temp = line;

        while (!Trim(temp).empty() && line.substr(11,3) == site_id){
            site_block << line << endl;
            getline(stream_block, line);
            temp = line;
        }
        PdbSiteCard* site = new PdbSiteCard(site_block);
        site_id = Trim(site_id);
        residue_site_cards_[site->GetSiteId()] = site;
    }
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////

string PdbSiteSection::GetRecordName(){
    return record_name_;
}

PdbSiteSection::PdbSiteCardMap PdbSiteSection::GetResidueSiteCards(){
    return residue_site_cards_;
}

//////////////////////////////////////////////////////////
//                       MUTATOR                        //
//////////////////////////////////////////////////////////

void PdbSiteSection::SetRecordName(const string record_name){
    record_name_ = record_name;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbSiteSection::Print(ostream &out)
{
    out << "Record Name: " << record_name_ << endl <<
           "==================== Residue Site ==============" << endl;
    for(PdbSiteSection::PdbSiteCardMap::iterator it = residue_site_cards_.begin(); it != residue_site_cards_.end(); it++)
    {
        out << "Site ID: " << (it)->first << endl;
        (it)->second->Print();
        out << endl;
    }
}
