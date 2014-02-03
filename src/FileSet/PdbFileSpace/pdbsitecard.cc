#include "../../../includes/FileSet/PdbFileSpace/pdbsite.hpp"
#include "../../../includes/FileSet/PdbFileSpace/pdbsitecard.hpp"
#include "../../../includes/utils.hpp"

using namespace std;
using namespace gmml;
using namespace PdbFileSpace;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbSiteCard::PdbSiteCard() {}

PdbSiteCard::PdbSiteCard(stringstream &stream_block)
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
        PdbSite* site = new PdbSite(site_block);
        residue_sites_[site->GetSiteId()] = site;
    }
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////

string PdbSiteCard::GetRecordName(){
    return record_name_;
}

PdbSiteCard::PdbSiteMap PdbSiteCard::GetResidueSites(){
    return residue_sites_;
}

//////////////////////////////////////////////////////////
//                       MUTATOR                        //
//////////////////////////////////////////////////////////

void PdbSiteCard::SetRecordName(const string record_name){
    record_name_ = record_name;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbSiteCard::Print(ostream &out)
{
    out << "Record Name: " << record_name_ << endl <<
           "==================== Residue Site ==============" << endl;
    for(PdbSiteCard::PdbSiteMap::iterator it = residue_sites_.begin(); it != residue_sites_.end(); it++)
    {
        out << "Site ID: " << (it)->first << endl;
        (it)->second->Print();
        out << endl;
    }
}
