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
    line = Trim(line);
    while (!Trim(line).empty())
    {
        if(!is_record_name_set){
            record_name_ = line.substr(0,6);
            is_record_name_set=true;
        }
        stringstream site_block;
        site_block << line << endl;
        string site_id = line.substr(11,3);

        getline(stream_block, line);

        while (!Trim(line).empty() && line.substr(11,3) == site_id){
            site_block << line << endl;
            getline(stream_block, line);
        }
        PdbSite* site = new PdbSite(site_block);
        residue_sites_[site_id] = site;
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



