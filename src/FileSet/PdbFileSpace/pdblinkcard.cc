#include "../../../includes/FileSet/PdbFileSpace/pdblinkcard.hpp"
#include "../../../includes/FileSet/PdbFileSpace/pdblink.hpp"
#include "../../../includes/utils.hpp"

using namespace std;
using namespace PdbFileSpace;
using namespace gmml;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbLinkCard::PdbLinkCard() {}

PdbLinkCard::PdbLinkCard(stringstream &stream_block)
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

        PdbLink* link = new PdbLink(line);
        AddResidueLink(link);
        getline(stream_block, line);
    }
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////

string PdbLinkCard::GetRecordName(){
    return record_name_;
}

PdbLinkCard::LinkVector PdbLinkCard::GetResidueLinks(){
    return residue_links_;
}

//////////////////////////////////////////////////////////
//                       MUTATOR                        //
//////////////////////////////////////////////////////////

void PdbLinkCard::SetRecordName(string record_name){
    record_name_ = record_name;
}

void PdbLinkCard::SetResidueLinks(const LinkVector residue_links){
    residue_links_ = residue_links;
}

void PdbLinkCard::AddResidueLink(PdbLink *residue_link)
{
    residue_links_.push_back(residue_link);
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////






