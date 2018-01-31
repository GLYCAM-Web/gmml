// Author: Alireza Khatamian

#include "../../../includes/InputSet/PdbFileSpace/pdbresiduemodificationcard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbresiduemodificationsection.hpp"
#include "../../../includes/utils.hpp"

using namespace std;
using namespace PdbFileSpace;
using namespace gmml;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbResidueModificationSection::PdbResidueModificationSection() : record_name_("MODRES") {}

PdbResidueModificationSection::PdbResidueModificationSection(const string &record_name) : record_name_(record_name) {}

PdbResidueModificationSection::PdbResidueModificationSection(stringstream &stream_block)
{
    string line;
    stringstream ss;
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

        ss << line << endl;
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
void PdbResidueModificationSection::SetRecordName(const string record_name)
{
    record_name_ = record_name;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbResidueModificationSection::Print(ostream &out)
{
    out << "Record Name: " << record_name_ << endl <<
           "============ Residue Modification ===========" << endl;
    for(PdbResidueModificationSection::ResidueModificationCardMap::iterator it = residue_modification_cards_.begin(); it != residue_modification_cards_.end(); it++)
    {
        out << "ID Code: " << (it)->first << endl;
        (it)->second->Print();
        out << endl;
    }
}
