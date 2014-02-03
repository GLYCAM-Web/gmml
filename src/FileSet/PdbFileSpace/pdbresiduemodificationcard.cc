// Author: Alireza Khatamian

#include "../../../includes/FileSet/PdbFileSpace/pdbresiduemodificationcard.hpp"
#include "../../../includes/FileSet/PdbFileSpace/pdbresiduemodification.hpp"
#include "../../../includes/utils.hpp"

using namespace std;
using namespace PdbFileSpace;
using namespace gmml;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbResidueModificationCard::PdbResidueModificationCard() : record_name_("MODRES") {}

PdbResidueModificationCard::PdbResidueModificationCard(const string &record_name) : record_name_(record_name) {}

PdbResidueModificationCard::PdbResidueModificationCard(stringstream &stream_block)
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
        PdbResidueModification* residue_modification = new PdbResidueModification(ss);
        residue_modifications_[residue_modification->GetIdCode()] = residue_modification;
        getline(stream_block, line);
        temp = line;
    }
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
std::string PdbResidueModificationCard::GetRecordName()
{
    return record_name_;
}

PdbResidueModificationCard::ResidueModificationMap PdbResidueModificationCard::GetResidueModifications()
{
    return residue_modifications_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbResidueModificationCard::SetRecordName(const string record_name)
{
    record_name_ = record_name;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbResidueModificationCard::Print(ostream &out)
{
    out << "Record Name: " << record_name_ << endl <<
           "============ Residue Modification ===========" << endl;
    for(PdbResidueModificationCard::ResidueModificationMap::iterator it = residue_modifications_.begin(); it != residue_modifications_.end(); it++)
    {
        out << "ID Code: " << (it)->first << endl;
        (it)->second->Print();
        out << endl;
    }
}
