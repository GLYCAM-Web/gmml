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
    line = Trim(line);
    while (!Trim(line).empty())
    {
        if(!is_record_name_set){
            record_name_ = line.substr(0,6);
            is_record_name_set=true;
        }

        ss << line;
        PdbResidueModification* residue_modification = new PdbResidueModification(ss);
        residue_modifications_[line.substr(7,4)] = residue_modification;
        getline(stream_block, line);
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

