// Author: Alireza Khatamian

#include "../../../includes/FileSet/PdbFileSpace/pdbresiduemodificationcard.hpp"

using namespace std;
using namespace PdbFileSpace;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbResidueModificationCard::PdbResidueModificationCard() : record_name_("MODRES") {}

PdbResidueModificationCard::PdbResidueModificationCard(const string &record_name) : record_name_(record_name) {}

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

