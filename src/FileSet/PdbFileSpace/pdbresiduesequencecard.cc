// Author: Alireza Khatamian

#include "../../../includes/FileSet/PdbFileSpace/pdbresiduesequencecard.hpp"

using namespace std;
using namespace PdbFileSpace;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbResidueSequenceCard::PdbResidueSequenceCard() : record_name_("SEQRES") {}

PdbResidueSequenceCard::PdbResidueSequenceCard(const string &record_name) : record_name_(record_name) {}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
string PdbResidueSequenceCard::GetRecordName()
{
    return record_name_;
}

PdbResidueSequenceCard::ResidueSequenceMap PdbResidueSequenceCard::GetResidueSequenceChain()
{
    return residue_sequence_chains_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbResidueSequenceCard::SetRecordName(const string record_name)
{
    record_name_ = record_name;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////

