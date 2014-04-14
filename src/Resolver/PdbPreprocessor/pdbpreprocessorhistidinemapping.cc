
#include "../../../includes/Resolver/PdbPreprocessor/pdbpreprocessorhistidinemapping.hpp"

using namespace std;
using namespace PdbPreprocessorSpace;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbPreprocessorHistidineMapping::PdbPreprocessorHistidineMapping() {}

PdbPreprocessorHistidineMapping::PdbPreprocessorHistidineMapping(char residue_chain_id, int residue_sequence_number, PdbPreprocessorHISMapping selected_mapping) :
    residue_chain_id_(residue_chain_id), residue_sequence_number_(residue_sequence_number), selected_mapping_(selected_mapping) {}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
char PdbPreprocessorHistidineMapping::GetResidueChainId()
{
    return residue_chain_id_;
}
int PdbPreprocessorHistidineMapping::GetResidueSequenceNumber()
{
    return residue_sequence_number_;
}
PdbPreprocessorHISMapping PdbPreprocessorHistidineMapping::GetSelectedMapping()
{
    return selected_mapping_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbPreprocessorHistidineMapping::SetResidueChainId(char residue_chain_id)
{
    residue_chain_id_ = residue_chain_id;
}
void PdbPreprocessorHistidineMapping::SetResidueSequenceNumber(int residue_sequence_number)
{
    residue_sequence_number_ = residue_sequence_number;
}
void PdbPreprocessorHistidineMapping::SetSelectedMapping(PdbPreprocessorHISMapping selected_mapping)
{
    selected_mapping_ = selected_mapping;
}

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbPreprocessorHistidineMapping::Print(ostream &out)
{
}







