
#include "../../../includes/Resolver/PdbPreprocessor/pdbpreprocessorhistidinemapping.hpp"

using namespace std;
using namespace PdbPreprocessorSpace;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbPreprocessorHistidineMapping::PdbPreprocessorHistidineMapping() {}

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
vector<string> PdbPreprocessorHistidineMapping::GetPossibleMappings()
{
    return possible_mappings_;
}
string PdbPreprocessorHistidineMapping::GetSelectedMapping()
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
void PdbPreprocessorHistidineMapping::SetPossibleMappings(vector<string> possible_mappings)
{
    possible_mappings_.clear();
    for(vector<string>::iterator it = possible_mappings.begin(); it != possible_mappings.end(); it++)
    {
        possible_mappings_.push_back(*it);
    }
}
void PdbPreprocessorHistidineMapping::SetSelectedMapping(string selected_mapping)
{
    selected_mapping_ = selected_mapping;
}

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbPreprocessorHistidineMapping::Print(ostream &out)
{
}







