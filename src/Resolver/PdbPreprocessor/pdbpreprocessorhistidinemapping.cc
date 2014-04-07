
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
int PdbPreprocessorHistidineMapping::GetResidueNumber()
{
    return residue_number_;
}
vector<string> PdbPreprocessorHistidineMapping::GetPossibleMappings()
{
    return possible_mappings_;
}
string PdbPreprocessorHistidineMapping::GetSelectedMappings()
{
    return selected_mappings_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbPreprocessorHistidineMapping::SetResidueChainId(char residue_chain_id)
{
    residue_chain_id_ = residue_chain_id;
}
void PdbPreprocessorHistidineMapping::SetResidueNumber(int residue_number)
{
    residue_number_ = residue_number;
}
void PdbPreprocessorHistidineMapping::SetPossibleMappings(vector<string> possible_mappings)
{
    possible_mappings_.clear();
    for(vector<string>::iterator it = possible_mappings.begin(); it != possible_mappings.end(); it++)
    {
        possible_mappings_.push_back(*it);
    }
}
void PdbPreprocessorHistidineMapping::SetSelectedMappings(string selected_mappings)
{
    selected_mappings_ = selected_mappings;
}

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbPreprocessorHistidineMapping::Print(ostream &out)
{
}







