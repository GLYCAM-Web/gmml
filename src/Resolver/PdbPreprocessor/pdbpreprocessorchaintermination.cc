
#include "../../../includes/Resolver/PdbPreprocessor/pdbpreprocessorchaintermination.hpp"

using namespace std;
using namespace PdbPreprocessorSpace;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbPreprocessorChainTermination::PdbPreprocessorChainTermination() {}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
char PdbPreprocessorChainTermination::GetResidueChainId()
{
    return residue_chain_id_;
}
int PdbPreprocessorChainTermination::GetResidueStartIndex()
{
    return residue_start_index_;
}
int PdbPreprocessorChainTermination::GetResidueEndIndex()
{
    return residue_end_index_;
}
vector<string> PdbPreprocessorChainTermination::GetPossibleNTerminations()
{
    return possible_n_terminations_;
}
vector<string> PdbPreprocessorChainTermination::GetPossibleCTerminations()
{
    return possible_c_terminations_;
}
string PdbPreprocessorChainTermination::GetSelectedNTerminations()
{
    return selected_n_terminations_;
}
string PdbPreprocessorChainTermination::GetSelectedCTerminations()
{
    return selected_c_terminations_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbPreprocessorChainTermination::SetResidueChainId(char residue_chain_id)
{
    residue_chain_id_ = residue_chain_id;
}
void PdbPreprocessorChainTermination::SetResidueStartIndex(int residue_start_index)
{
    residue_start_index_ = residue_start_index;
}
void PdbPreprocessorChainTermination::SetResidueEndIndex(int residue_end_index)
{
    residue_end_index_ = residue_end_index;
}
void PdbPreprocessorChainTermination::SetPossibleNTerminations(vector<string> possible_n_terminations)
{
    possible_n_terminations_.clear();
    for(vector<string>::iterator it = possible_n_terminations.begin(); it != possible_n_terminations.end(); it++)
    {
        possible_n_terminations_.push_back(*it);
    }
}
void PdbPreprocessorChainTermination::SetPossibleCTerminations(vector<string> possible_c_terminations)
{
    possible_c_terminations_.clear();
    for(vector<string>::iterator it = possible_c_terminations.begin(); it != possible_c_terminations.end(); it++)
    {
        possible_c_terminations_.push_back(*it);
    }
}
void PdbPreprocessorChainTermination::SetSelectedNTerminations(const string selected_n_terminations)
{
    selected_n_terminations_ = selected_n_terminations;
}
void PdbPreprocessorChainTermination::SetSelectedCTerminations(const string selected_c_terminations)
{
    selected_c_terminations_ = selected_c_terminations;
}

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbPreprocessorChainTermination::Print(ostream &out)
{
}







