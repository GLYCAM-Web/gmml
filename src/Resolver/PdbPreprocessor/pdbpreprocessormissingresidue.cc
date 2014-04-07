
#include "../../../includes/Resolver/PdbPreprocessor/pdbpreprocessormissingresidue.hpp"

using namespace std;
using namespace PdbPreprocessorSpace;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbPreprocessorMissingResidue::PdbPreprocessorMissingResidue() {}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
char PdbPreprocessorMissingResidue::GetResidueChainId()
{
    return residue_chain_id_;
}
int PdbPreprocessorMissingResidue::GetResidueStartIndex()
{
    return residue_start_index_;
}
int PdbPreprocessorMissingResidue::GetResidueEndIndex()
{
    return residue_end_index_;
}
int PdbPreprocessorMissingResidue::GetResidueBeforeGap()
{
    return residue_before_gap_;
}
int PdbPreprocessorMissingResidue::GetResidueAfterGap()
{
    return residue_after_gap_;
}
vector<string> PdbPreprocessorMissingResidue::GetPossibleNTerminations()
{
    return possible_n_terminations_;
}
vector<string> PdbPreprocessorMissingResidue::GetPossibleCTerminations()
{
    return possible_c_terminations_;
}
string PdbPreprocessorMissingResidue::GetSelectedNTerminations()
{
    return selected_n_terminations_;
}
string PdbPreprocessorMissingResidue::GetSelectedCTerminations()
{
    return selected_c_terminations_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbPreprocessorMissingResidue::SetResidueChainId(char residue_chain_id)
{
    residue_chain_id_ = residue_chain_id;
}
void PdbPreprocessorMissingResidue::SetResidueStartIndex(int residue_start_index)
{
    residue_start_index_ = residue_start_index;
}
void PdbPreprocessorMissingResidue::SetResidueEndIndex(int residue_end_index)
{
    residue_end_index_ = residue_end_index;
}
void PdbPreprocessorMissingResidue::SetResidueBeforeGap(int residue_before_gap)
{
    residue_before_gap_ = residue_before_gap;
}

void PdbPreprocessorMissingResidue::SetResidueAfterGap(int residue_after_gap)
{
    residue_after_gap_ = residue_after_gap;
}
void PdbPreprocessorMissingResidue::SetPossibleNTerminations(vector<string> possible_n_terminations)
{
    possible_n_terminations_.clear();
    for(vector<string>::iterator it = possible_n_terminations.begin(); it != possible_n_terminations.end(); it++)
    {
        possible_n_terminations_.push_back(*it);
    }
}
void PdbPreprocessorMissingResidue::SetPossibleCTerminations(vector<string> possible_c_terminations)
{
    possible_c_terminations_.clear();
    for(vector<string>::iterator it = possible_c_terminations.begin(); it != possible_c_terminations.end(); it++)
    {
        possible_c_terminations_.push_back(*it);
    }
}
void PdbPreprocessorMissingResidue::SetSelectedNTerminations(const string selected_n_terminations)
{
    selected_n_terminations_ = selected_n_terminations;
}
void PdbPreprocessorMissingResidue::SetSelectedCTerminations(const string selected_c_terminations)
{
    selected_c_terminations_ = selected_c_terminations;
}

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbPreprocessorMissingResidue::Print(ostream &out)
{
}







