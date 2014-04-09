
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
int PdbPreprocessorMissingResidue::GetStartingResidueSequenceNumber()
{
    return starting_residue_sequence_number_;
}
int PdbPreprocessorMissingResidue::GetEndingResidueSequenceNumber()
{
    return ending_residue_sequence_number_;
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
string PdbPreprocessorMissingResidue::GetSelectedNTermination()
{
    return selected_n_termination_;
}
string PdbPreprocessorMissingResidue::GetSelectedCTermination()
{
    return selected_c_termination_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbPreprocessorMissingResidue::SetResidueChainId(char residue_chain_id)
{
    residue_chain_id_ = residue_chain_id;
}
void PdbPreprocessorMissingResidue::SetStartingResidueSequenceNumber(int starting_residue_sequence_number)
{
    starting_residue_sequence_number_ = starting_residue_sequence_number;
}
void PdbPreprocessorMissingResidue::SetEndingResidueSequenceNumber(int ending_residue_sequence_number)
{
    ending_residue_sequence_number_ = ending_residue_sequence_number;
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
void PdbPreprocessorMissingResidue::SetSelectedNTermination(const string selected_n_termination)
{
    selected_n_termination_ = selected_n_termination;
}
void PdbPreprocessorMissingResidue::SetSelectedCTermination(const string selected_c_termination)
{
    selected_c_termination_ = selected_c_termination;
}

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbPreprocessorMissingResidue::Print(ostream &out)
{
}







