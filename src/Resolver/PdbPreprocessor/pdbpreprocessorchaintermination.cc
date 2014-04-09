
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
int PdbPreprocessorChainTermination::GetStartingResidueSequenceNumber()
{
    return starting_residue_sequence_number_;
}
int PdbPreprocessorChainTermination::GetEndingResidueSequenceNumber()
{
    return ending_residue_sequence_number_;
}
vector<string> PdbPreprocessorChainTermination::GetPossibleNTerminations()
{
    return possible_n_terminations_;
}
vector<string> PdbPreprocessorChainTermination::GetPossibleCTerminations()
{
    return possible_c_terminations_;
}
string PdbPreprocessorChainTermination::GetSelectedNTermination()
{
    return selected_n_termination_;
}
string PdbPreprocessorChainTermination::GetSelectedCTermination()
{
    return selected_c_termination_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbPreprocessorChainTermination::SetResidueChainId(char residue_chain_id)
{
    residue_chain_id_ = residue_chain_id;
}
void PdbPreprocessorChainTermination::SetStartingResidueSequenceNumber(int starting_residue_sequence_number)
{
    starting_residue_sequence_number_ = starting_residue_sequence_number;
}
void PdbPreprocessorChainTermination::SetEndingResidueSequenceNumber(int ending_residue_sequence_number)
{
    ending_residue_sequence_number_ = ending_residue_sequence_number;
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
void PdbPreprocessorChainTermination::SetSelectedNTermination(const string selected_n_termination)
{
    selected_n_termination_ = selected_n_termination;
}
void PdbPreprocessorChainTermination::SetSelectedCTermination(const string selected_c_termination)
{
    selected_c_termination_ = selected_c_termination;
}

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbPreprocessorChainTermination::Print(ostream &out)
{
}







