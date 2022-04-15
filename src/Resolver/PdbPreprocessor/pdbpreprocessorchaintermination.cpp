#include "../../../includes/Resolver/PdbPreprocessor/pdbpreprocessorchaintermination.hpp"

using PdbPreprocessorSpace::PdbPreprocessorChainTermination;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbPreprocessorChainTermination::PdbPreprocessorChainTermination() {}

PdbPreprocessorChainTermination::PdbPreprocessorChainTermination(
        char chain_id, int starting_sequence_number, int ending_sequence_number,
        char starting_residue_insertion_code, char ending_residue_insertion_code) :
    residue_chain_id_(chain_id), starting_residue_sequence_number_(starting_sequence_number),
    ending_residue_sequence_number_(ending_sequence_number),
    starting_residue_insertion_code_(starting_residue_insertion_code),
    ending_residue_insertion_code_(ending_residue_insertion_code)
{
    selected_n_termination_ = (gmml::PossibleNChainTermination)2;
    selected_c_termination_ = (gmml::PossibleCChainTermination)3;
}

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
gmml::PossibleNChainTermination PdbPreprocessorChainTermination::GetSelectedNTermination()
{
    return selected_n_termination_;
}
gmml::PossibleCChainTermination PdbPreprocessorChainTermination::GetSelectedCTermination()
{
    return selected_c_termination_;
}
char PdbPreprocessorChainTermination::GetStartingResidueInsertionCode()
{
    return starting_residue_insertion_code_;
}
char PdbPreprocessorChainTermination::GetEndingResidueInsertionCode()
{
    return ending_residue_insertion_code_;
}

std::string PdbPreprocessorChainTermination::GetStringFormatOfSelectedNTermination()
{
    switch (selected_n_termination_)
    {
        case 1:
            return "ACE";
        case 2:
            return "gmml::NH3+";
        default:
            return "";
    }
}

std::string PdbPreprocessorChainTermination::GetStringFormatOfSelectedCTermination()
{
    switch(selected_c_termination_)
    {
        case 1:
            return "NHE";
        case 2:
            return "NME";
        case 3:
            return "gmml::CO2-";
        default:
            return "";
    }
}

std::string PdbPreprocessorChainTermination::GetStringFormatOfNTermination(gmml::PossibleNChainTermination n_termination)
{
    switch (n_termination)
    {
        case 1:
            return "ACE";
        case 2:
            return "gmml::NH3+";
        default:
            return "";
    }
}

std::string PdbPreprocessorChainTermination::GetStringFormatOfCTermination(gmml::PossibleCChainTermination c_termination)
{
    switch(c_termination)
    {
        case 1:
            return "NHE";
        case 2:
            return "NME";
        case 3:
            return "gmml::CO2-";
        default:
            return "";
    }
}

std::vector<std::string> PdbPreprocessorChainTermination::GetAllPossibleNChainTerminationAsString()
{
    std::vector<std::string> all_possible_n_chian_termination_as_string;
    for(int possible_n_chian_termination = gmml::COCH3; possible_n_chian_termination != gmml::NH3; possible_n_chian_termination++)
        all_possible_n_chian_termination_as_string.push_back(GetStringFormatOfNTermination((gmml::PossibleNChainTermination)possible_n_chian_termination));
    return all_possible_n_chian_termination_as_string;
}

std::vector<std::string> PdbPreprocessorChainTermination::GetAllPossibleCChainTerminationAsString()
{
    std::vector<std::string> all_possible_c_chian_termination_as_string;
    for(int possible_c_chian_termination = gmml::NH2; possible_c_chian_termination != gmml::CO2; possible_c_chian_termination++)
        all_possible_c_chian_termination_as_string.push_back(GetStringFormatOfCTermination((gmml::PossibleCChainTermination)possible_c_chian_termination));
    return all_possible_c_chian_termination_as_string;
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
void PdbPreprocessorChainTermination::SetSelectedNTermination(gmml::PossibleNChainTermination selected_n_termination)
{
    selected_n_termination_ = selected_n_termination;
}
void PdbPreprocessorChainTermination::SetSelectedCTermination(gmml::PossibleCChainTermination selected_c_termination)
{
    selected_c_termination_ = selected_c_termination;
}
void PdbPreprocessorChainTermination::SetStartingResidueInsertionCode(char starting_residue_insertion_code)
{
    starting_residue_insertion_code_ = starting_residue_insertion_code;
}
void PdbPreprocessorChainTermination::SetEndingResidueInsertionCode(char ending_residue_insertion_code)
{
    ending_residue_insertion_code_ = ending_residue_insertion_code;
}

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbPreprocessorChainTermination::Print(std::ostream &out)
{
    out << "Chain id " << residue_chain_id_
         << ", Start sequence number: " << starting_residue_sequence_number_
         << ", End sequence number: " << ending_residue_sequence_number_
         << ", Start insertion code: " << starting_residue_insertion_code_
         << ", End insertion code: " << ending_residue_insertion_code_
         << ", Selected n termination: " << GetStringFormatOfSelectedNTermination()
         << ", Selected c termination: " << GetStringFormatOfSelectedCTermination()
         << std::endl;
}
