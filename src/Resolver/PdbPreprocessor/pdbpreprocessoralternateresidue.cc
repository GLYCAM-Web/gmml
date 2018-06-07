#include "../../../includes/Resolver/PdbPreprocessor/pdbpreprocessoralternateresidue.hpp"

using PdbPreprocessorSpace::PdbPreprocessorAlternateResidue;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbPreprocessorAlternateResidue::PdbPreprocessorAlternateResidue() {}

PdbPreprocessorAlternateResidue::PdbPreprocessorAlternateResidue(
        std::string residue_name, char chain_id, int sequence_number, char residue_insertion_code,
        std::vector<char> residue_alternate_location, std::vector<bool> selected_alternate_location) :
    residue_chain_id_(chain_id), residue_sequence_number_(sequence_number), residue_name_(residue_name), 
    residue_insertion_code_(residue_insertion_code)
{
    residue_alternate_location_ = std::vector<char>();
    for(std::vector<char>::iterator it = residue_alternate_location.begin(); it != residue_alternate_location.end(); it++)
    {
        residue_alternate_location_.push_back(*it);
    }
    selected_alternate_location_ = std::vector<bool>();
    for(std::vector<bool>::iterator it = selected_alternate_location.begin(); it != selected_alternate_location.end(); it++)
    {
        selected_alternate_location_.push_back(*it);
    }
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
char PdbPreprocessorAlternateResidue::GetResidueChainId()
{
    return residue_chain_id_;
}
int PdbPreprocessorAlternateResidue::GetResidueSequenceNumber()
{
    return residue_sequence_number_;
}
std::string PdbPreprocessorAlternateResidue::GetResidueName()
{
    return residue_name_;
}
char PdbPreprocessorAlternateResidue::GetResidueInsertionCode()
{
    return residue_insertion_code_;
}
std::vector<char> PdbPreprocessorAlternateResidue::GetResidueAlternateLocation()
{
    return residue_alternate_location_;
}
std::vector<bool> PdbPreprocessorAlternateResidue::GetSelectedAlternateLocation()
{
    return selected_alternate_location_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbPreprocessorAlternateResidue::SetResidueChainId(char residue_chain_id)
{
    residue_chain_id_ = residue_chain_id;
}
void PdbPreprocessorAlternateResidue::SetResidueSequenceNumber(int residue_sequence_number)
{
    residue_sequence_number_ = residue_sequence_number;
}
void PdbPreprocessorAlternateResidue::SetResidueName(const std::string residue_name)
{
    residue_name_ = residue_name;
}
void PdbPreprocessorAlternateResidue::SetResidueInsertionCode(char residue_insertion_code)
{
    residue_insertion_code_ = residue_insertion_code;
}
void PdbPreprocessorAlternateResidue::SetResidueAlternateLocation(std::vector<char> residue_alternate_location)
{
    residue_alternate_location_.clear();
    for(std::vector<char>::iterator it = residue_alternate_location.begin(); it != residue_alternate_location.end(); it++)
    {
        residue_alternate_location_.push_back(*it);
    }
}
void PdbPreprocessorAlternateResidue::SetSelectedAlternateLocation(std::vector<bool> selected_alternate_location)
{
    selected_alternate_location_.clear();
    for(std::vector<bool>::iterator it = selected_alternate_location.begin(); it != selected_alternate_location.end(); it++)
    {
        selected_alternate_location_.push_back(*it);
    }
}

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbPreprocessorAlternateResidue::Print(std::ostream &out)
{
    out << "Residue name: " << residue_name_
         << ", Chain id: " << residue_chain_id_
         << ", Sequence number: " << residue_sequence_number_
         << ", Insertion code: " << residue_insertion_code_ << std::endl;
         for(unsigned int i = 0; i < residue_alternate_location_.size(); i++)
         {
             out << "\t Alternate Location: " << residue_alternate_location_[i] << " " << selected_alternate_location_[i] << std::endl;
         }

         out << std::endl;
}
