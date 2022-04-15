#include "../../../includes/Resolver/PdbPreprocessor/pdbpreprocessorunrecognizedresidue.hpp"

using PdbPreprocessorSpace::PdbPreprocessorUnrecognizedResidue;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbPreprocessorUnrecognizedResidue::PdbPreprocessorUnrecognizedResidue() {}

PdbPreprocessorUnrecognizedResidue::PdbPreprocessorUnrecognizedResidue(
        std::string residue_name, char chain_id, int sequence_number, char residue_insertion_code,
        char residue_alternate_location, bool middle_of_chain) :
    residue_chain_id_(chain_id), residue_sequence_number_(sequence_number),
    residue_name_(residue_name), residue_insertion_code_(residue_insertion_code),
    residue_alternate_location_(residue_alternate_location), middle_of_chain_(middle_of_chain) {}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
char PdbPreprocessorUnrecognizedResidue::GetResidueChainId()
{
    return residue_chain_id_;
}
int PdbPreprocessorUnrecognizedResidue::GetResidueSequenceNumber()
{
    return residue_sequence_number_;
}
std::string PdbPreprocessorUnrecognizedResidue::GetResidueName()
{
    return residue_name_;
}
char PdbPreprocessorUnrecognizedResidue::GetResidueInsertionCode()
{
    return residue_insertion_code_;
}
char PdbPreprocessorUnrecognizedResidue::GetResidueAlternateLocation()
{
    return residue_alternate_location_;
}
bool PdbPreprocessorUnrecognizedResidue::GetMiddleOfChain()
{
    return middle_of_chain_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbPreprocessorUnrecognizedResidue::SetResidueChainId(char residue_chain_id)
{
    residue_chain_id_ = residue_chain_id;
}
void PdbPreprocessorUnrecognizedResidue::SetResidueSequenceNumber(int residue_sequence_number)
{
    residue_sequence_number_ = residue_sequence_number;
}
void PdbPreprocessorUnrecognizedResidue::SetResidueName(const std::string residue_name)
{
    residue_name_ = residue_name;
}
void PdbPreprocessorUnrecognizedResidue::SetResidueInsertionCode(char residue_insertion_code)
{
    residue_insertion_code_ = residue_insertion_code;
}
void PdbPreprocessorUnrecognizedResidue::SetResidueAlternateLocation(char residue_alternate_location)
{
    residue_alternate_location_ = residue_alternate_location;
}
void PdbPreprocessorUnrecognizedResidue::SetMiddleOfChain(bool middle_of_chain)
{
    middle_of_chain_ = middle_of_chain;
}

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbPreprocessorUnrecognizedResidue::Print(std::ostream &out)
{
    out << "Residue name: " << residue_name_
         << ", Chain id: " << residue_chain_id_
         << ", Sequence number: " << residue_sequence_number_
         << ", Insertion code: " << residue_insertion_code_
         << ", Alternate location: " << residue_alternate_location_
         << ", In the middle of chain: " << middle_of_chain_
         << std::endl;
}
