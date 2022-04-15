#include "../../../includes/InputSet/PdbFileSpace/pdbresidue.hpp"

using PdbFileSpace::PdbResidue;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbResidue::PdbResidue() {}

PdbResidue::PdbResidue( std::string residue_name, char chain_id,
                        int sequence_number, char insertion_code,
                        char residue_alternate_location) :
    residue_chain_id_(chain_id), residue_name_(residue_name), 
    residue_sequence_number_(sequence_number), residue_insertion_code_(insertion_code),
    residue_alternate_location_(residue_alternate_location) {}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
char PdbResidue::GetResidueChainId()
{
    return residue_chain_id_;
}
std::string PdbResidue::GetResidueName()
{
    return residue_name_;
}
int PdbResidue::GetResidueSequenceNumber()
{
    return residue_sequence_number_;
}
char PdbResidue::GetResidueInsertionCode()
{
    return residue_insertion_code_;
}
char PdbResidue::GetResidueAlternateLocation()
{
    return residue_alternate_location_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbResidue::SetResidueChainId(char residue_chain_id)
{
    residue_chain_id_ = residue_chain_id;
}
void PdbResidue::SetResidueName(const std::string residue_name)
{
    residue_name_ = residue_name;
}
void PdbResidue::SetResidueSequenceNumber(int residue_sequence_number)
{
    residue_sequence_number_ = residue_sequence_number;
}
void PdbResidue::SetResidueInsertionCode(char residue_insertion_code)
{
    residue_insertion_code_ = residue_insertion_code;
}
void PdbResidue::SetResidueAlternateLocation(char residue_alternate_location)
{
    residue_alternate_location_ = residue_alternate_location;
}

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbResidue::Print(std::ostream &out)
{
    out << "Residue name: " << residue_name_
         << ", Chain id: " << residue_chain_id_
         << ", Sequence number: " << residue_sequence_number_
         << ", Insertion code: " << residue_insertion_code_
         << ", Alternate location: " << residue_alternate_location_ << std::endl;
}
