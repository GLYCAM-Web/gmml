#include "../../../includes/Resolver/PdbPreprocessor/pdbpreprocessorreplacedhydrogen.hpp"

using PdbPreprocessorSpace::PdbPreprocessorReplacedHydrogen;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbPreprocessorReplacedHydrogen::PdbPreprocessorReplacedHydrogen() {}

PdbPreprocessorReplacedHydrogen::PdbPreprocessorReplacedHydrogen(
        char residue_chain_id, int atom_serial_number, std::string atom_name,
        std::string residue_name, int residue_sequence_number, char residue_insertion_code,
        char residue_alternate_location) :
    atom_serial_number_(atom_serial_number), atom_name_(atom_name), residue_name_(residue_name),
    residue_sequence_number_(residue_sequence_number), residue_chain_id_(residue_chain_id),
    residue_insertion_code_(residue_insertion_code), residue_alternate_location_(residue_alternate_location) {}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
int PdbPreprocessorReplacedHydrogen::GetAtomSerialNumber()
{
    return atom_serial_number_;
}
std::string PdbPreprocessorReplacedHydrogen::GetAtomName()
{
    return atom_name_;
}
std::string PdbPreprocessorReplacedHydrogen::GetResidueName()
{
    return residue_name_;
}
int PdbPreprocessorReplacedHydrogen::GetResidueSequenceNumber()
{
    return residue_sequence_number_;
}
char PdbPreprocessorReplacedHydrogen::GetResidueChainId()
{
    return residue_chain_id_;
}
char PdbPreprocessorReplacedHydrogen::GetResidueInsertionCode()
{
    return residue_insertion_code_;
}
char PdbPreprocessorReplacedHydrogen::GetResidueAlternateLocation()
{
    return residue_alternate_location_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbPreprocessorReplacedHydrogen::SetAtomSerialNumber(int atom_serial_number)
{
    atom_serial_number_ = atom_serial_number;
}
void PdbPreprocessorReplacedHydrogen::SetAtomName(const std::string atom_name)
{
    atom_name_ = atom_name;
}
void PdbPreprocessorReplacedHydrogen::SetResidueName(const std::string residue_name)
{
    residue_name_ = residue_name;
}
void PdbPreprocessorReplacedHydrogen::SetResidueSequenceNumber(int residue_sequence_number)
{
    residue_sequence_number_ = residue_sequence_number;
}
void PdbPreprocessorReplacedHydrogen::SetResidueChainId(char residue_chain_id)
{
    residue_chain_id_ = residue_chain_id;
}
void PdbPreprocessorReplacedHydrogen::SetResidueInsertionCode(char residue_insertion_code)
{
    residue_insertion_code_ = residue_insertion_code;
}
void PdbPreprocessorReplacedHydrogen::SetResidueAlternateLocation(char residue_alternate_location)
{
    residue_alternate_location_ = residue_alternate_location;
}

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbPreprocessorReplacedHydrogen::Print(std::ostream &out)
{
    out << "Atom name: " << atom_name_
         << ", Serial number: " << atom_serial_number_
         << ", Residue name: " << residue_name_
         << ", Sequence number: " << residue_sequence_number_
         << ", chain id: " << residue_chain_id_
         << ", insertion code: " << residue_insertion_code_
         << ", Alternate location: " << residue_alternate_location_
         << std::endl;
}
