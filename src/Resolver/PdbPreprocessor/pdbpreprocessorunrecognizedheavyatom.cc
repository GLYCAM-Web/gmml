#include "../../../includes/Resolver/PdbPreprocessor/pdbpreprocessorunrecognizedheavyatom.hpp"

using PdbPreprocessorSpace::PdbPreprocessorUnrecognizedHeavyAtom;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbPreprocessorUnrecognizedHeavyAtom::PdbPreprocessorUnrecognizedHeavyAtom() {}

PdbPreprocessorUnrecognizedHeavyAtom::PdbPreprocessorUnrecognizedHeavyAtom(
        char residue_chain_id, int atom_serial_number, std::string atom_name, std::string residue_name,
        int residue_sequence_number, char residue_insertion_code, char residue_alternate_location) :
    residue_chain_id_(residue_chain_id), atom_serial_number_(atom_serial_number),
    atom_name_(atom_name), residue_name_(residue_name), residue_sequence_number_(residue_sequence_number),
    residue_insertion_code_(residue_insertion_code), residue_alternate_location_(residue_alternate_location) {}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
char PdbPreprocessorUnrecognizedHeavyAtom::GetResidueChainId()
{
    return residue_chain_id_;
}
int PdbPreprocessorUnrecognizedHeavyAtom::GetAtomSerialNumber()
{
    return atom_serial_number_;
}
std::string PdbPreprocessorUnrecognizedHeavyAtom::GetAtomName()
{
    return atom_name_;
}
std::string PdbPreprocessorUnrecognizedHeavyAtom::GetResidueName()
{
    return residue_name_;
}
int PdbPreprocessorUnrecognizedHeavyAtom::GetResidueSequenceNumber()
{
    return residue_sequence_number_;
}
char PdbPreprocessorUnrecognizedHeavyAtom::GetResidueInsertionCode()
{
    return residue_insertion_code_;
}
char PdbPreprocessorUnrecognizedHeavyAtom::GetResidueAlternateLocation()
{
    return residue_alternate_location_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbPreprocessorUnrecognizedHeavyAtom::SetResidueChainId(char residue_chain_id)
{
    residue_chain_id_ = residue_chain_id;
}
void PdbPreprocessorUnrecognizedHeavyAtom::SetAtomSerialNumber(int atom_serial_number)
{
    atom_serial_number_ = atom_serial_number;
}
void PdbPreprocessorUnrecognizedHeavyAtom::SetAtomName(const std::string atom_name)
{
    atom_name_ = atom_name;
}
void PdbPreprocessorUnrecognizedHeavyAtom::SetResidueName(const std::string residue_name)
{
    residue_name_ = residue_name;
}
void PdbPreprocessorUnrecognizedHeavyAtom::SetResidueSequenceNumber(int residue_sequence_number)
{
    residue_sequence_number_ = residue_sequence_number;
}
void PdbPreprocessorUnrecognizedHeavyAtom::SetResidueInsertionCode(char residue_insertion_code)
{
    residue_insertion_code_ = residue_insertion_code;
}
void PdbPreprocessorUnrecognizedHeavyAtom::SetResidueAlternateLocation(char residue_alternate_location)
{
    residue_alternate_location_ = residue_alternate_location;
}

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbPreprocessorUnrecognizedHeavyAtom::Print(std::ostream &out)
{
    out << "Atom name: " << atom_name_
         << ", Serial number: " << atom_serial_number_
         << ", Chain id: " << residue_chain_id_
         << ", Residue name: " << residue_name_
         << ", Sequence number: " << residue_sequence_number_
         << ", Insertion code: " << residue_insertion_code_
         << ", Alternate location: " << residue_alternate_location_
         << std::endl;
}
