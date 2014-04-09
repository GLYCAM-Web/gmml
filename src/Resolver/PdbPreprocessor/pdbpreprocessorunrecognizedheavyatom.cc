
#include "../../../includes/Resolver/PdbPreprocessor/pdbpreprocessorunrecognizedheavyatom.hpp"

using namespace std;
using namespace PdbPreprocessorSpace;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbPreprocessorUnrecognizedHeavyAtom::PdbPreprocessorUnrecognizedHeavyAtom() {}

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
string PdbPreprocessorUnrecognizedHeavyAtom::GetAtomName()
{
    return atom_name_;
}
string PdbPreprocessorUnrecognizedHeavyAtom::GetResidueName()
{
    return residue_name_;
}
int PdbPreprocessorUnrecognizedHeavyAtom::GetResidueSequenceNumber()
{
    return residue_sequence_number_;
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
void PdbPreprocessorUnrecognizedHeavyAtom::SetAtomName(const string atom_name)
{
    atom_name_ = atom_name;
}
void PdbPreprocessorUnrecognizedHeavyAtom::SetResidueName(const string residue_name)
{
    residue_name_ = residue_name;
}
void PdbPreprocessorUnrecognizedHeavyAtom::SetResidueSequenceNumber(int residue_sequence_number)
{
    residue_sequence_number_ = residue_sequence_number;
}

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbPreprocessorUnrecognizedHeavyAtom::Print(ostream &out)
{
}








