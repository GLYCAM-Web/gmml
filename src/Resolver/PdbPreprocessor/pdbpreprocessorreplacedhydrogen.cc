
#include "../../../includes/Resolver/PdbPreprocessor/pdbpreprocessorreplacedhydrogen.hpp"

using namespace std;
using namespace PdbPreprocessorSpace;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbPreprocessorReplacedHydrogen::PdbPreprocessorReplacedHydrogen() {}

PdbPreprocessorReplacedHydrogen::PdbPreprocessorReplacedHydrogen(int atom_serial_number, string atom_name, string residue_name, int residue_sequence_number) :
    atom_serial_number_(atom_serial_number), atom_name_(atom_name), residue_name_(residue_name), residue_sequence_number_(residue_sequence_number) {}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
int PdbPreprocessorReplacedHydrogen::GetAtomSerialNumber()
{
    return atom_serial_number_;
}
string PdbPreprocessorReplacedHydrogen::GetAtomName()
{
    return atom_name_;
}
string PdbPreprocessorReplacedHydrogen::GetResidueName()
{
    return residue_name_;
}
int PdbPreprocessorReplacedHydrogen::GetResidueSequenceNumber()
{
    return residue_sequence_number_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbPreprocessorReplacedHydrogen::SetAtomSerialNumber(int atom_serial_number)
{
    atom_serial_number_ = atom_serial_number;
}
void PdbPreprocessorReplacedHydrogen::SetAtomName(const string atom_name)
{
    atom_name_ = atom_name;
}
void PdbPreprocessorReplacedHydrogen::SetResidueName(const string residue_name)
{
    residue_name_ = residue_name;
}
void PdbPreprocessorReplacedHydrogen::SetResidueSequenceNumber(int residue_sequence_number)
{
    residue_sequence_number_ = residue_sequence_number;
}

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbPreprocessorReplacedHydrogen::Print(ostream &out)
{
    cout << "Atom name: " << atom_name_
         << ", Serial number: " << atom_serial_number_
         << ", Residue name: " << residue_name_
         << ", Sequence number: " << residue_sequence_number_
         << endl;
}









