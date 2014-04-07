
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
int PdbPreprocessorUnrecognizedHeavyAtom::GetAtomIndex()
{
    return atom_index_;
}
string PdbPreprocessorUnrecognizedHeavyAtom::GetAtomName()
{
    return atom_name_;
}
string PdbPreprocessorUnrecognizedHeavyAtom::GetResidueName()
{
    return residue_name_;
}
int PdbPreprocessorUnrecognizedHeavyAtom::GetResidueNumber()
{
    return residue_number_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbPreprocessorUnrecognizedHeavyAtom::SetResidueChainId(char residue_chain_id)
{
    residue_chain_id_ = residue_chain_id;
}
void PdbPreprocessorUnrecognizedHeavyAtom::SetAtomIndex(int atom_index)
{
    atom_index_ = atom_index;
}
void PdbPreprocessorUnrecognizedHeavyAtom::SetAtomName(const string atom_name)
{
    atom_name_ = atom_name;
}
void PdbPreprocessorUnrecognizedHeavyAtom::SetResidueName(const string residue_name)
{
    residue_name_ = residue_name;
}
void PdbPreprocessorUnrecognizedHeavyAtom::SetResidueNumber(int residue_number)
{
    residue_number_ = residue_number;
}

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbPreprocessorUnrecognizedHeavyAtom::Print(ostream &out)
{
}








