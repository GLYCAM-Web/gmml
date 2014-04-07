
#include "../../../includes/Resolver/PdbPreprocessor/pdbpreprocessorreplacedhydrogen.hpp"

using namespace std;
using namespace PdbPreprocessorSpace;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbPreprocessorReplacedHydrogen::PdbPreprocessorReplacedHydrogen() {}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
int PdbPreprocessorReplacedHydrogen::GetAtomIndex()
{
    return atom_index_;
}
string PdbPreprocessorReplacedHydrogen::GetAtomName()
{
    return atom_name_;
}
string PdbPreprocessorReplacedHydrogen::GetResidueName()
{
    return residue_name_;
}
int PdbPreprocessorReplacedHydrogen::GetResidueNumber()
{
    return residue_number_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbPreprocessorReplacedHydrogen::SetAtomIndex(int atom_index)
{
    atom_index_ = atom_index;
}
void PdbPreprocessorReplacedHydrogen::SetAtomName(const string atom_name)
{
    atom_name_ = atom_name;
}
void PdbPreprocessorReplacedHydrogen::SetResidueName(const string residue_name)
{
    residue_name_ = residue_name;
}
void PdbPreprocessorReplacedHydrogen::SetResidueNumber(int residue_number)
{
    residue_number_ = residue_number;
}

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbPreprocessorReplacedHydrogen::Print(ostream &out)
{
}









