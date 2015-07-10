#include "../../../includes/FileSet/PdbqtFileSpace/pdbqtrootcard.hpp"
#include "../../../includes/FileSet/PdbqtFileSpace/pdbqtatomcard.hpp"
#include "../../../includes/utils.hpp"
#include "../../../includes/common.hpp"

using namespace std;
using namespace gmml;
using namespace PdbqtFileSpace;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbqtRootCard::PdbqtRootCard(){}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
PdbqtRootCard::AtomCardVector PdbqtRootCard::GetRootAtoms()
{
    return root_atoms_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbqtRootCard::SetRootAtoms(AtomCardVector root_atoms)
{
    root_atoms_.clear();
    for(AtomCardVector::iterator it = root_atoms.begin(); it != root_atoms.end(); it++)
    {
        root_atoms_.push_back(*it);
    }
}
void PdbqtRootCard::AddRootAtom(PdbqtAtomCard *root_atom)
{
    root_atoms_.push_back(root_atom);
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbqtRootCard::Print(ostream &out)
{
}


