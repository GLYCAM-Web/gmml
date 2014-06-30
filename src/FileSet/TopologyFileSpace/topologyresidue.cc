
#include "../../../includes/FileSet/TopologyFileSpace/topologyresidue.hpp"
#include "../../../includes/FileSet/TopologyFileSpace/topologyatom.hpp"
#include "../../../includes/FileSet/TopologyFileSpace/topologybond.hpp"
#include "../../../includes/FileSet/TopologyFileSpace/topologyangle.hpp"
#include "../../../includes/FileSet/TopologyFileSpace/topologydihedral.hpp"

using namespace std;
using namespace TopologyFileSpace;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
TopologyResidue::TopologyResidue() {}

TopologyResidue::TopologyResidue(string residue_name, TopologyAtomMap atoms, int index, int starting_atom_index) :
    residue_name_(residue_name), index_(index), starting_atom_index_(starting_atom_index)
{
    atoms_.clear();
    for(TopologyAtomMap::iterator it = atoms.begin(); it != atoms.end(); it++)
    {
        string atom_name = (*it).first;
        TopologyAtom* atom = (*it).second;
        atoms_[atom_name] = atom;
    }
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
string TopologyResidue::GetResidueName()
{
    return residue_name_;
}
TopologyResidue::TopologyAtomMap TopologyResidue::GetAtoms()
{
    return atoms_;
}
TopologyResidue::TopologyBondMap TopologyResidue::GetBonds()
{
    return bonds_;
}
TopologyResidue::TopologyAngleMap TopologyResidue::GetAngles()
{
    return angles_;
}
TopologyResidue::TopologyDihedralMap TopologyResidue::GetDihedrals()
{
    return dihedrals_;
}
int TopologyResidue::GetIndex()
{
    return index_;
}
int TopologyResidue::GetStartingAtomIndex()
{
    return starting_atom_index_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void TopologyResidue::SetResidueName(const string residue_name)
{
    residue_name_ = residue_name;
}
void TopologyResidue::SetIndex(int index)
{
    index_ = index;
}
void TopologyResidue::SetStartingAtomIndex(int starting_atom_index)
{
    starting_atom_index_ = starting_atom_index;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void TopologyResidue::Print(ostream &out)
{
}





