
#include "../../../includes/InputSet/TopologyFileSpace/topologyresidue.hpp"
#include "../../../includes/InputSet/TopologyFileSpace/topologyatom.hpp"
#include "../../../includes/common.hpp"

using namespace std;
using namespace gmml;
using namespace TopologyFileSpace;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
TopologyResidue::TopologyResidue()
{
    residue_name_ = "";
    atoms_ = TopologyAtomVector();
    index_ = iNotSet;
    starting_atom_index_ = iNotSet;
    is_residue_solvent_=false;
}

TopologyResidue::TopologyResidue(string residue_name, TopologyAtomVector atoms, int index, int starting_atom_index) :
    residue_name_(residue_name), index_(index), starting_atom_index_(starting_atom_index)
{
    atoms_.clear();
    for(TopologyAtomVector::iterator it = atoms.begin(); it != atoms.end(); it++)
    {
        TopologyAtom* atom = (*it);
        atoms_.push_back(atom);
    }
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
string TopologyResidue::GetResidueName()
{
    return residue_name_;
}
TopologyResidue::TopologyAtomVector TopologyResidue::GetAtoms()
{
    return atoms_;
}
TopologyAtom* TopologyResidue::GetAtomByIndex(int index)
{
    for(TopologyAtomVector::iterator it = atoms_.begin(); it != atoms_.end(); it++)
    {
        TopologyAtom* atom = (*it);
        if(atom->GetIndex() == index)
            return atom;
    }
    return NULL;
}
int TopologyResidue::GetIndex()
{
    return index_;
}
int TopologyResidue::GetStartingAtomIndex()
{
    return starting_atom_index_;
}
bool TopologyResidue::GetIsResidueSolvent()
{
    return is_residue_solvent_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void TopologyResidue::SetResidueName(const string residue_name)
{
    residue_name_ = residue_name;
}
void TopologyResidue::AddAtom(TopologyAtom *atom)
{
    atoms_.push_back(atom);
}
void TopologyResidue::SetIndex(int index)
{
    index_ = index;
}
void TopologyResidue::SetStartingAtomIndex(int starting_atom_index)
{
    starting_atom_index_ = starting_atom_index;
}
void TopologyResidue::SetIsResidueSolvent(bool is_residue_solvent)
{
    is_residue_solvent_ = is_residue_solvent;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void TopologyResidue::Print(ostream &out)
{
    out << "Residue Name: " << residue_name_ << ", Residue Index: " << index_ << ", Starting Atom Index: " << starting_atom_index_ << endl << endl;
//    out << "------------------------------------- " << residue_name_ << " --------------------------------------" << endl;
    for(TopologyAtomVector::iterator it = atoms_.begin(); it != atoms_.end(); it++)
    {
//        string atom_name = (*it).first;
        TopologyAtom* atom = (*it);
        atom->Print(out);
    }
}





