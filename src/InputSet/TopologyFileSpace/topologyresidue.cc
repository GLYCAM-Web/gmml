
#include "../../../includes/InputSet/TopologyFileSpace/topologyresidue.hpp"
#include "../../../includes/InputSet/TopologyFileSpace/topologyatom.hpp"
#include "../../../includes/common.hpp"

using TopologyFileSpace::TopologyResidue;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
TopologyResidue::TopologyResidue()
{
    residue_name_ = "";
    atoms_ = TopologyAtomVector();
    index_ = gmml::iNotSet;
    starting_atom_index_ = gmml::iNotSet;
    is_residue_solvent_=false;
}

TopologyResidue::TopologyResidue(std::string residue_name, TopologyAtomVector atoms, int index, int starting_atom_index) :
    residue_name_(residue_name), index_(index), starting_atom_index_(starting_atom_index)
{
    atoms_.clear();
    for(TopologyAtomVector::iterator it = atoms.begin(); it != atoms.end(); it++)
    {
        TopologyFileSpace::TopologyAtom* atom = (*it);
        atoms_.push_back(atom);
    }
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
std::string TopologyResidue::GetResidueName()
{
    return residue_name_;
}
TopologyResidue::TopologyAtomVector TopologyResidue::GetAtoms()
{
    return atoms_;
}
TopologyFileSpace::TopologyAtom* TopologyResidue::GetAtomByIndex(int index)
{
    for(TopologyAtomVector::iterator it = atoms_.begin(); it != atoms_.end(); it++)
    {
        TopologyFileSpace::TopologyAtom* atom = (*it);
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
void TopologyResidue::SetResidueName(const std::string residue_name)
{
    residue_name_ = residue_name;
}
void TopologyResidue::AddAtom(TopologyFileSpace::TopologyAtom *atom)
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
void TopologyResidue::Print(std::ostream &out)
{
    out << "Residue Name: " << residue_name_ << ", Residue Index: " << index_ << ", Starting Atom Index: " << starting_atom_index_ << std::endl << std::endl;
//    out << "------------------------------------- " << residue_name_ << " --------------------------------------" << std::endl;
    for(TopologyAtomVector::iterator it = atoms_.begin(); it != atoms_.end(); it++)
    {
//        std::string atom_name = (*it).first;
        TopologyFileSpace::TopologyAtom* atom = (*it);
        atom->Print(out);
    }
}
