#include <iostream>

#include "../../../includes/ParameterSet/OffFileSpace/offfileresidue.hpp"
#include "../../../includes/ParameterSet/OffFileSpace/offfileatom.hpp"
#include "../../../includes/common.hpp"

using OffFileSpace::OffFileResidue;

//////////////////////////////////////////////////////////
//                       Constructor                    //
//////////////////////////////////////////////////////////
OffFileResidue::OffFileResidue() :  name_(""), atoms_(), box_angle_(0.0),
                                            box_length_(0.0), box_width_(0.0), box_height_(0.0),
                                            head_atom_index_(-1), tail_atom_index_(-1) {}

///////////////////////////////////////////////////////////
//                           ACCESSOR                   //
//////////////////////////////////////////////////////////
/// Return the residue name
std::string OffFileResidue::GetName()
{
    return name_;
}

OffFileResidue::OffFileAtomVector OffFileResidue::GetAtoms()
{
    return atoms_;
}

/// Return an atom belonging to the residue by a given index number
OffFileSpace::OffFileAtom* OffFileResidue::GetAtomByIndex(int index)
{
    return atoms_[index];
}

/*
>>>>>>> 59b93388cbe59d7a2ba06f22b8a076e9d6382f44
/// Return an atom belonging to the residue by a given order number
OffFileSpace::OffFileAtom* OffFileResidue::GetAtomByOrder(int order)
{
    for(OffFileResidue::AtomMap::iterator it = atoms_.begin(); it != atoms_.end(); it++)
    {
        if(it->second->GetAtomOrder() == order)
            return it->second;
    }
    return NULL;
}
<<<<<<< HEAD
=======
*/

/// Return the box angle of the residue
double OffFileResidue::GetBoxAngle()
{
    return box_angle_;
}

/// Return the box length of the residue
double OffFileResidue::GetBoxLength()
{
    return box_length_;
}

/// Return the box width of the residue
double OffFileResidue::GetBoxWidth()
{
    return box_width_;
}

/// Return the box hight of the residue
double OffFileResidue::GetBoxHeight()
{
    return box_height_;
}

/// Return the head atom index of the residue
int OffFileResidue::GetHeadAtomIndex()
{
    return head_atom_index_;
}

/// Return the tail atom index of the residue
int OffFileResidue::GetTailAtomIndex()
{
    return tail_atom_index_;
}

/*
>>>>>>> 59b93388cbe59d7a2ba06f22b8a076e9d6382f44
OffFileSpace::OffFileAtom* OffFileResidue::GetOffAtomByAtomName(std::string atom_name)
{
    for(OffFileResidue::AtomMap::iterator it = atoms_.begin(); it != atoms_.end(); it++)
    {
        OffFileSpace::OffFileAtom* atom = (*it).second;
        if(atom->GetName().compare(atom_name) == 0)
            return atom;
    }
    return NULL;
}
<<<<<<< HEAD
=======
*/
int OffFileResidue::GetListingIndex()
{
    return listing_index_;
}

//////////////////////////////////////////////////////////
//                           MUTATOR                    //
//////////////////////////////////////////////////////////
/// Set the residue name
void OffFileResidue::SetName(std::string name)
{
    name_ = name;
}

/*
>>>>>>> 59b93388cbe59d7a2ba06f22b8a076e9d6382f44
/// Set the atom map belonging to the residue mapped to their insertion order
void OffFileResidue::SetAtoms(OffFileResidue::AtomMap& atoms)
{
    atoms_.clear();
    for(OffFileResidue::AtomMap::iterator it = atoms.begin(); it != atoms.end(); it++)
    {
        atoms_[it->second->GetAtomOrder()] = it->second;
    }
}
<<<<<<< HEAD

/// Add a new atom to map of the residue
void OffFileResidue::AddAtom(OffFileSpace::OffFileAtom* atom)
{
    atoms_[atom->GetAtomOrder()] = atom;
=======
*/
/// Add a new atom to map of the residue
void OffFileResidue::AddAtom(OffFileSpace::OffFileAtom* atom)
{
    atoms_.push_back(atom);
}

/// Set the box angle of the residue
void OffFileResidue::SetBoxAngle(double box_angle)
{
    box_angle_ = box_angle;
}

/// Set the box length of the residue
void OffFileResidue::SetBoxLength(double box_length)
{
    box_length_ = box_length;
}

/// Set the box width of the residue
void OffFileResidue::SetBoxWidth(double box_width)
{
    box_width_ = box_width;
}

/// Set the box height of the residue
void OffFileResidue::SetBoxHeight(double box_height)
{
    box_height_ = box_height;
}

/// Set the head atom index of the residue
void OffFileResidue::SetHeadAtomIndex(int head_atom_index)
{
    head_atom_index_ = head_atom_index;
}

/// Set the tail atom index of the residue
void OffFileResidue::SetTailAtomIndex(int tail_atom_index)
{
    tail_atom_index_ = tail_atom_index;
}
void OffFileResidue::SetListingIndex(int listing_index)
{
    listing_index_ = listing_index;
}

//////////////////////////////////////////////////////////
//                     DISPLAY FUNCTIONS                //
//////////////////////////////////////////////////////////
void OffFileResidue::Print(std::ostream& out)
{}
