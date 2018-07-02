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

OffFileResidue::OffFileResidue(std::string &name, int listing_index) :
    name_(name), box_angle_(0.0), box_length_(0.0), box_width_(0.0), box_height_(0.0),
    head_atom_index_(-1), tail_atom_index_(-1), listing_index_(listing_index)
{
    atoms_ = AtomMap();
}

OffFileResidue::OffFileResidue( std::string& name, int listing_index,
                                        std::vector<OffFileSpace::OffFileAtom*>& atoms,
                                        int head_atom_index, int tail_atom_index,
                                        double box_angle, double box_length,
                                        double box_width, double box_height) :
    name_(name), box_angle_(box_angle), box_length_(box_length), box_width_(box_width),
    box_height_(box_height), head_atom_index_(head_atom_index), tail_atom_index_(tail_atom_index),
    listing_index_(listing_index)
{
    atoms_.clear();
    for(std::vector<OffFileSpace::OffFileAtom*>::iterator it = atoms.begin(); it != atoms.end(); it++)
    {
        atoms_[(*it)->GetAtomOrder()] = *it;
    }
}

///////////////////////////////////////////////////////////
//                           ACCESSOR                   //
//////////////////////////////////////////////////////////
/// Return the residue name
std::string OffFileResidue::GetName()
{
    return name_;
}

/// Return a map of all atoms belonging to the residue mapped to their order of insertion into the residue
OffFileResidue::AtomMap OffFileResidue::GetAtoms()
{
    return atoms_;
}

/// Return an atom belonging to the residue by a given index number
OffFileSpace::OffFileAtom* OffFileResidue::GetAtomByIndex(int index)
{
    return atoms_[index];
}

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

/// Set the atom map belonging to the residue mapped to their insertion order
void OffFileResidue::SetAtoms(OffFileResidue::AtomMap& atoms)
{
    atoms_.clear();
    for(OffFileResidue::AtomMap::iterator it = atoms.begin(); it != atoms.end(); it++)
    {
        atoms_[it->second->GetAtomOrder()] = it->second;
    }
}

/// Add a new atom to map of the residue
void OffFileResidue::AddAtom(OffFileSpace::OffFileAtom* atom)
{
    atoms_[atom->GetAtomOrder()] = atom;
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
{
    out << std::setw(80) << "***************************** " + name_ + " *****************************" << std::endl;
    out << std::setw(60) << "======================== ATOMS ========================" << std::endl;
    out << std::setw(2) << "#"
        << std::setw(6) << "NAME"
        << std::setw(6) << "TYPE"
        << std::setw(6) << "RES_I"
        << std::setw(6) << "ATM_I"
        << std::setw(10) << "CHG"
        << std::setw(45) << "CRD"
        << std::setw(20) << "BND_ATMS"
        << std::endl;

    for(OffFileResidue::AtomMap::iterator it = atoms_.begin(); it != atoms_.end(); it++)
    {
        it->second->Print(out);
        for(unsigned int j = 0; j < it->second->GetBondedAtomsIndices().size(); j++)
            out << std::setw(4) << it -> second -> GetBondedAtomsIndices()[j] << ", ";
        out << std::endl;
    }

    out << std::endl << std::setw(60) << "======================= Bound Box =======================" << std::endl;
    out << std::setw(6) << "ANGL"
        << std::setw(6) << "LEN"
        << std::setw(6) << "WDTH"
        << std::setw(6) << "HGHT"
        << std::endl;

    if(box_angle_ == gmml::dNotSet)
        out << std::setw(6) << "--";
    else
        out << std::setw(6) << box_angle_;
    if(box_length_ == gmml::dNotSet)
        out << std::setw(6) << "--";
    else
        out << std::setw(6) << box_length_;
    if(box_width_ == gmml::dNotSet)
        out << std::setw(6) << "--";
    else
        out << std::setw(6) << box_width_;
    if(box_height_ == gmml::dNotSet)
        out << std::setw(6) << "--";
    else
        out << std::setw(6) << box_height_;
    out << std::endl;

    out << "HEAD ATOM INDEX: " << head_atom_index_ << std::endl;
    out << "TAIL ATOM INDEX: " << tail_atom_index_ << std::endl;

}
