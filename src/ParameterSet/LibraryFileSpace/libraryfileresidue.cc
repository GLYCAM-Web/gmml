#include <iostream>

#include "../../../includes/ParameterSet/LibraryFileSpace/libraryfileresidue.hpp"
#include "../../../includes/ParameterSet/LibraryFileSpace/libraryfileatom.hpp"
#include "../../../includes/common.hpp"

using LibraryFileSpace::LibraryFileResidue;

//////////////////////////////////////////////////////////
//                       Constructor                    //
//////////////////////////////////////////////////////////
LibraryFileResidue::LibraryFileResidue() :  name_(""), atoms_(), box_angle_(0.0),
                                            box_length_(0.0), box_width_(0.0), box_height_(0.0),
                                            head_atom_index_(-1), tail_atom_index_(-1) {}

LibraryFileResidue::LibraryFileResidue(std::string &name, int listing_index) :
    name_(name), box_angle_(0.0), box_length_(0.0), box_width_(0.0), box_height_(0.0),
    head_atom_index_(-1), tail_atom_index_(-1), listing_index_(listing_index)
{
    atoms_ = AtomMap();
}

LibraryFileResidue::LibraryFileResidue( std::string& name, int listing_index,
                                        std::vector<LibraryFileSpace::LibraryFileAtom*>& atoms,
                                        int head_atom_index, int tail_atom_index,
                                        double box_angle, double box_length,
                                        double box_width, double box_height) :
    name_(name), box_angle_(box_angle), box_length_(box_length), box_width_(box_width),
    box_height_(box_height), head_atom_index_(head_atom_index), tail_atom_index_(tail_atom_index),
    listing_index_(listing_index)
{
    atoms_.clear();
    for(std::vector<LibraryFileSpace::LibraryFileAtom*>::iterator it = atoms.begin(); it != atoms.end(); it++)
    {
        atoms_[(*it)->GetAtomOrder()] = *it;
    }
}

///////////////////////////////////////////////////////////
//                           ACCESSOR                   //
//////////////////////////////////////////////////////////
/// Return the residue name
std::string LibraryFileResidue::GetName()
{
    return name_;
}

/// Return a map of all atoms belonging to the residue mapped to their order of insertion into the residue
LibraryFileResidue::AtomMap LibraryFileResidue::GetAtoms()
{
    return atoms_;
}

/// Return an atom belonging to the residue by a given index number
LibraryFileSpace::LibraryFileAtom* LibraryFileResidue::GetAtomByIndex(int index)
{
    return atoms_[index];
}

/// Return an atom belonging to the residue by a given order number
LibraryFileSpace::LibraryFileAtom* LibraryFileResidue::GetAtomByOrder(int order)
{
    for(LibraryFileResidue::AtomMap::iterator it = atoms_.begin(); it != atoms_.end(); it++)
    {
        if(it->second->GetAtomOrder() == order)
            return it->second;
    }
    return NULL;
}

/// Return the box angle of the residue
double LibraryFileResidue::GetBoxAngle()
{
    return box_angle_;
}

/// Return the box length of the residue
double LibraryFileResidue::GetBoxLength()
{
    return box_length_;
}

/// Return the box width of the residue
double LibraryFileResidue::GetBoxWidth()
{
    return box_width_;
}

/// Return the box hight of the residue
double LibraryFileResidue::GetBoxHeight()
{
    return box_height_;
}

/// Return the head atom index of the residue
int LibraryFileResidue::GetHeadAtomIndex()
{
    return head_atom_index_;
}

/// Return the tail atom index of the residue
int LibraryFileResidue::GetTailAtomIndex()
{
    return tail_atom_index_;
}
LibraryFileSpace::LibraryFileAtom* LibraryFileResidue::GetLibraryAtomByAtomName(std::string atom_name)
{
    for(LibraryFileResidue::AtomMap::iterator it = atoms_.begin(); it != atoms_.end(); it++)
    {
        LibraryFileSpace::LibraryFileAtom* atom = (*it).second;
        if(atom->GetName().compare(atom_name) == 0)
            return atom;
    }
    return NULL;
}
int LibraryFileResidue::GetListingIndex()
{
    return listing_index_;
}

//////////////////////////////////////////////////////////
//                           MUTATOR                    //
//////////////////////////////////////////////////////////
/// Set the residue name
void LibraryFileResidue::SetName(std::string name)
{
    name_ = name;
}

/// Set the atom map belonging to the residue mapped to their insertion order
void LibraryFileResidue::SetAtoms(LibraryFileResidue::AtomMap& atoms)
{
    atoms_.clear();
    for(LibraryFileResidue::AtomMap::iterator it = atoms.begin(); it != atoms.end(); it++)
    {
        atoms_[it->second->GetAtomOrder()] = it->second;
    }
}

/// Add a new atom to map of the residue
void LibraryFileResidue::AddAtom(LibraryFileSpace::LibraryFileAtom* atom)
{
    atoms_[atom->GetAtomOrder()] = atom;
}

/// Set the box angle of the residue
void LibraryFileResidue::SetBoxAngle(double box_angle)
{
    box_angle_ = box_angle;
}

/// Set the box length of the residue
void LibraryFileResidue::SetBoxLength(double box_length)
{
    box_length_ = box_length;
}

/// Set the box width of the residue
void LibraryFileResidue::SetBoxWidth(double box_width)
{
    box_width_ = box_width;
}

/// Set the box height of the residue
void LibraryFileResidue::SetBoxHeight(double box_height)
{
    box_height_ = box_height;
}

/// Set the head atom index of the residue
void LibraryFileResidue::SetHeadAtomIndex(int head_atom_index)
{
    head_atom_index_ = head_atom_index;
}

/// Set the tail atom index of the residue
void LibraryFileResidue::SetTailAtomIndex(int tail_atom_index)
{
    tail_atom_index_ = tail_atom_index;
}
void LibraryFileResidue::SetListingIndex(int listing_index)
{
    listing_index_ = listing_index;
}

//////////////////////////////////////////////////////////
//                     DISPLAY FUNCTIONS                //
//////////////////////////////////////////////////////////
void LibraryFileResidue::Print(std::ostream& out)
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

    for(LibraryFileResidue::AtomMap::iterator it = atoms_.begin(); it != atoms_.end(); it++)
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
