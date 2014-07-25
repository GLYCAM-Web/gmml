#include <iostream>

#include "../../../includes/ParameterSet/LibraryFileSpace/libraryfileresidue.hpp"
#include "../../../includes/ParameterSet/LibraryFileSpace/libraryfileatom.hpp"
#include "../../../includes/common.hpp"

using namespace std;
using namespace LibraryFileSpace;
using namespace gmml;

//////////////////////////////////////////////////////////
//                       Constructor                    //
//////////////////////////////////////////////////////////
LibraryFileResidue::LibraryFileResidue() : name_(""), atoms_(), box_angle_(0.0), box_length_(0.0), box_width_(0.0), box_height_(0.0),
    head_atom_index_(-1), tail_atom_index_(-1) {}

LibraryFileResidue::LibraryFileResidue(string &name) : name_(name), atoms_(), box_angle_(0.0), box_length_(0.0), box_width_(0.0), box_height_(0.0),
    head_atom_index_(-1), tail_atom_index_(-1) {}

LibraryFileResidue::LibraryFileResidue(string &name, vector<LibraryFileAtom*> &atoms, int head_atom_index, int tail_atom_index, double box_angle, double box_length,
                                       double box_width, double box_height) :
    name_(name), box_angle_(box_angle), box_length_(box_length), box_width_(box_width), box_height_(box_height), head_atom_index_(head_atom_index),
    tail_atom_index_(tail_atom_index)
{
    atoms_.clear();
    for(vector<LibraryFileAtom*>::iterator it = atoms.begin(); it != atoms.end(); it++)
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
LibraryFileAtom* LibraryFileResidue::GetAtomByIndex(int index)
{
    return atoms_[index];
}

/// Return an atom belonging to the residue by a given order number
LibraryFileAtom* LibraryFileResidue::GetAtomByOrder(int order)
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
LibraryFileAtom* LibraryFileResidue::GetLibraryAtomByAtomName(string atom_name)
{
    LibraryFileAtom* library_file_atom = new LibraryFileAtom();
    for(LibraryFileResidue::AtomMap::iterator it = atoms_.begin(); it != atoms_.end(); it++)
     if(it->second->GetName() == atom_name)
        library_file_atom = it->second;
    else
        library_file_atom = NULL;
    return library_file_atom;
}

//////////////////////////////////////////////////////////
//                           MUTATOR                    //
//////////////////////////////////////////////////////////
/// Set the residue name
void LibraryFileResidue::SetName(std::string& name)
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
void LibraryFileResidue::AddAtom(LibraryFileAtom* atom)
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

//////////////////////////////////////////////////////////
//                     DISPLAY FUNCTIONS                //
//////////////////////////////////////////////////////////
void LibraryFileResidue::Print(std::ostream& out)
{
    out << setw(80) << "***************************** " + name_ + " *****************************" << endl;
    out << setw(60) << "======================== ATOMS ========================" << endl;
    out << setw(2) << "#"
        << setw(6) << "NAME"
        << setw(6) << "TYPE"
        << setw(6) << "RES_I"
        << setw(6) << "ATM_I"
        << setw(10) << "CHG"
        << setw(45) << "CRD"
        << setw(20) << "BND_ATMS"
        << endl;

    for(LibraryFileResidue::AtomMap::iterator it = atoms_.begin(); it != atoms_.end(); it++)
    {
        it->second->Print(out);
        for(unsigned int j = 0; j < it->second->GetBondedAtomsIndices().size(); j++)
            out << setw(4) << it -> second -> GetBondedAtomsIndices()[j] << ", ";
        out << endl;
    }

    out << endl << setw(60) << "======================= Bound Box =======================" << endl;
    out << setw(6) << "ANGL"
        << setw(6) << "LEN"
        << setw(6) << "WDTH"
        << setw(6) << "HGHT"
        << endl;

    if(box_angle_ == dNotSet)
        out << setw(6) << "--";
    else
        out << setw(6) << box_angle_;
    if(box_length_ == dNotSet)
        out << setw(6) << "--";
    else
        out << setw(6) << box_length_;
    if(box_width_ == dNotSet)
        out << setw(6) << "--";
    else
        out << setw(6) << box_width_;
    if(box_height_ == dNotSet)
        out << setw(6) << "--";
    else
        out << setw(6) << box_height_;
    out << endl;

    out << "HEAD ATOM INDEX: " << head_atom_index_ << endl;
    out << "TAIL ATOM INDEX: " << tail_atom_index_ << endl;

}
