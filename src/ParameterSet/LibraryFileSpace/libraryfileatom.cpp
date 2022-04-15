#include "../../../includes/ParameterSet/LibraryFileSpace/libraryfileatom.hpp"
#include "../../../includes/GeometryTopology/coordinate.hpp"
#include <iostream>

using LibraryFileSpace::LibraryFileAtom;

//////////////////////////////////////////////////////////
//                       Constructor                    //
//////////////////////////////////////////////////////////
LibraryFileAtom::LibraryFileAtom() : type_(""), name_(""), residue_index_(-1), atom_index_(-1), atomic_number_(0), charge_(0.0), coordinate_(),
    atom_order_(0) {}

LibraryFileAtom::LibraryFileAtom(std::string type, std::string name, int residue_index, int atom_index, int atomic_number, double charge) :
    type_(type), name_(name), residue_index_(residue_index), atom_index_(atom_index), atomic_number_(atomic_number), charge_(charge), coordinate_(),
    atom_order_(0) {}

LibraryFileAtom::LibraryFileAtom(std::string type, std::string name, int residue_index, int atom_index, int atomic_number, double charge,
                                 GeometryTopology::Coordinate coordinate, std::vector<int> bonded_atoms_indices, int atom_order) :
    type_(type), name_(name), residue_index_(residue_index), atom_index_(atom_index), atomic_number_(atomic_number), charge_(charge),
    coordinate_(coordinate), atom_order_(atom_order)
{
    bonded_atoms_indices_.clear();
    for(std::vector<int>::const_iterator it = bonded_atoms_indices.begin(); it != bonded_atoms_indices.end(); it++)
    {
        bonded_atoms_indices_.push_back(*it);
    }
}

//////////////////////////////////////////////////////////
//                           ACCESSOR                   //
//////////////////////////////////////////////////////////
/// Return atom type
std::string LibraryFileAtom::GetType()
{
    return type_;
}

/// Return atom name
std::string LibraryFileAtom::GetName()
{
    return name_;
}

/// Return the residue index in the list of residues in the library file that the current atom belongs to it
int LibraryFileAtom::GetResidueIndex()
{
    return residue_index_;
}

/// return atom index in the belonging residue
int LibraryFileAtom::GetAtomIndex()
{
    return atom_index_;
}

/// Return atomic number of the atom
int LibraryFileAtom::GetAtomicNumber()
{
    return atomic_number_;
}

/// Return charge of the atom
double LibraryFileAtom::GetCharge()
{
    return charge_;
}

/// Return position of the atom
GeometryTopology::Coordinate LibraryFileAtom::GetCoordinate()
{
    return coordinate_;
}

/// Return all bonded atom indices of the atom
std::vector<int> LibraryFileAtom::GetBondedAtomsIndices()
{
    return bonded_atoms_indices_;
}

/// Return inserting order of the atom in a residue
int LibraryFileAtom::GetAtomOrder()
{
    return atom_order_;
}

//////////////////////////////////////////////////////////
//                           MUTATOR                    //
//////////////////////////////////////////////////////////
/// Set atom type
void LibraryFileAtom::SetType(const std::string type)
{
    type_ = type;
}

/// Set atom name
void LibraryFileAtom::SetName(const std::string name)
{
    name_ = name;
}

/// Set the index of the residue in the list of residues in the library file that the current atom belongs to it
void LibraryFileAtom::SetResidueIndex(int residue_index)
{
    residue_index_ = residue_index;
}

/// Set the atom index in its residue
void LibraryFileAtom::SetAtomIndex(int atom_index)
{
    atom_index_ = atom_index;
}

/// Set the atomic number of the atom
void LibraryFileAtom::SetAtomicNumber(int atomic_number)
{
    atomic_number_ = atomic_number;
}

/// Set the postion of the atom
void LibraryFileAtom::SetCoordinate(GeometryTopology::Coordinate& coordinate)
{
    coordinate_.SetX(coordinate.GetX());
    coordinate_.SetY(coordinate.GetY());
    coordinate_.SetZ(coordinate.GetZ());
}

/// Set all the bonded atom indices to the current atom
void LibraryFileAtom::SetBondedAtomsIndices(const std::vector<int> bonded_atoms_indices)
{
    bonded_atoms_indices_.clear();
    for(std::vector<int>::const_iterator it = bonded_atoms_indices.begin(); it != bonded_atoms_indices.end(); it++)
    {
        bonded_atoms_indices_.push_back(*it);
    }
}

/// Add a new bonded atom index to the list of bonded atom indices
void LibraryFileAtom::AddBondedAtomIndex(int index)
{
    bonded_atoms_indices_.push_back(index);
}

/// Set the order of insertion of the current atom into its residue
void LibraryFileAtom::SetAtomOrder(int atom_order)
{
    atom_order_ = atom_order;
}

//////////////////////////////////////////////////////////
//                     DISPLAY FUNCTIONS                //
//////////////////////////////////////////////////////////
void LibraryFileAtom::Print(std::ostream& out)
{
    out << std::setw(2) << atom_order_
        << std::setw(6) << name_
        << std::setw(6) << type_
        << std::setw(6) << residue_index_
        << std::setw(6) << atom_index_
        << std::setw(10) << charge_
        << std::setw(15) << coordinate_.GetX() << ", " << std::setw(15) << coordinate_.GetY() << ", " << std::setw(15) << coordinate_.GetZ();
}
