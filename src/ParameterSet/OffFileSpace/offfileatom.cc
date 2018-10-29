#include "../../../includes/ParameterSet/OffFileSpace/offfileatom.hpp"
#include "../../../includes/GeometryTopology/coordinate.hpp"
#include <iostream>

using OffFileSpace::OffFileAtom;

//////////////////////////////////////////////////////////
//                       Constructor                    //
//////////////////////////////////////////////////////////
OffFileAtom::OffFileAtom() : type_(""), name_(""), residue_index_(-1), atom_index_(-1), atomic_number_(0), charge_(0.0), coordinate_(),
    atom_order_(0) {}

//////////////////////////////////////////////////////////
//                           ACCESSOR                   //
//////////////////////////////////////////////////////////

/// Return atom name
std::string OffFileAtom::GetName()
{
    return name_;
}

/// Return atom type
std::string OffFileAtom::GetType()
{
    return type_;
}

/// Return the residue index in the list of residues in the off file that the current atom belongs to it
int OffFileAtom::GetResidueIndex()
{
    return residue_index_;
}

/// return atom index in the belonging residue
int OffFileAtom::GetAtomIndex()
{
    return atom_index_;
}

/// Return atomic number of the atom
int OffFileAtom::GetAtomicNumber()
{
    return atomic_number_;
}

/// Return charge of the atom
double OffFileAtom::GetCharge()
{
    return charge_;
}

/// Return position of the atom
GeometryTopology::Coordinate OffFileAtom::GetCoordinate()
{
    return coordinate_;
}

/// Return all bonded atom indices of the atom
std::vector<int> OffFileAtom::GetBondedAtomsIndices()
{
    return bonded_atoms_indices_;
}

/// Return inserting order of the atom in a residue
int OffFileAtom::GetAtomOrder()
{
    return atom_order_;
}

//////////////////////////////////////////////////////////
//                           MUTATOR                    //
//////////////////////////////////////////////////////////
/// Set atom name
void OffFileAtom::SetName(const std::string name)
{
    name_ = name;
}

/// Set atom type
void OffFileAtom::SetType(const std::string type)
{
    type_ = type;
}

/// Set the index of the residue in the list of residues in the off file that the current atom belongs to it
void OffFileAtom::SetResidueIndex(int residue_index)
{
    residue_index_ = residue_index;
}

/// Set the atom index in its residue
void OffFileAtom::SetAtomIndex(int atom_index)
{
    atom_index_ = atom_index;
}

/// Set the atomic number of the atom
void OffFileAtom::SetAtomicNumber(int atomic_number)
{
    atomic_number_ = atomic_number;
}

/// Set the postion of the atom
void OffFileAtom::SetCoordinate(GeometryTopology::Coordinate* coordinate)
{
    coordinate_.SetX(coordinate->GetX());
    coordinate_.SetY(coordinate->GetY());
    coordinate_.SetZ(coordinate->GetZ());
}

/*
>>>>>>> 59b93388cbe59d7a2ba06f22b8a076e9d6382f44
/// Set all the bonded atom indices to the current atom
void OffFileAtom::SetBondedAtomsIndices(const std::vector<int> bonded_atoms_indices)
{
    bonded_atoms_indices_.clear();
    for(std::vector<int>::const_iterator it = bonded_atoms_indices.begin(); it != bonded_atoms_indices.end(); it++)
    {
        bonded_atoms_indices_.push_back(*it);
    }
}
<<<<<<< HEAD

=======
*/
/// Add a new bonded atom index to the list of bonded atom indices
void OffFileAtom::AddBondedAtomIndex(int index)
{
    bonded_atoms_indices_.push_back(index);
}

/// Set the order of insertion of the current atom into its residue
void OffFileAtom::SetAtomOrder(int atom_order)
{
    atom_order_ = atom_order;
}

void OffFileAtom::SetAtomCharge(double atom_charge)
{
    charge_=atom_charge;
}
//////////////////////////////////////////////////////////
//                     DISPLAY FUNCTIONS                //
//////////////////////////////////////////////////////////
void OffFileAtom::Print(std::ostream& out)
{
    out << std::setw(2) << name_
	<< std::setw(6) << atom_order_
        << std::setw(6) << type_
        << std::setw(6) << residue_index_
        << std::setw(6) << atom_index_
        << std::setw(10) << charge_
        << std::setw(15) << coordinate_.GetX() << ", " << std::setw(15) << coordinate_.GetY() << ", " << std::setw(15) << coordinate_.GetZ();
}
