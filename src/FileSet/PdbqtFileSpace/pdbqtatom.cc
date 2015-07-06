#include "../../../includes/FileSet/PdbqtFileSpace/pdbqtatom.hpp"
#include "../../../includes/utils.hpp"
#include "../../../includes/common.hpp"

using namespace std;
using namespace PdbqtFileSpace;
using namespace gmml;
using namespace Geometry;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbqtAtom::PdbqtAtom() {}

//////////////////////////////////////////////////////////
//                       ACCESSOR                       //
//////////////////////////////////////////////////////////
string PdbqtAtom::GetCardType()
{
    return card_type_;
}

int PdbqtAtom::GetAtomSerialNumber()
{
    return atom_serial_number_;
}

string PdbqtAtom::GetAtomName()
{
    return atom_name_;
}

char PdbqtAtom::GetAtomAlternateLocation()
{
    return atom_alternate_location_;
}

string PdbqtAtom::GetAtomResidueName()
{
    return atom_residue_name_;
}

char PdbqtAtom::GetAtomChainId()
{
    return atom_chain_id_;
}

int PdbqtAtom::GetAtomResidueSequenceNumber()
{
    return atom_residue_sequence_number_;
}

char PdbqtAtom::GetAtomInsertionCode()
{
    return atom_insertion_code_;
}

Geometry::Coordinate PdbqtAtom::GetAtomOrthogonalCoordinate()
{
    return atom_orthogonal_coordinate_;
}

double PdbqtAtom::GetAtomOccupancy()
{
    return atom_occupancy_;
}

double PdbqtAtom::GetAtomTempretureFactor()
{
    return atom_temperature_factor_;
}

double PdbqtAtom::GetAtomCharge()
{
    return atom_charge_;
}
string PdbqtAtom::GetAtomType()
{
    return atom_type_;
}

//////////////////////////////////////////////////////////
//                       MUTATOR                        //
//////////////////////////////////////////////////////////
void PdbqtAtom::SetCardType(const string card_type)
{
    card_type_ = card_type;
}

void PdbqtAtom::SetAtomSerialNumber(int atom_serial_number)
{
    atom_serial_number_ = atom_serial_number;
}

void PdbqtAtom::SetAtomName(const string atom_name)
{
    atom_name_ = atom_name;
}

void PdbqtAtom::SetAtomAlternateLocation(char atom_alternate_location)
{
    atom_alternate_location_ = atom_alternate_location;
}

void PdbqtAtom::SetAtomResidueName(const string atom_residue_name)
{
    atom_residue_name_ = atom_residue_name;
}

void PdbqtAtom::SetAtomChainId(char atom_chain_id)
{
    atom_chain_id_ = atom_chain_id;
}

void PdbqtAtom::SetAtomResidueSequenceNumber(int atom_residue_sequence_number)
{
    atom_residue_sequence_number_ = atom_residue_sequence_number;
}

void PdbqtAtom::SetAtomInsertionCode(char atom_insertion_code)
{
    atom_insertion_code_ = atom_insertion_code;
}

void PdbqtAtom::SetAtomOrthogonalCoordinate(Geometry::Coordinate atom_orthogonal_coordinate)
{
    atom_orthogonal_coordinate_ = atom_orthogonal_coordinate;
}

void PdbqtAtom::SetAtomOccupancy(double atom_occupancy)
{
    atom_occupancy_ = atom_occupancy;
}

void PdbqtAtom::SetAtomTempretureFactor(double atom_temperature_factor)
{
    atom_temperature_factor_ = atom_temperature_factor;
}

void PdbqtAtom::SetAtomCharge(double atom_charge)
{
    atom_charge_ = atom_charge;
}
void PdbqtAtom::SetAtomType(string atom_type)
{
    atom_type_ = atom_type;
}

//////////////////////////////////////////////////////////
//                       DISPLAY FUNCTION               //
//////////////////////////////////////////////////////////
void PdbqtAtom::Print(ostream &out)
{
}

