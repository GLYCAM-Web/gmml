#include "../../../includes/FileSet/PdbFileSpace/pdbatom.hpp"
//#include "../../../includes/Geometry/coordinate.hpp"

using namespace std;
using namespace PdbFileSpace;
using namespace Geometry;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbAtom::PdbAtom():atom_orthogonal_coordinate_() {}


//////////////////////////////////////////////////////////
//                       ACCESSOR                       //
//////////////////////////////////////////////////////////
int PdbAtom::GetAtomSerialNumber(){
    return atom_serial_number_;
}

string PdbAtom::GetAtomName(){
    return atom_name_;
}

char PdbAtom::GetAtomAlternateLocation(){
    return atom_alternate_location_;
}

string PdbAtom::GetAtomResidueName(){
    return atom_residue_name_;
}

char PdbAtom::GetAtomChainId(){
    return atom_chain_id_;
}

int PdbAtom::GetAtomResidueSequenceNumber(){
    return atom_residue_sequence_number_;
}

char PdbAtom::GetAtomInsertionCode(){
    return atom_insertion_code_;
}

Geometry::Coordinate PdbAtom::GetAtomOrthogonalCoordinate(){
    return atom_orthogonal_coordinate_;
}

double PdbAtom::GetAtomOccupancy(){
    return atom_occupancy_;
}

double PdbAtom::GetAtomTempretureFactor(){
    return atom_tempreture_factor_;
}

string PdbAtom::GetAtomElementSymbol(){
    return atom_element_symbol_;
}

string PdbAtom::GetAtomCharge(){
    return atom_charge_;
}

//////////////////////////////////////////////////////////
//                       MUTATOR                        //
//////////////////////////////////////////////////////////
void PdbAtom::SetAtomSerialNumber(int atom_serial_number){
    atom_serial_number_ = atom_serial_number;
}

void PdbAtom::SetAtomName(const string atom_name){
    atom_name_ = atom_name;
}

void PdbAtom::SetAtomAlternateLocation(char atom_alternate_location){
    atom_alternate_location_ = atom_alternate_location;
}

void PdbAtom::SetAtomResidueName(const string atom_residue_name){
    atom_residue_name_ = atom_residue_name;
}

void PdbAtom::SetAtomChainId(char atom_chain_id){
    atom_chain_id_ = atom_chain_id;
}

void PdbAtom::SetAtomResidueSequenceNumber(int atom_residue_sequence_number){
    atom_residue_sequence_number_ = atom_residue_sequence_number;
}

void PdbAtom::SetAtomInsertionCode(char atom_insertion_code){
    atom_insertion_code_ = atom_insertion_code;
}

void PdbAtom::SetAtomOrthogonalCoordinate(Geometry::Coordinate atom_orthogonal_coordinate){
    atom_orthogonal_coordinate_ = atom_orthogonal_coordinate;
}

void PdbAtom::SetAtomOccupancy(double atom_occupancy){
    atom_occupancy_ = atom_occupancy;
}

void PdbAtom::SetAtomTempretureFactor(double atom_tempreture_factor){
    atom_tempreture_factor_ = atom_tempreture_factor;
}

void PdbAtom::SetAtomElementSymbol(const string atom_element_symbol){
    atom_element_symbol_ = atom_element_symbol;
}

void PdbAtom::SetAtomCharge(const string atom_charge){
    atom_charge_ = atom_charge;
}

//////////////////////////////////////////////////////////
//                       DISPLAY FUNCTION               //
//////////////////////////////////////////////////////////


