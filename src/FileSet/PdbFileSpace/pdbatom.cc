#include "../../../includes/FileSet/PdbFileSpace/pdbatom.hpp"
#include "../../../includes/utils.hpp"
#include "../../../includes/common.hpp"

using namespace std;
using namespace PdbFileSpace;
using namespace gmml;
using namespace Geometry;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbAtom::PdbAtom():atom_orthogonal_coordinate_() {}

PdbAtom::PdbAtom(string &line)
{
    if(line.substr(6,5) == "     ")
    {
        atom_serial_number_ =  iNotSet;
    }
    else
    {
        atom_serial_number_ = ConvertString<int>(line.substr(6,5));
    }

    atom_name_ = line.substr(12, 4);
    Trim(atom_name_);

    if(line.substr(16,1) == " ")
    {
        atom_alternate_location_ = ' ';
    }
    else
    {
        atom_alternate_location_ = ConvertString<char>(line.substr(16,1));
    }

    if(line.substr(17,3) == "   ")
    {
        atom_residue_name_ = " ";
    }
    else
    {
        atom_residue_name_ = line.substr(17,3);
        atom_residue_name_ = Trim(atom_residue_name_);
    }

    if(line.substr(21,1) == " ")
    {
        atom_chain_id_ = ' ';
    }
    else
    {
        atom_chain_id_ = ConvertString<char>(line.substr(21, 1));
    }

    if(line.substr(22,4) == "    ")
    {
        atom_residue_sequence_number_ = iNotSet;
    }
    else
    {
        atom_residue_sequence_number_ = ConvertString<int>(line.substr(22, 4));
    }

    if(line.substr(26,1) == " ")
    {
        atom_insertion_code_ = ' ';
    }
    else
    {
        atom_insertion_code_ = ConvertString<char>(line.substr(26, 1));
    }

    if(line.substr(30,8) == "        ")
    {
        atom_orthogonal_coordinate_.SetX(dNotSet);
    }
    else
    {
        atom_orthogonal_coordinate_.SetX(ConvertString<double>(line.substr(30, 8)));
    }

    if(line.substr(38,8) == "        ")
    {
        atom_orthogonal_coordinate_.SetY(dNotSet);
    }
    else
    {
        atom_orthogonal_coordinate_.SetY( ConvertString<double>(line.substr(38,8)));
    }

    if(line.substr(46,8) == "        ")
    {
        atom_orthogonal_coordinate_.SetZ(dNotSet);
    }
    else
    {
        atom_orthogonal_coordinate_.SetZ( ConvertString<double>(line.substr(46,8)));
    }

    if(line.substr(54, 6) == "      ")
    {
        atom_occupancy_ = dNotSet;
    }
    else
    {
        atom_occupancy_ = ConvertString<double>(line.substr(54, 6));
    }

    if(line.substr(60, 6) == "      ")
    {
        atom_temperature_factor_ = dNotSet;
    }
    else
    {
        atom_temperature_factor_ = ConvertString<double>(line.substr(60, 6));
    }

    atom_element_symbol_ = line.substr(76, 2);
    Trim(atom_element_symbol_);

    atom_charge_ = line.substr(78, 2);
    Trim(atom_charge_);
}
PdbAtom::PdbAtom(char residue_chain_id, string atom_name, string residue_name, int residue_sequence_number, char residue_insertion_code, char atom_alternate_location) :
    atom_chain_id_(residue_chain_id), atom_name_(atom_name), atom_residue_name_(residue_name), atom_residue_sequence_number_(residue_sequence_number),
    atom_insertion_code_(residue_insertion_code), atom_alternate_location_(atom_alternate_location) {}

PdbAtom::PdbAtom(int atom_serial_number, string atom_name, char atom_alternate_location, string residue_name, char chain_id,
                 int residue_sequence_number, char insertion_code, Coordinate coordinate, double occupancy, double tempreture_factor,
                 string element_symbol, string charge) :
    atom_serial_number_(atom_serial_number), atom_name_(atom_name), atom_alternate_location_(atom_alternate_location), atom_residue_name_(residue_name),
    atom_chain_id_(chain_id), atom_residue_sequence_number_(residue_sequence_number), atom_insertion_code_(insertion_code), atom_orthogonal_coordinate_(coordinate),
    atom_occupancy_(occupancy), atom_temperature_factor_(tempreture_factor), atom_element_symbol_(element_symbol), atom_charge_(charge) {}

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
    return atom_temperature_factor_;
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

void PdbAtom::SetAtomTempretureFactor(double atom_temperature_factor){
    atom_temperature_factor_ = atom_temperature_factor;
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
void PdbAtom::Print(ostream &out)
{
    out << "Atom Serial Number: ";
    if(atom_serial_number_ == iNotSet)
        out << " ";
    else
        out << atom_serial_number_;
    out << ", Atom Name: " << atom_name_
        << ", Atom Alternate Location: " << atom_alternate_location_
        << ", Atom Residue Name: " << atom_residue_name_
        << ", Atom Chain ID: " << atom_chain_id_
        << ", Atom Residue Sequence Number: ";
    if(atom_residue_sequence_number_ == iNotSet)
        out << " ";
    else
        out << atom_residue_sequence_number_;
    out << ", Atom Inserion Code: " << atom_insertion_code_
        << ", Atom Orthogonal Coordinate: ";
    atom_orthogonal_coordinate_.Print(out);
    out << ", Atom Occupancy: ";
    if(atom_occupancy_ == dNotSet)
        out << " ";
    else
        out << atom_occupancy_;
    out << ", Atom Tempreture Factor: ";
    if(atom_temperature_factor_ == dNotSet)
        out << " ";
    else
        out << atom_temperature_factor_;
    out << ", Atom Element Symbol: " << atom_element_symbol_
        << ", Atom Charge: " << atom_charge_ << endl;
}
