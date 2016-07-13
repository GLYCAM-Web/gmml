#include "../../../includes/InputSet/PdbFileSpace/pdbatom.hpp"
#include "../../../includes/utils.hpp"
#include "../../../includes/common.hpp"

using namespace std;
using namespace PdbFileSpace;
using namespace gmml;
using namespace GeometryTopology;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbAtom::PdbAtom():atom_orthogonal_coordinate_() {}

PdbAtom::PdbAtom(string &line)
{
    string temp = line.substr(6, 5);
    temp = Trim(temp);
    if(temp.empty())
    {
        atom_serial_number_ =  iNotSet;
    }
    else
    {
        atom_serial_number_ = ConvertString<int>(temp);
    }

    atom_name_ = line.substr(12, 4);
    Trim(atom_name_);

    temp = line.substr(16,1);
    temp = Trim(temp);
    if(temp.empty())
    {
        atom_alternate_location_ = BLANK_SPACE;
    }
    else
    {
        atom_alternate_location_ = ConvertString<char>(temp);
    }

    temp = line.substr(17,3);
    temp = Trim(temp);
    if(temp.empty())
    {
        atom_residue_name_ = " ";
    }
    else
    {
        atom_residue_name_ = temp;
    }

    temp = line.substr(21,1);
    temp = Trim(temp);
    if(temp.empty() || temp.compare("_") == 0)
    {
        atom_chain_id_ = BLANK_SPACE;
    }
    else
    {
        atom_chain_id_ = ConvertString<char>(temp);
    }

    temp = line.substr(22,4);
    temp = Trim(temp);
    if(temp.empty())
    {
        atom_residue_sequence_number_ = iNotSet;
    }
    else
    {
        atom_residue_sequence_number_ = ConvertString<int>(temp);
    }

    temp = line.substr(26,1);
    temp = Trim(temp);
    if(temp.empty())
    {
        atom_insertion_code_ = BLANK_SPACE;
    }
    else
    {
        atom_insertion_code_ = ConvertString<char>(temp);
    }

    temp = line.substr(30, 8);
    temp = Trim(temp);
    if(temp.empty())
    {
        atom_orthogonal_coordinate_.SetX(dNotSet);
    }
    else
    {
        atom_orthogonal_coordinate_.SetX(ConvertString<double>(temp));
    }

    temp = line.substr(38,8);
    temp = Trim(temp);
    if(temp.empty())
    {
        atom_orthogonal_coordinate_.SetY(dNotSet);
    }
    else
    {
        atom_orthogonal_coordinate_.SetY( ConvertString<double>(line.substr(38,8)));
    }

    temp = line.substr(46,8);
    temp = Trim(temp);
    if(temp.empty())
    {
        atom_orthogonal_coordinate_.SetZ(dNotSet);
    }
    else
    {
        atom_orthogonal_coordinate_.SetZ( ConvertString<double>(temp));
    }

    temp = line.substr(54, 6);
    temp = Trim(temp);
    if(temp.empty())
    {
        atom_occupancy_ = dNotSet;
    }
    else
    {
        atom_occupancy_ = ConvertString<double>(temp);
    }

    temp = line.substr(60, 6);
    temp = Trim(temp);
    if(temp.empty())
    {
        atom_temperature_factor_ = dNotSet;
    }
    else
    {
        atom_temperature_factor_ = ConvertString<double>(temp);
    }

    temp = line.substr(76, 2);
    temp = Trim(temp);
    atom_element_symbol_ = temp;

    temp = line.substr(78, 2);
    temp = Trim(temp);
    atom_charge_ = temp;
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

GeometryTopology::Coordinate PdbAtom::GetAtomOrthogonalCoordinate(){
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
string PdbAtom::GetAtomCardIndexInResidueSet()
{
    return atom_card_index_in_residue_sequence_;
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

void PdbAtom::SetAtomOrthogonalCoordinate(GeometryTopology::Coordinate atom_orthogonal_coordinate){
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
void PdbAtom::SetAtomCardIndexInResidueSet(string atom_card_index_in_residue_sequence)
{
    atom_card_index_in_residue_sequence_ = atom_card_index_in_residue_sequence;
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
