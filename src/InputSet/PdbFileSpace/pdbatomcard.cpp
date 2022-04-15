#include "../../../includes/InputSet/PdbFileSpace/pdbatomcard.hpp"
#include "../../../includes/utils.hpp"
#include "../../../includes/common.hpp"

using PdbFileSpace::PdbAtomCard;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbAtomCard::PdbAtomCard():atom_orthogonal_coordinate_() {}

PdbAtomCard::PdbAtomCard(std::string &line)
{
    std::string temp = line.substr(6, 5);
    temp = gmml::Trim(temp);
    if(temp.empty())
    {
        atom_serial_number_ =  gmml::iNotSet;
    }
    else
    {
        atom_serial_number_ = gmml::ConvertString<int>(temp);
    }

    atom_name_ = line.substr(12, 4);
    gmml::Trim(atom_name_);

    temp = line.substr(16,1);
    temp = gmml::Trim(temp);
    if(temp.empty())
    {
        atom_alternate_location_ = gmml::BLANK_SPACE;
    }
    else
    {
        atom_alternate_location_ = gmml::ConvertString<char>(temp);
    }

    temp = line.substr(17,3);
    temp = gmml::Trim(temp);
    if(temp.empty())
    {
        atom_residue_name_ = " ";
    }
    else
    {
        atom_residue_name_ = temp;
    }

    temp = line.substr(21,1);
    temp = gmml::Trim(temp);
    if(temp.empty() || temp.compare("_") == 0)
    {
        atom_chain_id_ = gmml::BLANK_SPACE;
    }
    else
    {
        atom_chain_id_ = gmml::ConvertString<char>(temp);
    }

    temp = line.substr(22,4);
    temp = gmml::Trim(temp);
    if(temp.empty())
    {
        atom_residue_sequence_number_ = gmml::iNotSet;
    }
    else
    {
        atom_residue_sequence_number_ = gmml::ConvertString<int>(temp);
    }

    temp = line.substr(26,1);
    temp = gmml::Trim(temp);
    if(temp.empty())
    {
        atom_insertion_code_ = gmml::BLANK_SPACE;
    }
    else
    {
        atom_insertion_code_ = gmml::ConvertString<char>(temp);
    }

    temp = line.substr(30, 8);
    temp = gmml::Trim(temp);
    if(temp.empty())
    {
        atom_orthogonal_coordinate_.SetX(gmml::dNotSet);
    }
    else
    {
        atom_orthogonal_coordinate_.SetX(gmml::ConvertString<double>(temp));
    }

    temp = line.substr(38,8);
    temp = gmml::Trim(temp);
    if(temp.empty())
    {
        atom_orthogonal_coordinate_.SetY(gmml::dNotSet);
    }
    else
    {
        atom_orthogonal_coordinate_.SetY( gmml::ConvertString<double>(line.substr(38,8)));
    }

    temp = line.substr(46,8);
    temp = gmml::Trim(temp);
    if(temp.empty())
    {
        atom_orthogonal_coordinate_.SetZ(gmml::dNotSet);
    }
    else
    {
        atom_orthogonal_coordinate_.SetZ( gmml::ConvertString<double>(temp));
    }

    temp = line.substr(54, 6);
    temp = gmml::Trim(temp);
    if(temp.empty())
    {
        atom_occupancy_ = gmml::dNotSet;
    }
    else
    {
        atom_occupancy_ = gmml::ConvertString<double>(temp);
    }

    temp = line.substr(60, 6);
    temp = gmml::Trim(temp);
    if(temp.empty())
    {
        atom_temperature_factor_ = gmml::dNotSet;
    }
    else
    {
        atom_temperature_factor_ = gmml::ConvertString<double>(temp);
    }

    temp = line.substr(76, 2);
    temp = gmml::Trim(temp);
    if(temp.size() == 2)
    {
      temp[1] = std::tolower(temp[1]);
    }
    atom_element_symbol_ = temp;

    temp = line.substr(78, 2);
    temp = gmml::Trim(temp);
    atom_charge_ = temp;
}
PdbAtomCard::PdbAtomCard(char residue_chain_id, std::string atom_name,
                        std::string residue_name, int residue_sequence_number,
                        char residue_insertion_code, char atom_alternate_location, std::vector<PdbAtomCard*> alternate_atom_locations) :
    atom_name_(atom_name), atom_alternate_location_(atom_alternate_location),
    atom_residue_name_(residue_name), atom_chain_id_(residue_chain_id),
    atom_residue_sequence_number_(residue_sequence_number), atom_insertion_code_(residue_insertion_code),
    alternate_atom_locations_(alternate_atom_locations) {}

PdbAtomCard::PdbAtomCard(int atom_serial_number, std::string atom_name, char atom_alternate_location,
                        std::string residue_name, char chain_id, int residue_sequence_number,
                        char insertion_code, GeometryTopology::Coordinate coordinate, double occupancy,
                        double tempreture_factor, std::string element_symbol, std::string charge, std::vector<PdbAtomCard*> alternate_atom_locations) :
    atom_serial_number_(atom_serial_number), atom_name_(atom_name), atom_alternate_location_(atom_alternate_location), atom_residue_name_(residue_name),
    atom_chain_id_(chain_id), atom_residue_sequence_number_(residue_sequence_number), atom_insertion_code_(insertion_code),
    atom_orthogonal_coordinate_(coordinate), atom_occupancy_(occupancy), atom_temperature_factor_(tempreture_factor), atom_element_symbol_(element_symbol),
    atom_charge_(charge), alternate_atom_locations_(alternate_atom_locations) {}

//////////////////////////////////////////////////////////
//                       ACCESSOR                       //
//////////////////////////////////////////////////////////
int PdbAtomCard::GetAtomSerialNumber(){
    return atom_serial_number_;
}

std::string PdbAtomCard::GetAtomName(){
    return atom_name_;
}

char PdbAtomCard::GetAtomAlternateLocation(){
    return atom_alternate_location_;
}

std::string PdbAtomCard::GetAtomResidueName(){
    return atom_residue_name_;
}

char PdbAtomCard::GetAtomChainId(){
    return atom_chain_id_;
}

int PdbAtomCard::GetAtomResidueSequenceNumber(){
    return atom_residue_sequence_number_;
}

char PdbAtomCard::GetAtomInsertionCode(){
    return atom_insertion_code_;
}

GeometryTopology::Coordinate PdbAtomCard::GetAtomOrthogonalCoordinate(){
    return atom_orthogonal_coordinate_;
}

double PdbAtomCard::GetAtomOccupancy(){
    return atom_occupancy_;
}

double PdbAtomCard::GetAtomTempretureFactor(){
    return atom_temperature_factor_;
}

std::string PdbAtomCard::GetAtomElementSymbol(){
    return atom_element_symbol_;
}

std::string PdbAtomCard::GetAtomCharge(){
    return atom_charge_;
}
std::string PdbAtomCard::GetAtomCardIndexInResidueSet()
{
    return atom_card_index_in_residue_sequence_;
}

std::vector<PdbAtomCard*> PdbAtomCard::GetAlternateAtomCards()
{
  return alternate_atom_locations_;
}

//////////////////////////////////////////////////////////
//                       MUTATOR                        //
//////////////////////////////////////////////////////////
void PdbAtomCard::SetAtomSerialNumber(int atom_serial_number){
    atom_serial_number_ = atom_serial_number;
}

void PdbAtomCard::SetAtomName(const std::string atom_name){
    atom_name_ = atom_name;
}

void PdbAtomCard::SetAtomAlternateLocation(char atom_alternate_location){
    atom_alternate_location_ = atom_alternate_location;
}

void PdbAtomCard::SetAtomResidueName(const std::string atom_residue_name){
    atom_residue_name_ = atom_residue_name;
}

void PdbAtomCard::SetAtomChainId(char atom_chain_id){
    atom_chain_id_ = atom_chain_id;
}

void PdbAtomCard::SetAtomResidueSequenceNumber(int atom_residue_sequence_number){
    atom_residue_sequence_number_ = atom_residue_sequence_number;
}

void PdbAtomCard::SetAtomInsertionCode(char atom_insertion_code){
    atom_insertion_code_ = atom_insertion_code;
}

void PdbAtomCard::SetAtomOrthogonalCoordinate(GeometryTopology::Coordinate atom_orthogonal_coordinate){
    atom_orthogonal_coordinate_ = atom_orthogonal_coordinate;
}

void PdbAtomCard::SetAtomOccupancy(double atom_occupancy){
    atom_occupancy_ = atom_occupancy;
}

void PdbAtomCard::SetAtomTempretureFactor(double atom_temperature_factor){
    atom_temperature_factor_ = atom_temperature_factor;
}

void PdbAtomCard::SetAtomElementSymbol(const std::string atom_element_symbol){
    atom_element_symbol_ = atom_element_symbol;
}

void PdbAtomCard::SetAtomCharge(const std::string atom_charge){
    atom_charge_ = atom_charge;
}
void PdbAtomCard::SetAtomCardIndexInResidueSet(std::string atom_card_index_in_residue_sequence)
{
    atom_card_index_in_residue_sequence_ = atom_card_index_in_residue_sequence;
}

void PdbAtomCard::AddAlternateLocation(PdbAtomCard* alternate_atom)
{
  alternate_atom_locations_.push_back(alternate_atom);
}

//////////////////////////////////////////////////////////
//                       DISPLAY FUNCTION               //
//////////////////////////////////////////////////////////
void PdbAtomCard::Print(std::ostream &out)
{
    out << "Atom Serial Number: ";
    if(atom_serial_number_ == gmml::iNotSet)
        out << " ";
    else
        out << atom_serial_number_;
    out << ", Atom Name: " << atom_name_
        << ", Atom Alternate Location: " << atom_alternate_location_
        << ", Atom Residue Name: " << atom_residue_name_
        << ", Atom Chain ID: " << atom_chain_id_
        << ", Atom Residue Sequence Number: ";
    if(atom_residue_sequence_number_ == gmml::iNotSet)
        out << " ";
    else
        out << atom_residue_sequence_number_;
    out << ", Atom Inserion Code: " << atom_insertion_code_
        << ", Atom Orthogonal Coordinate: ";
    atom_orthogonal_coordinate_.Print(out);
    out << ", Atom Occupancy: ";
    if(atom_occupancy_ == gmml::dNotSet)
        out << " ";
    else
        out << atom_occupancy_;
    out << ", Atom Tempreture Factor: ";
    if(atom_temperature_factor_ == gmml::dNotSet)
        out << " ";
    else
        out << atom_temperature_factor_;
    out << ", Atom Element Symbol: " << atom_element_symbol_
        << ", Atom Charge: " << atom_charge_ << std::endl;
}
