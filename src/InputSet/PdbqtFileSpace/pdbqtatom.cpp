#include "../../../includes/InputSet/PdbqtFileSpace/pdbqtatom.hpp"
#include "../../../includes/utils.hpp"
#include "../../../includes/common.hpp"

using PdbqtFileSpace::PdbqtAtom;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbqtAtom::PdbqtAtom():atom_orthogonal_coordinate_() {}

PdbqtAtom::PdbqtAtom(std::string &line)
{
    std::string temp = line.substr(0, 6);
    temp = gmml::Trim(temp);
    if(temp.compare("ATOM") == 0)
        type_ = "ATOM";
    else if(temp.compare("HETATM") == 0)
        type_ = "HETATM";
    temp =line.substr(6, 5);
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
    if(temp.empty())
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

    temp = line.substr(70, 6);
    temp = gmml::Trim(temp);
    if(temp.empty())
    {
        atom_charge_ = gmml::dNotSet;
    }
    else
    {
        atom_charge_ = gmml::ConvertString<double>(temp);
    }

    temp = line.substr(77, 2);
    temp = gmml::Trim(temp);
    atom_type_ = temp;
}
PdbqtAtom::PdbqtAtom(int atom_serial_number, std::string atom_name, char atom_alternate_location, std::string residue_name, char chain_id,
                     int residue_sequence_number, char insertion_code, GeometryTopology::Coordinate coordinate, double occupancy, double temperature_factor,
                     double charge, std::string atom_type, std::string type) :
    atom_serial_number_(atom_serial_number), atom_name_(atom_name), atom_alternate_location_(atom_alternate_location), atom_residue_name_(residue_name),
    atom_chain_id_(chain_id), atom_residue_sequence_number_(residue_sequence_number), atom_insertion_code_(insertion_code),
    atom_orthogonal_coordinate_(coordinate), atom_occupancy_(occupancy), atom_temperature_factor_(temperature_factor), atom_charge_(charge),
    atom_type_(atom_type), type_(type){}

//////////////////////////////////////////////////////////
//                       ACCESSOR                       //
//////////////////////////////////////////////////////////
std::string PdbqtAtom::GetCardType()
{
    return card_type_;
}

int PdbqtAtom::GetAtomSerialNumber()
{
    return atom_serial_number_;
}

std::string PdbqtAtom::GetAtomName()
{
    return atom_name_;
}

char PdbqtAtom::GetAtomAlternateLocation()
{
    return atom_alternate_location_;
}

std::string PdbqtAtom::GetAtomResidueName()
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

GeometryTopology::Coordinate PdbqtAtom::GetAtomOrthogonalCoordinate()
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
std::string PdbqtAtom::GetAtomType()
{
    return atom_type_;
}
std::string PdbqtAtom::GetType()
{
    return type_;
}

//////////////////////////////////////////////////////////
//                       MUTATOR                        //
//////////////////////////////////////////////////////////
void PdbqtAtom::SetCardType(const std::string card_type)
{
    card_type_ = card_type;
}

void PdbqtAtom::SetAtomSerialNumber(int atom_serial_number)
{
    atom_serial_number_ = atom_serial_number;
}

void PdbqtAtom::SetAtomName(const std::string atom_name)
{
    atom_name_ = atom_name;
}

void PdbqtAtom::SetAtomAlternateLocation(char atom_alternate_location)
{
    atom_alternate_location_ = atom_alternate_location;
}

void PdbqtAtom::SetAtomResidueName(const std::string atom_residue_name)
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

void PdbqtAtom::SetAtomOrthogonalCoordinate(GeometryTopology::Coordinate atom_orthogonal_coordinate)
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
void PdbqtAtom::SetAtomType(std::string atom_type)
{
    atom_type_ = atom_type;
}
void PdbqtAtom::SetType(std::string type)
{
    type_ = type;
}

//////////////////////////////////////////////////////////
//                       DISPLAY FUNCTION               //
//////////////////////////////////////////////////////////
void PdbqtAtom::Print(std::ostream &out)
{
    out << "Record Name: " << type_
        << ", Atom Serial Number: ";
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
    out << ", Atom Charge: " << atom_charge_
        << ", Atom Type: " << atom_type_ << std::endl;
}
