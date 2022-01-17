#include "includes/InputSet/PdbFile/atomRecord.hpp"
#include "includes/utils.hpp"
#include "includes/common.hpp"

using pdb::AtomRecord;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////

AtomRecord::AtomRecord(const std::string &line, int modelNumber)
{
    this->SetModelNumber(modelNumber);
    std::string temp = line.substr(6, 5);
    temp = gmml::Trim(temp);
    if(temp.empty())
    {
        serialNumber_ =  gmml::iNotSet;
    }
    else
    {
        serialNumber_ = gmml::ConvertString<int>(temp);
    }

    name_ = line.substr(12, 4);
    gmml::Trim(name_);

    // OG Dec2021, we don't want alt locations in gmml for now. Messes with the preprocessor.
//    temp = line.substr(16,1);
//    temp = gmml::Trim(temp);
//    if(temp.empty())
//    {
    alternateLocation_ = gmml::BLANK_SPACE;
//    }
//    else
//    {
//        alternateLocation_ = gmml::ConvertString<char>(temp);
//    }

    temp = line.substr(17,3);
    temp = gmml::Trim(temp);
    if(temp.empty())
    {
        residueName_ = " ";
    }
    else
    {
        residueName_ = temp;
    }

    temp = line.substr(21,1);
    temp = gmml::Trim(temp);
    if(temp.empty() || temp.compare("_") == 0)
    {
        chainId_ = gmml::BLANK_SPACE;
    }
    else
    {
        chainId_ = gmml::ConvertString<char>(temp);
    }

    temp = line.substr(22,4);
    temp = gmml::Trim(temp);
    if(temp.empty())
    {
        residueSequenceNumber_ = gmml::iNotSet;
    }
    else
    {
        residueSequenceNumber_ = gmml::ConvertString<int>(temp);
    }

    temp = line.substr(26,1);
    temp = gmml::Trim(temp);
    if(temp.empty())
    {
        insertionCode_ = gmml::BLANK_SPACE;
    }
    else
    {
        insertionCode_ = gmml::ConvertString<char>(temp);
    }

    temp = line.substr(30, 8);
    temp = gmml::Trim(temp);
    if(temp.empty())
    {
        coordinate_.SetX(gmml::dNotSet);
    }
    else
    {
        coordinate_.SetX(gmml::ConvertString<double>(temp));
    }

    temp = line.substr(38,8);
    temp = gmml::Trim(temp);
    if(temp.empty())
    {
        coordinate_.SetY(gmml::dNotSet);
    }
    else
    {
        coordinate_.SetY( gmml::ConvertString<double>(line.substr(38,8)));
    }

    temp = line.substr(46,8);
    temp = gmml::Trim(temp);
    if(temp.empty())
    {
        coordinate_.SetZ(gmml::dNotSet);
    }
    else
    {
        coordinate_.SetZ( gmml::ConvertString<double>(temp));
    }

    temp = line.substr(54, 6);
    temp = gmml::Trim(temp);
    if(temp.empty())
    {
        occupancy_ = gmml::dNotSet;
    }
    else
    {
        occupancy_ = gmml::ConvertString<double>(temp);
    }

    temp = line.substr(60, 6);
    temp = gmml::Trim(temp);
    if(temp.empty())
    {
        temperatureFactor_ = gmml::dNotSet;
    }
    else
    {
        temperatureFactor_ = gmml::ConvertString<double>(temp);
    }

    temp = line.substr(76, 2);
    temp = gmml::Trim(temp);
    if(temp.size() == 2)
    {
      temp[1] = std::tolower(temp[1]);
    }
    element_ = temp;

    temp = line.substr(78, 2);
    temp = gmml::Trim(temp);
    charge_ = temp;
}
/////////////////////////////////////////////////////////
//                       MUTATOR                        //
//////////////////////////////////////////////////////////
void AtomRecord::SetModelNumber(const int i)
{
    modelNumber_ = i;
}
void AtomRecord::SetAtomSerialNumber(int atom_serial_number){
    serialNumber_ = atom_serial_number;
}

void AtomRecord::SetAtomName(const std::string atom_name){
    name_ = atom_name;
}

void AtomRecord::SetAtomAlternateLocation(char atom_alternate_location){
    alternateLocation_ = atom_alternate_location;
}

void AtomRecord::SetAtomResidueName(const std::string atom_residue_name){
    residueName_ = atom_residue_name;
}

void AtomRecord::SetAtomChainId(char atom_chain_id){
    chainId_ = atom_chain_id;
}

void AtomRecord::SetAtomResidueSequenceNumber(int atom_residue_sequence_number){
    residueSequenceNumber_ = atom_residue_sequence_number;
}

void AtomRecord::SetAtomInsertionCode(char atom_insertion_code){
    insertionCode_ = atom_insertion_code;
}

void AtomRecord::SetAtomOrthogonalCoordinate(GeometryTopology::Coordinate atom_orthogonal_coordinate){
    coordinate_ = atom_orthogonal_coordinate;
}

void AtomRecord::SetAtomOccupancy(double atom_occupancy){
    occupancy_ = atom_occupancy;
}

void AtomRecord::SetAtomTempretureFactor(double atom_temperature_factor){
    temperatureFactor_ = atom_temperature_factor;
}

void AtomRecord::SetAtomElementSymbol(const std::string atom_element_symbol){
    element_ = atom_element_symbol;
}

void AtomRecord::SetAtomCharge(const std::string atom_charge){
    charge_ = atom_charge;
}
void AtomRecord::SetAtomCardIndexInResidueSet(std::string atom_card_index_in_residue_sequence)
{
    residueSequenceIndex_ = atom_card_index_in_residue_sequence;
}

void AtomRecord::AddAlternateLocation(AtomRecord* alternate_atom)
{
  alternateLocations_.push_back(alternate_atom);
}

//////////////////////////////////////////////////////////
//                       DISPLAY FUNCTION               //
//////////////////////////////////////////////////////////
void AtomRecord::Print(std::ostream &out) const
{
    out << "Atom Serial Number: ";
    if(serialNumber_ == gmml::iNotSet)
        out << " ";
    else
        out << serialNumber_;
    out << ", Atom Name: " << name_
        << ", Atom Alternate Location: " << alternateLocation_
        << ", Atom Residue Name: " << residueName_
        << ", Atom Chain ID: " << chainId_
        << ", Atom Residue Sequence Number: ";
    if(residueSequenceNumber_ == gmml::iNotSet)
        out << " ";
    else
        out << residueSequenceNumber_;
    out << ", Atom Inserion Code: " << insertionCode_
        << ", Atom Orthogonal Coordinate: ";
    coordinate_.Print(out);
    out << ", Atom Occupancy: ";
    if(occupancy_ == gmml::dNotSet)
        out << " ";
    else
        out << occupancy_;
    out << ", Atom Tempreture Factor: ";
    if(temperatureFactor_ == gmml::dNotSet)
        out << " ";
    else
        out << temperatureFactor_;
    out << ", Atom Element Symbol: " << element_
        << ", Atom Charge: " << charge_ << std::endl;
}
