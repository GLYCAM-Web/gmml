#include "includes/InputSet/PdbFile/atomRecord.hpp"
#include "includes/common.hpp"
#include "includes/CodeUtils/strings.hpp"
#include "includes/CodeUtils/logging.hpp"

using pdb::AtomRecord;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////

AtomRecord::AtomRecord(const std::string &line, int modelNumber) : serialNumber_(gmml::iNotSet), modelNumber_(gmml::iNotSet), atomName_(""), alternateLocation_(' '), residueName_(""), chainId_(' '), residueSequenceNumber_(gmml::iNotSet), insertionCode_(' '), occupancy_(gmml::dNotSet), temperatureFactor_(gmml::dNotSet), element_(""), charge_(""), residueSequenceIndex_("")
{
    this->SetModelNumber(modelNumber);
    this->SetRecordName(codeUtils::RemoveWhiteSpace(line.substr(0,6)));
    try
    {
        serialNumber_ = std::stoi(codeUtils::RemoveWhiteSpace(line.substr(6, 5)));
    }
    catch (...)
    {
        gmml::log(__LINE__, __FILE__, gmml::ERR, "Error converting to serialNumber from: " + line.substr(6,5));
    }
    atomName_ = codeUtils::RemoveWhiteSpace(line.substr(12, 4));
    // OG Dec2021, we don't want alt locations in gmml for now. Messes with the preprocessor.
    alternateLocation_ = gmml::BLANK_SPACE;
    residueName_ = codeUtils::RemoveWhiteSpace(line.substr(17,3));
    std::string temp = codeUtils::RemoveWhiteSpace(line.substr(21,1));
    chainId_ = temp[0];
    residueSequenceNumber_ = std::stoi(codeUtils::RemoveWhiteSpace(line.substr(22,4)));
    temp = codeUtils::RemoveWhiteSpace(line.substr(26,1));
    insertionCode_ = temp[0];

    temp = codeUtils::RemoveWhiteSpace(line.substr(30,8));
    if(!temp.empty())
    {
        coordinate_.SetX(std::stod(temp));
    }
    temp = codeUtils::RemoveWhiteSpace(line.substr(38,8));
    if(!temp.empty())
    {
        coordinate_.SetY(std::stod(temp));
    }
    temp = codeUtils::RemoveWhiteSpace(line.substr(46,8));
    if(!temp.empty())
    {
        coordinate_.SetZ(std::stod(temp));
    }
    occupancy_ = std::stod(codeUtils::RemoveWhiteSpace(line.substr(54,6)));
    temperatureFactor_ = std::stod(codeUtils::RemoveWhiteSpace(line.substr(60,6)));
    temp = codeUtils::RemoveWhiteSpace(line.substr(76, 2));
    if(temp.size() == 2)
    {
      temp[1] = std::tolower(temp[1]);
    }
    element_ = temp;
    charge_ = codeUtils::RemoveWhiteSpace(line.substr(78, 2));
}
/////////////////////////////////////////////////////////
//                       MUTATOR                        //
//////////////////////////////////////////////////////////
void AtomRecord::SetModelNumber(const int i)
{
    modelNumber_ = i;
}
void AtomRecord::SetRecordName(const std::string s)
{
    recordName_ = s;
}
void AtomRecord::SetAtomSerialNumber(int atom_serial_number){
    serialNumber_ = atom_serial_number;
}

void AtomRecord::SetAtomName(const std::string atom_name){
    atomName_ = atom_name;
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
    out << "Serial Number: ";
    if(serialNumber_ == gmml::iNotSet)
        out << " ";
    else
        out << serialNumber_;
    out << ", Atom Name: " << atomName_
        << ", Alternate Location: " << alternateLocation_
        << ", Residue Name: " << residueName_
        << ", Chain ID: " << chainId_
        << ", Residue Sequence Number: ";
    if(residueSequenceNumber_ == gmml::iNotSet)
        out << " ";
    else
        out << residueSequenceNumber_;
    out << ", Inserion Code: " << insertionCode_
        << ", Coordinate: ";
    coordinate_.Print(out);
    out << ", Occupancy: ";
    if(occupancy_ == gmml::dNotSet)
        out << " ";
    else
        out << occupancy_;
    out << ", Temperature Factor: ";
    if(temperatureFactor_ == gmml::dNotSet)
        out << " ";
    else
        out << temperatureFactor_;
    out << ", Element: " << element_
        << ", Charge: " << charge_ << std::endl;
}
