#include "includes/InputSet/PdbFile/atomRecord.hpp"
#include "includes/common.hpp"
#include "includes/CodeUtils/strings.hpp"
#include "includes/CodeUtils/logging.hpp"

using pdb::AtomRecord;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
AtomRecord::AtomRecord() : modelNumber_(1), recordName_(""), serialNumber_(gmml::iNotSet), atomName_(""), alternateLocation_(""), residueName_(""), chainId_(""), residueSequenceNumber_(gmml::iNotSet), insertionCode_(""), occupancy_(gmml::dNotSet), temperatureFactor_(gmml::dNotSet), element_(""), charge_("") {}

AtomRecord::AtomRecord(const std::string &line, int modelNumber) : modelNumber_(modelNumber), recordName_(""), serialNumber_(gmml::iNotSet), atomName_(""), alternateLocation_(""), residueName_(""), chainId_(""), residueSequenceNumber_(gmml::iNotSet), insertionCode_(""), occupancy_(gmml::dNotSet), temperatureFactor_(gmml::dNotSet), element_(""), charge_("")
{
    //gmml::log(__LINE__, __FILE__, gmml::INF, "Parsing " + line);
    // In the PDB file the residue number overruns after 9999 and serial number overruns after 99999. First overun for serial doesn't matter as there should be a space between the number and the name. So the problem is above 999999
    this->SetModelNumber(modelNumber);
    this->SetRecordName(codeUtils::RemoveWhiteSpace(line.substr(0,6)));
    // Dealing with number overruns for serialNumber and residueNumber
    int shift = 0;
    shift = codeUtils::GetSizeOfIntInString(line.substr(12));
    try
    {
        serialNumber_ = std::stoi(codeUtils::RemoveWhiteSpace(line.substr(6, 6 + shift)));
    }
    catch (...)
    {
        gmml::log(__LINE__, __FILE__, gmml::ERR, "Error converting to serialNumber from: " + line.substr(6, 6 + shift));
    }
    atomName_ = codeUtils::RemoveWhiteSpace(line.substr(12 + shift, 4));
    if (atomName_.empty())
    {
        atomName_ = gmml::BLANK_SPACE;
    }
    //alternateLocation_ = ""; // OG Dec2021, we don't want alt locations in gmml for now. Messes with the preprocessor.
    residueName_ = codeUtils::RemoveWhiteSpace(line.substr(17 + shift, 3));
    if (residueName_.empty())
    {
        residueName_ = gmml::BLANK_SPACE;
    }
    chainId_ = codeUtils::RemoveWhiteSpace(line.substr(21 + shift, 1));
    if (chainId_.empty())
    {
        chainId_ = gmml::BLANK_SPACE;
    }
    int secondShift = codeUtils::GetSizeOfIntInString(line.substr(26 + shift));
    try
    {
        residueSequenceNumber_ = std::stoi(codeUtils::RemoveWhiteSpace(line.substr(22 + shift, 4 + secondShift)));
    }
    catch (...)
    {
        gmml::log(__LINE__, __FILE__, gmml::ERR, "Error setting residue number from this line:\n" + line);
    }
    // Insertion code gets shifted right by every overrun in residue number.
    insertionCode_ = codeUtils::RemoveWhiteSpace(line.substr(26 + shift + secondShift, 1));
    if (insertionCode_.empty())
    {
        insertionCode_ = gmml::BLANK_SPACE;
    }
    // Coordinates etc don't get shifted by first overrun in residue number, but do by the rest.
    if (secondShift > 1)
    {
        shift += (secondShift - 1); // Combine the shifts, but ignore the first shift in residue sequence number from here on
    }
    try
    {
        coordinate_ = GeometryTopology::Coordinate(codeUtils::RemoveWhiteSpace(line.substr(30 + shift, 8)), codeUtils::RemoveWhiteSpace(line.substr(38 + shift, 8)), codeUtils::RemoveWhiteSpace(line.substr(46 + shift, 8)));
    }
    catch (...)
    {
        gmml::log(__LINE__, __FILE__, gmml::ERR, "Error setting coordinate from this line:\n" + line);
    }
    try
    {
        occupancy_ = std::stod(codeUtils::RemoveWhiteSpace(line.substr(54 + shift, 6)));
    }
    catch (...)
    {
        gmml::log(__LINE__, __FILE__, gmml::ERR, "Error converting to occupancy from: " + line.substr(54 + shift, 6));
    }
    try
    {
        temperatureFactor_ = std::stod(codeUtils::RemoveWhiteSpace(line.substr(60 + shift, 6)));
    }
    catch (...)
    {
        gmml::log(__LINE__, __FILE__, gmml::ERR, "Error converting to temperatureFactor_ from: " + line.substr(60 + shift, 6));
    }
    // Worried about case where shift is 2 and these don't exist, so line length is 80. More Ifs: Ugh.
    if (shift <= 2)
    {
        element_ = codeUtils::RemoveWhiteSpace(line.substr(76 + shift, 2));; // this used to be for element: temp[1] = std::tolower(temp[1]);
    }
    if (shift == 0)
    {
        charge_ = codeUtils::RemoveWhiteSpace(line.substr(78, 2));
    }
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
void AtomRecord::SetSerialNumber(const int atom_serial_number)
{
    serialNumber_ = atom_serial_number;
}
void AtomRecord::SetAtomName(const std::string atom_name)
{
    atomName_ = atom_name;
}
void AtomRecord::SetAlternateLocation(const std::string atom_alternate_location)
{
    alternateLocation_ = atom_alternate_location;
}
void AtomRecord::SetResidueName(const std::string atom_residue_name)
{
    residueName_ = atom_residue_name;
}
void AtomRecord::SetChainId(const std::string atom_chain_id)
{
    chainId_ = atom_chain_id;
}
void AtomRecord::SetResidueSequenceNumber(const int atom_residue_sequence_number)
{
    residueSequenceNumber_ = atom_residue_sequence_number;
}
void AtomRecord::SetInsertionCode(const std::string atom_insertion_code)
{
    insertionCode_ = atom_insertion_code;
}
void AtomRecord::SetCoordinate(const GeometryTopology::Coordinate c)
{
    coordinate_ = c;
}
void AtomRecord::SetOccupancy(const double atom_occupancy)
{
    occupancy_ = atom_occupancy;
}
void AtomRecord::SetTempretureFactor(const double atom_temperature_factor)
{
    temperatureFactor_ = atom_temperature_factor;
}
void AtomRecord::SetElement(const std::string atom_element_symbol)
{
    element_ = atom_element_symbol;
}
void AtomRecord::SetCharge(const std::string atom_charge)
{
    charge_ = atom_charge;
}
//void AtomRecord::AddAlternateLocation(AtomRecord* alternate_atom)
//{
//  alternateLocations_.push_back(alternate_atom);
//}
//////////////////////////////////////////////////////////
//                       FUNCTION                       //
//////////////////////////////////////////////////////////
std::string AtomRecord::GetId() const
{
    std::stringstream ss;
    ss << this->GetName() << "_" << this->GetSerialNumber() << "_" << this->GetResidueId();
    return ss.str();
}

std::string AtomRecord::GetResidueId() const
{
    std::stringstream ss;
    ss << this->GetResidueName() << "_" << this->GetInsertionCode() << "_" << this->GetResidueSequenceNumber() << "_" + this->GetChainId();
    return ss.str();
}
double AtomRecord::CalculateDistance(const AtomRecord* otherAtom) const
{
    return this->GetCoordinate().Distance(otherAtom->GetCoordinate());
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
