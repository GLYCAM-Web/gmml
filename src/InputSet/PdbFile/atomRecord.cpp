#include "includes/InputSet/PdbFile/atomRecord.hpp"
#include "includes/common.hpp"
#include "includes/CodeUtils/strings.hpp"
#include "includes/CodeUtils/logging.hpp"

using pdb::AtomRecord;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
AtomRecord::AtomRecord(const std::string &line, int modelNumber) : modelNumber_(modelNumber), recordName_(""), serialNumber_(gmml::iNotSet), atomName_(""), alternateLocation_(""), residueName_(""), chainId_(""), residueSequenceNumber_(gmml::iNotSet), insertionCode_(""), occupancy_(gmml::dNotSet), temperatureFactor_(gmml::dNotSet), element_(""), charge_("")
{
//    gmml::log(__LINE__, __FILE__, gmml::INF, "Parsing " + line);
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
    if (atomName_.empty())
    {
        atomName_ = gmml::BLANK_SPACE;
    }
    //alternateLocation_ = ""; // OG Dec2021, we don't want alt locations in gmml for now. Messes with the preprocessor.
    residueName_ = codeUtils::RemoveWhiteSpace(line.substr(17,3));
    if (residueName_.empty())
    {
        residueName_ = gmml::BLANK_SPACE;
    }
    chainId_ = codeUtils::RemoveWhiteSpace(line.substr(21,1));
    if (chainId_.empty())
    {
        chainId_ = gmml::BLANK_SPACE;
    }
    try
    {
        residueSequenceNumber_ = std::stoi(codeUtils::RemoveWhiteSpace(line.substr(22,4)));
    }
    catch (...)
    {
        gmml::log(__LINE__, __FILE__, gmml::ERR, "Error setting residue number from this line:\n" + line);
    }
    insertionCode_ = codeUtils::RemoveWhiteSpace(line.substr(26,1));
    if (insertionCode_.empty())
    {
        insertionCode_ = gmml::BLANK_SPACE;
    }
    try
    {
        coordinate_ = GeometryTopology::Coordinate(codeUtils::RemoveWhiteSpace(line.substr(30,8)), codeUtils::RemoveWhiteSpace(line.substr(38,8)), codeUtils::RemoveWhiteSpace(line.substr(46,8)));
    }
    catch (...)
    {
        gmml::log(__LINE__, __FILE__, gmml::ERR, "Error setting coordinate from this line:\n" + line);
    }
    try
    {
        occupancy_ = std::stod(codeUtils::RemoveWhiteSpace(line.substr(54,6)));
    }
    catch (...)
    {
        gmml::log(__LINE__, __FILE__, gmml::ERR, "Error converting to occupancy from: " + line.substr(54,6));
    }
    try
    {
        temperatureFactor_ = std::stod(codeUtils::RemoveWhiteSpace(line.substr(60,6)));
    }
    catch (...)
    {
        gmml::log(__LINE__, __FILE__, gmml::ERR, "Error converting to temperatureFactor_ from: " + line.substr(60,6));
    }
    element_ = codeUtils::RemoveWhiteSpace(line.substr(76, 2));; // this used to be for element: temp[1] = std::tolower(temp[1]);
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
std::string AtomRecord::GetResidueId() const
{
    std::stringstream ss;
    ss << this->GetResidueSequenceNumber() << "_" << this->GetInsertionCode() << "_" + this->GetChainId();
    return ss.str();
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
