#include "includes/InputSet/PdbFile/pdbAtom.hpp"
#include "includes/common.hpp"
#include "includes/CodeUtils/strings.hpp"
#include "includes/CodeUtils/logging.hpp"

using pdb::pdbAtom;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
pdbAtom::pdbAtom(const std::string& name, const Coordinate& coord)
: cdsAtom(name, coord)
{

}
pdbAtom::pdbAtom(const std::string &line)
{
    //gmml::log(__LINE__, __FILE__, gmml::INF, "Parsing " + line);
    // In the PDB file the residue number overruns after 9999 and serial number overruns after 99999. First overun for serial doesn't matter as there should be a space between the number and the name. So the problem is above 999999
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
        serialNumber_ = gmml::iNotSet;
    }
    std::string atomName = codeUtils::RemoveWhiteSpace(line.substr(12 + shift, 4));
    if (atomName.empty())
    {
        this->setName("    ");
    }
    else
    {
        this->setName(atomName);
    }
    gmml::log(__LINE__, __FILE__, gmml::INF, "Hi, my name is " + this->getName());
    //alternateLocation_ = ""; // OG Dec2021, we don't want alt locations in gmml for now. Messes with the preprocessor.
    residueName_ = codeUtils::RemoveWhiteSpace(line.substr(17 + shift, 3));
    if (residueName_.empty())
    {
        residueName_ = "   ";
    }
    chainId_ = codeUtils::RemoveWhiteSpace(line.substr(21 + shift, 1));
    if (chainId_.empty())
    {
        chainId_ = " ";
    }
    int secondShift = codeUtils::GetSizeOfIntInString(line.substr(26 + shift));
    try
    {
        residueSequenceNumber_ = std::stoi(codeUtils::RemoveWhiteSpace(line.substr(22 + shift, 4 + secondShift)));
    }
    catch (...)
    {
        gmml::log(__LINE__, __FILE__, gmml::ERR, "Error setting residue number from this line:\n" + line);
        residueSequenceNumber_ = gmml::iNotSet;
    }
    // Insertion code gets shifted right by every overrun in residue number.
    insertionCode_ = codeUtils::RemoveWhiteSpace(line.substr(26 + shift + secondShift, 1));
    if (insertionCode_.empty())
    {
        insertionCode_ = " ";
    }
    // Coordinates etc don't get shifted by first overrun in residue number, but do by the rest.
    if (secondShift > 1)
    {
        shift += (secondShift - 1); // Combine the shifts, but ignore the first shift in residue sequence number from here on
    }
    try
    {
        this->setCoordinate(GeometryTopology::Coordinate(codeUtils::RemoveWhiteSpace(line.substr(30 + shift, 8)), codeUtils::RemoveWhiteSpace(line.substr(38 + shift, 8)), codeUtils::RemoveWhiteSpace(line.substr(46 + shift, 8))));
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
        gmml::log(__LINE__, __FILE__, gmml::WAR, "Problem converting to occupancy from: " + line.substr(54 + shift, 6));
    }
    try
    {
        temperatureFactor_ = std::stod(codeUtils::RemoveWhiteSpace(line.substr(60 + shift, 6)));
    }
    catch (...)
    {
        gmml::log(__LINE__, __FILE__, gmml::WAR, "Problem converting to temperatureFactor_ from: " + line.substr(60 + shift, 6));
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
void pdbAtom::SetRecordName(const std::string s)
{
    recordName_ = s;
}
void pdbAtom::SetSerialNumber(const int atom_serial_number)
{
    serialNumber_ = atom_serial_number;
}
void pdbAtom::SetAlternateLocation(const std::string atom_alternate_location)
{
    alternateLocation_ = atom_alternate_location;
}
void pdbAtom::SetChainId(const std::string atom_chain_id)
{
    chainId_ = atom_chain_id;
}
void pdbAtom::SetResidueSequenceNumber(const int atom_residue_sequence_number)
{
    residueSequenceNumber_ = atom_residue_sequence_number;
}
void pdbAtom::SetInsertionCode(const std::string atom_insertion_code)
{
    insertionCode_ = atom_insertion_code;
}
void pdbAtom::SetOccupancy(const double atom_occupancy)
{
    occupancy_ = atom_occupancy;
}
void pdbAtom::SetTempretureFactor(const double atom_temperature_factor)
{
    temperatureFactor_ = atom_temperature_factor;
}
void pdbAtom::SetElement(const std::string atom_element_symbol)
{
    element_ = atom_element_symbol;
}
void pdbAtom::SetCharge(const std::string atom_charge)
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
std::string pdbAtom::GetId() const
{
    std::stringstream ss;
    ss << this->getName() << "_" << this->GetSerialNumber();
    return ss.str();
}

std::string pdbAtom::GetId(const std::string &residueId) const
{
    std::stringstream ss;
    ss << this->GetId() << "_" << residueId;
    return ss.str();
}

//////////////////////////////////////////////////////////
//                       DISPLAY FUNCTION               //
//////////////////////////////////////////////////////////
void pdbAtom::Print(std::ostream &out) const
{
    out << "Serial Number: ";
    if(serialNumber_ == gmml::iNotSet)
    {
        out << " ";
    }
    else
    {
        out << serialNumber_;
    }
    out << ", Atom Name: " << this->getName()
            << ", Alternate Location: " << alternateLocation_
            << ", Residue Name: " << residueName_
            << ", Chain ID: " << chainId_
            << ", Residue Sequence Number: ";
    if(residueSequenceNumber_ == gmml::iNotSet)
    {
        out << " ";
    }
    else
    {
        out << residueSequenceNumber_;
    }
    out << ", Inserion Code: " << insertionCode_
            << ", Coordinate: ";
    this->getCoordinate()->Print(out);
    out << ", Occupancy: ";
    if(occupancy_ == gmml::dNotSet)
    {
        out << " ";
    }
    else
    {
        out << occupancy_;
    }
    out << ", Temperature Factor: ";
    if(temperatureFactor_ == gmml::dNotSet)
    {
        out << " ";
    }
    else
    {
        out << temperatureFactor_;
    }
    out << ", Element: " << element_
            << ", Charge: " << charge_ << std::endl;
}

void pdbAtom::Write(std::ostream& stream) const
{ // Just copied from original class.
    stream << std::left << std::setw(6) << this->GetRecordName();
    if(this->GetSerialNumber() != gmml::iNotSet)
        stream << std::right << std::setw(5) << this->GetSerialNumber();
    else
        stream << std::right << std::setw(5) << " ";
    stream << std::left << std::setw(1) << " "
            << std::left << std::setw(4) << this->getName();
    if(this->GetAlternateLocation() == gmml::sNotSet)
        stream << std::left << std::setw(1) << " ";
    else
        stream << std::left << std::setw(1) << this->GetAlternateLocation();
    stream << std::right << std::setw(3) << this->GetResidueName()
                                       << std::left << std::setw(1) << " ";
    if(this->GetChainId() == gmml::sNotSet)
        stream << std::left << std::setw(1) << " ";
    else
        stream << std::left << std::setw(1) << this->GetChainId();
    if(this->GetResidueSequenceNumber() != gmml::iNotSet)
        stream << std::right << std::setw(4) << this->GetResidueSequenceNumber();
    else
        stream << std::right << std::setw(4) << " ";
    if(this->GetInsertionCode() == gmml::sNotSet)
        stream << std::left << std::setw(1) <<  " ";
    else
        stream << std::left << std::setw(1) << this->GetInsertionCode();
    stream << std::left << std::setw(3) << " ";
    if(this->getCoordinate()->CompareTo(GeometryTopology::Coordinate(gmml::dNotSet, gmml::dNotSet, gmml::dNotSet)) == false)
    {
        stream << std::right << std::setw(8) << std::fixed << std::setprecision(3) << this->getCoordinate()->GetX()
                                           << std::right << std::setw(8) << std::fixed << std::setprecision(3) << this->getCoordinate()->GetY()
                                           << std::right << std::setw(8) << std::fixed << std::setprecision(3) << this->getCoordinate()->GetZ();
    }
    else
    {
        stream << std::right << std::setw(8) << " "
                << std::right << std::setw(8) << " "
                << std::right << std::setw(8) << " ";
    }
    if(this->GetOccupancy() != gmml::dNotSet)
        stream << std::right << std::setw(6) << std::fixed << std::setprecision(2) << this->GetOccupancy();
    else
        stream << std::right << std::setw(6) << " ";
    if(this->GetTemperatureFactor() != gmml::dNotSet)
        stream << std::right << std::setw(6) << std::fixed << std::setprecision(2) << this->GetTemperatureFactor();
    else
        stream << std::right << std::setw(6) << " ";
    stream << std::left << std::setw(10) << " "
            << std::right << std::setw(2) << this->GetElementSymbol()
            << std::left << std::setw(2) << this->GetCharge()
            << std::endl;
}
