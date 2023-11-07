#include "includes/CentralDataStructure/Readers/Pdb/pdbAtom.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbFunctions.hpp"
#include "includes/CodeUtils/constants.hpp" // gmml::iNotSet
#include "includes/CodeUtils/strings.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CentralDataStructure/Writers/pdbWriter.hpp"

using pdb::PdbAtom;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
// pdbAtom::pdbAtom(const std::string& name, const Coordinate& coord)
//{
//	this->setCoordinate(coord);
//	this->setName(name);
// }
PdbAtom::PdbAtom(const std::string& line)
{
    // gmml::log(__LINE__, __FILE__, gmml::INF, "Parsing " + line);
    //  In the PDB file the residue number overruns after 9999 and serial number overruns after 99999. First overun for
    //  serial doesn't matter as there should be a space between the number and the name. So the problem is above 999999
    this->SetRecordName(codeUtils::RemoveWhiteSpace(line.substr(0, 6)));
    // Dealing with number overruns for serialNumber and residueNumber
    int shift = pdb::checkShiftFromSerialNumberOverrun(line);
    try
    {
        this->setNumber(std::stoi(codeUtils::RemoveWhiteSpace(line.substr(6, 6 + shift))));
    }
    catch (...)
    {
        gmml::log(__LINE__, __FILE__, gmml::ERR,
                  "Error converting to atom's serial number from: " + line.substr(6, 6 + shift));
        this->setNumber(constants::iNotSet);
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
    // gmml::log(__LINE__, __FILE__, gmml::INF, "Hi, my name is " + this->getName());
    // alternateLocation_ = ""; // OG Dec2021, we don't want alt locations in gmml for now. Messes with the
    // preprocessor.
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
    int secondShift = pdb::checkSecondShiftFromResidueNumberOverrun(line, shift);
    // int secondShift = codeUtils::GetSizeOfIntInString(line.substr(26 + shift));
    try
    {
        residueSequenceNumber_ = std::stoi(codeUtils::RemoveWhiteSpace(line.substr(22 + shift, 4 + secondShift)));
    }
    catch (...)
    {
        gmml::log(__LINE__, __FILE__, gmml::ERR, "Error setting residue number from this line:\n" + line);
        residueSequenceNumber_ = constants::iNotSet;
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
        shift +=
            (secondShift - 1); // Combine the shifts, but ignore the first shift in residue sequence number from here on
    }
    try
    {
        this->setCoordinate(cds::Coordinate(codeUtils::RemoveWhiteSpace(line.substr(30 + shift, 8)),
                                            codeUtils::RemoveWhiteSpace(line.substr(38 + shift, 8)),
                                            codeUtils::RemoveWhiteSpace(line.substr(46 + shift, 8))));
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
        gmml::log(__LINE__, __FILE__, gmml::WAR, "Problematic line is:" + line);
        occupancy_ = 1.0;
    }
    try
    {
        temperatureFactor_ = std::stod(codeUtils::RemoveWhiteSpace(line.substr(60 + shift, 6)));
    }
    catch (...)
    {
        gmml::log(__LINE__, __FILE__, gmml::WAR,
                  "Problem converting to temperatureFactor_ from: " + line.substr(60 + shift, 6));
        gmml::log(__LINE__, __FILE__, gmml::WAR, "Problematic line is:" + line);
        temperatureFactor_ = 0.0;
    }
    // Worried about case where shift is 2 and these don't exist, so line length is 80. More Ifs: Ugh.
    if (shift <= 2)
    {
        element_ = codeUtils::RemoveWhiteSpace(line.substr(76 + shift, 2));
        ; // this used to be for element: temp[1] = std::tolower(temp[1]);
    }
    if (shift == 0)
    {
        charge_ = codeUtils::RemoveWhiteSpace(line.substr(78, 2));
    }
}

// ToDo Is this necessary? Won't the base class one be called?
PdbAtom::PdbAtom(const std::string& name, const Coordinate& coord) : cds::Atom(name, coord)
{
    this->SetRecordName("ATOM");
    this->SetTempretureFactor(0.0);
    this->SetOccupancy(1.0);
}

/////////////////////////////////////////////////////////
//                       MUTATOR                        //
//////////////////////////////////////////////////////////
void PdbAtom::SetRecordName(const std::string s)
{
    recordName_ = s;
}

void PdbAtom::SetAlternateLocation(const std::string atom_alternate_location)
{
    alternateLocation_ = atom_alternate_location;
}

void PdbAtom::SetChainId(const std::string atom_chain_id)
{
    chainId_ = atom_chain_id;
}

void PdbAtom::SetResidueSequenceNumber(const int atom_residue_sequence_number)
{
    residueSequenceNumber_ = atom_residue_sequence_number;
}

void PdbAtom::SetInsertionCode(const std::string atom_insertion_code)
{
    insertionCode_ = atom_insertion_code;
}

void PdbAtom::SetOccupancy(const double atom_occupancy)
{
    occupancy_ = atom_occupancy;
}

void PdbAtom::SetTempretureFactor(const double atom_temperature_factor)
{
    temperatureFactor_ = atom_temperature_factor;
}

void PdbAtom::SetElement(const std::string atom_element_symbol)
{
    element_ = atom_element_symbol;
}

void PdbAtom::SetCharge(const std::string atom_charge)
{
    charge_ = atom_charge;
}

// void AtomRecord::AddAlternateLocation(AtomRecord* alternate_atom)
//{
//   alternateLocations_.push_back(alternate_atom);
// }
//////////////////////////////////////////////////////////
//                       FUNCTION                       //
//////////////////////////////////////////////////////////
std::string PdbAtom::GetId() const
{
    std::stringstream ss;
    ss << this->getName() << "_" << this->getNumber();
    return ss.str();
}

std::string PdbAtom::GetId(const std::string& residueId) const
{
    std::stringstream ss;
    ss << this->GetId() << "_" << residueId;
    return ss.str();
}

//////////////////////////////////////////////////////////
//                       DISPLAY FUNCTION               //
//////////////////////////////////////////////////////////
void PdbAtom::Print(std::ostream& out) const
{
    out << "Serial Number: ";
    if (this->getNumber() == constants::iNotSet)
    {
        out << " ";
    }
    else
    {
        out << this->getNumber();
    }
    out << ", Atom Name: " << this->getName() << ", Alternate Location: " << alternateLocation_
        << ", Residue Name: " << residueName_ << ", Chain ID: " << chainId_ << ", Residue Sequence Number: ";
    if (residueSequenceNumber_ == constants::iNotSet)
    {
        out << " ";
    }
    else
    {
        out << residueSequenceNumber_;
    }
    out << ", Inserion Code: " << insertionCode_ << ", Coordinate: ";
    this->getCoordinate()->Print(out);
    out << ", Occupancy: ";
    if (occupancy_ == constants::dNotSet)
    {
        out << " ";
    }
    else
    {
        out << occupancy_;
    }
    out << ", Temperature Factor: ";
    if (temperatureFactor_ == constants::dNotSet)
    {
        out << " ";
    }
    else
    {
        out << temperatureFactor_;
    }
    out << ", Element: " << element_ << ", Charge: " << charge_ << std::endl;
}

void PdbAtom::Write(std::ostream& stream, const std::string residueName, const unsigned int residueNumber,
                    const std::string chainId, const std::string insertionCode) const
{
    cds::writeAtomToPdb(stream, this, this->GetRecordName(), residueName, residueNumber, chainId, insertionCode,
                        this->GetOccupancy(), this->GetTemperatureFactor());
}
