// See http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM for an explanation of atom formats in PDB files
#ifndef GMML_INCLUDES_INPUTSET_PDBFILE_ATOMRECORD_HPP
#define GMML_INCLUDES_INPUTSET_PDBFILE_ATOMRECORD_HPP

#include <string>
#include <iostream>
#include "includes/GeometryTopology/coordinate.hpp"
#include "includes/common.hpp" // gmml::iNotSet
#include "includes/MolecularModeling/TemplateGraph/AbstractObject/includes/Labels.hpp"
#include "includes/CentralDataStructure/cdsAtom.hpp"

namespace pdb
{
class AtomRecord : public cds::cdsAtom
{
public:
    //////////////////////////////////////////////////////////
    //                    CONSTRUCTOR                       //
    //////////////////////////////////////////////////////////
    // Default constructor:
    AtomRecord(const AtomRecord &tempAtom);
    AtomRecord(const std::string& atomName = "", const std::string& residueName = "", const int& residueSequenceNumber = gmml::iNotSet, const std::string& insertionCode = gmml::sNotSet, const GeometryTopology::Coordinate& coord = GeometryTopology::Coordinate(), const std::string& chainId = gmml::sNotSet, const int& modelNumber = 1, const int& serialNumber = gmml::iNotSet, const std::string& recordName = "ATOM",  const std::string& alternateLocation = gmml::sNotSet, const double& occupancy = gmml::dNotSet, const double& temperatureFactor = gmml::dNotSet, const std::string& element = "", const std::string& charge = "");
    // Constructor when reading lines:
    AtomRecord(const std::string& line, int modelNumber = 1);
    // Handy constructor that copies info from sisterAtom.
    AtomRecord(const std::string& name, const GeometryTopology::Coordinate& coord, AtomRecord *sisterAtom);
    //////////////////////////////////////////////////////////
    //                       ACCESSOR                       //
    //////////////////////////////////////////////////////////
    inline const std::string& GetRecordName() const {return recordName_;}
    inline const int& GetSerialNumber() const {return serialNumber_;}
    inline const std::string& GetName() const {return atomName_;}
    inline const std::string& GetResidueName() const {return residueName_;}
    inline const std::string& GetChainId() const {return chainId_;}
    inline const int& GetResidueSequenceNumber() const {return residueSequenceNumber_;}
    inline const std::string& GetInsertionCode() const {return insertionCode_;}
    inline const std::string& GetAlternateLocation() const {return alternateLocation_;}
    inline const GeometryTopology::Coordinate& GetCoordinate() const {return coordinate_;}
    inline const double& GetOccupancy() const {return occupancy_;}
    inline const double& GetTemperatureFactor() const {return temperatureFactor_;}
    inline const std::string& GetElementSymbol() const {return element_;}
    inline const std::string& GetCharge() const {return charge_;}
    inline const int& GetModelNumber() const {return modelNumber_;}
    //////////////////////////////////////////////////////////
    //                       MUTATOR                        //
    //////////////////////////////////////////////////////////
    void SetResidueName(const std::string atom_residue_name); // Make friend of pdb::Residue?
    //////////////////////////////////////////////////////////
    //                       FUNCTION                       //
    //////////////////////////////////////////////////////////
    std::string GetId() const;
    std::string GetResidueId() const;
    double CalculateDistance(const AtomRecord* otherAtom) const;
    //////////////////////////////////////////////////////////
    //                       DISPLAY FUNCTION               //
    //////////////////////////////////////////////////////////
    void Print(std::ostream& out = std::cerr) const;
    void Write(std::ostream& stream) const;
private:
    //////////////////////////////////////////////////////////
    //                       MUTATOR                        //
    //////////////////////////////////////////////////////////
    void SetModelNumber(const int i);
    void SetRecordName(const std::string s);
    void SetSerialNumber(const int atom_serial_number);
    void SetAtomName(const std::string atom_name);
    void SetAlternateLocation(const std::string atom_alternate_location);
    void SetChainId(const std::string atom_chain_id);
    void SetResidueSequenceNumber(const int atom_residue_sequence_number);
    void SetInsertionCode(const std::string atom_insertion_code);
    void SetCoordinate(const GeometryTopology::Coordinate c);
    void SetOccupancy(double atom_occupancy);
    void SetTempretureFactor(const double atom_temperature_factor);
    void SetElement(const std::string atom_element_symbol);
    void SetCharge(const std::string atom_charge);
    //void AddAlternateLocation(AtomRecord* alternate_atom);
    //////////////////////////////////////////////////////////
    //                       ATTRIBUTES                     //
    //////////////////////////////////////////////////////////
    int modelNumber_;                             // Model number that this card belongs to. Note not included in line.
    std::string recordName_;                      // Can be HETATM or ATOM
    int serialNumber_;                            // Atom serial number in a model card of a pdb file
    std::string atomName_;                        // Atom name in a single atom record in a model card of a pdb file
    std::string alternateLocation_;               // Atom residue name in a single atom record in a model card of a pdb file
    std::string residueName_;                     // Residue name that the atom is assigned to
    std::string chainId_;                         // Chain id that the atom belongs to
    int residueSequenceNumber_;                   // Sequence number of the residue that the atom is assigned to
    std::string insertionCode_;                   // Insertion code for the atom in it belonging residue
    double occupancy_;                            // Atom occupancy
    double temperatureFactor_;                    // Atom temperature factor
    std::string element_;                         // Atom element symbol
    std::string charge_;                          // Atom charge
    GeometryTopology::Coordinate coordinate_;     // Atom coordinate
    //std::vector<AtomRecord*> alternateLocations_; // Alternate atom locations, as a vector of atom cards, as there is   more information that may be needed (such as ID (A,B,C, etc), %occupancy, etc)
};
}
#endif// GMML_INCLUDES_INPUTSET_PDBFILE_ATOMRECORD_HPP
