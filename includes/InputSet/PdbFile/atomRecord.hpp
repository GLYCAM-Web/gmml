// See http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM for an explanation of atom formats in PDB files

#ifndef GMML_INCLUDES_INPUTSET_PDBFILE_ATOMRECORD_HPP
#define GMML_INCLUDES_INPUTSET_PDBFILE_ATOMRECORD_HPP

#include <string>
#include <iostream>

#include "includes/GeometryTopology/coordinate.hpp"

namespace pdb
{
    class AtomRecord
    {
        public:
            //////////////////////////////////////////////////////////
            //                    CONSTRUCTOR                       //
            //////////////////////////////////////////////////////////
            AtomRecord(std::string& line, int modelNumber = 1);
            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
            int GetSerialNumber();
            std::string GetName();
            std::string GetResidueName();
            char GetChain();
            int GetResidueSequenceNumber();
            char GetInsertionCode();
            GeometryTopology::Coordinate GetCoordinate();
            double GetOccupancy();
            double GetTemperatureFactor();
            std::string GetElementSymbol();
            std::string GetCharge();
            int GetModelNumber();
            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
            void SetAtomSerialNumber(int atom_serial_number);
            void SetAtomName(const std::string atom_name);
            void SetAtomAlternateLocation(char atom_alternate_location);
            void SetAtomResidueName(const std::string atom_residue_name);
            void SetAtomChainId(char atom_chain_id);
            void SetAtomResidueSequenceNumber(int atom_residue_sequence_number);
            void SetAtomInsertionCode(char atom_insertion_code);
            void SetAtomOrthogonalCoordinate(GeometryTopology::Coordinate atom_orthogonal_coordinate);
            void SetAtomOccupancy(double atom_occupancy);
            void SetAtomTempretureFactor(double atom_temperature_factor);
            void SetAtomElementSymbol(const std::string atom_element_symbol);
            void SetAtomCharge(const std::string atom_charge);
            void SetAtomCardIndexInResidueSet(std::string atom_card_index_in_residue_sequence);
            void AddAlternateLocation(AtomRecord* alternate_atom);
            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            void Print(std::ostream& out = std::cerr);
        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            int serialNumber_;                    /*!< Atom serial number in a model card of a pdb file */
            std::string name_;                     /*!< Atom name in a single atom record in a model card of a pdb file */
            char alternateLocation_;              /*!< Atom residue name in a single atom record in a model card of a pdb file */
            std::string atom_residue_name_;             /*!< Residue name that the atom is assigned to */
            char atom_chain_id_;                        /*!< Chain id that the atom belongs to */
            int atom_residue_sequence_number_;          /*!< Sequence number of the residue that the atom is assigned to */
            char atom_insertion_code_;                                  /*!< Insertion code for the atom in it belonging residue */
            GeometryTopology::Coordinate atom_orthogonal_coordinate_;   /*!< Atom coordinate */
            double atom_occupancy_;                                     /*!< Atom occupancy */
            double atom_temperature_factor_;                            /*!< Atom temperature factor */
            std::string atom_element_symbol_;                           /*!< Atom element symbol */
            std::string atom_charge_;                                   /*!< Atom charge */
            std::string atom_card_index_in_residue_sequence_;
            std::vector<AtomRecord*> alternate_atom_locations_;              /*!< Alternate atom locations, as a vector of atom cards, as there is
                                                                  more information that may be needed (such as ID (A,B,C, etc), %occupancy, etc) */
    };
}

#endif// GMML_INCLUDES_INPUTSET_PDBFILE_ATOMRECORD_HPP
