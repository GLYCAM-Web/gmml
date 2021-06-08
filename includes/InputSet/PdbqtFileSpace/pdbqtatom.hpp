#ifndef PDBQTATOM_HPP
#define PDBQTATOM_HPP

#include <string>
#include "../../GeometryTopology/coordinate.hpp"

namespace PdbqtFileSpace
{
    class PdbqtAtom
    {
        public:

            //////////////////////////////////////////////////////////
            //                       Constructor                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            PdbqtAtom();
            /*! \fn
              * Constructor with required parameters
              * @param line A single line in a pdb file that represents pdbqt atom in model card
              */
            PdbqtAtom(std::string& line);
            /*! \fn
              * Constructor with required parameters
              * @param atom_serial_number
              * @param atom_name
              * @param atom_alternate_location
              * @param residue_name
              * @param chain_id
              * @param residue_sequence_number
              * @param insertion_code
              * @param coordinate
              * @param occupancy
              * @param temperature_factor
              * @param charge
              * @param atom_type
              */
            PdbqtAtom(int atom_serial_number, std::string atom_name, char atom_alternate_location, std::string residue_name, char chain_id,
                    int residue_sequence_number, char insertion_code, GeometryTopology::Coordinate coordinate, double occupancy, double temperature_factor,
                    double charge, std::string atom_type, std::string type);

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
            /*! \fn
              * An accessor function in order to access to the card type in a pdbqt atom
              * @return card_type_ attribute of the current object of this class
              */
            std::string GetCardType();
            /*! \fn
              * An accessor function in order to access to the atom serial number in a pdbqt atom
              * @return atom_serial_number_ attribute of the current object of this class
              */
            int GetAtomSerialNumber();
            /*! \fn
              * An accessor function in order to access to the atom name in a pdbqt atom
              * @return atom_name_ attribute of the current object of this class
              */
            std::string GetAtomName();
            /*! \fn
              * An accessor function in order to access to the atom alternate location in a pdbqt atom
              * @return atom_alternate_location_ attribute of the current object of this class
              */
            char GetAtomAlternateLocation();
            /*! \fn
              * An accessor function in order to access to the atom residue name in a pdbqt atom
              * @return atom_residue_name_ attribute of the current object of this class
              */
            std::string GetAtomResidueName();
            /*! \fn
              * An accessor function in order to access to the atom chain id in a pdbqt atom
              * @return atom_chain_id_ attribute of the current object of this class
              */
            char GetAtomChainId();
            /*! \fn
              * An accessor function in order to access to the atom residue sequence number in a pdbqt atom
              * @return atom_residue_sequence_number_ attribute of the current object of this class
              */
            int GetAtomResidueSequenceNumber();
            /*! \fn
              * An accessor function in order to access to the atom insertion code in a pdbqt atom
              * @return atom_insertion_code_ attribute of the current object of this class
              */
            char GetAtomInsertionCode();
            /*! \fn
              * An accessor function in order to access to the atom orthogonal coordinate in a pdbqt atom
              * @return atom_orthogonal_coordinate_ attribute of the current object of this class
              */
            GeometryTopology::Coordinate GetAtomOrthogonalCoordinate();
            /*! \fn
              * An accessor function in order to access to the atom occupancy in a pdbqt atom
              * @return atom_occupancy_ attribute of the current object of this class
              */
            double GetAtomOccupancy();
            /*! \fn
              * An accessor function in order to access to the atom temperature factor in a pdbqt atom
              * @return atom_temperature_factor_ attribute of the current object of this class
              */
            double GetAtomTempretureFactor();
            /*! \fn
              * An accessor function in order to access to the atom charge in a pdb atom
              * @return atom_charge_ attribute of the current object of this class
              */
            double GetAtomCharge();
            /*! \fn
              * An accessor function in order to access to the atom type in a pdbqt atom
              * @return atom_type_ attribute of the current object of this class
              */
            std::string GetAtomType();
            /*! \fn
              * An accessor function in order to access to the type of the pdbqt atom
              * @return type_ attribute of the current object of this class
              */
            std::string GetType();

            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A mutator function in order to set the card type of the current object
              * Set the card_type_ attribute of the current pdbqt atom
              * @param card_type The card type of the current object
              */
            void SetCardType(const std::string card_type);
            /*! \fn
              * A mutator function in order to set the atom serial number of the current object
              * Set the atom_serial_number_ attribute of the current pdbqt atom
              * @param atom_serial_number The atom serial number of the current object
              */
            void SetAtomSerialNumber(int atom_serial_number);
            /*! \fn
              * A mutator function in order to set the record name of the current object
              * Set the atom_name_ attribute of the current pdbqt atom
              * @param atom_name The atom name of the current object
              */
            void SetAtomName(const std::string atom_name);
            /*! \fn
              * A mutator function in order to set the atom alternate location of the current object
              * Set the atom_alternate_location_ attribute of the current pdb atom
              * @param atom_alternate_location The atom alternate location of the current object
              */
            void SetAtomAlternateLocation(char atom_alternate_location);
            /*! \fn
              * A mutator function in order to set the atom residue name of the current object
              * Set the atom_residue_name_ attribute of the current pdbqt atom
              * @param atom_residue_name The atom residue name of the current object
              */
            void SetAtomResidueName(const std::string atom_residue_name);
            /*! \fn
              * A mutator function in order to set the atom chain id of the current object
              * Set the atom_chain_id_ attribute of the current pdbqt atom
              * @param atom_chain_id The atom chain id of the current object
              */
            void SetAtomChainId(char atom_chain_id);
            /*! \fn
              * A mutator function in order to set the atom residue sequence number of the current object
              * Set the atom_residue_sequence_number_ attribute of the current pdbqt atom
              * @param atom_residue_sequence_number The atom residue sequence number of the current object
              */
            void SetAtomResidueSequenceNumber(int atom_residue_sequence_number);
            /*! \fn
              * A mutator function in order to set the atom insertion code of the current object
              * Set the atom_insertion_code_ attribute of the current pdbqt atom
              * @param atom_insertion_code The atom insertion code of the current object
              */
            void SetAtomInsertionCode(char atom_insertion_code);
            /*! \fn
              * A mutator function in order to set the atom orthogonal coordinate of the current object
              * Set the atom_orthogonal_coordinate_ attribute of the current pdbqt atom
              * @param atom_orthogonal_coordinate The atom orthogonal coordinate of the current object
              */
            void SetAtomOrthogonalCoordinate(GeometryTopology::Coordinate atom_orthogonal_coordinate);
            /*! \fn
              * A mutator function in order to set the atom occupancy of the current object
              * Set the atom_occupancy_ attribute of the current pdbqt atom
              * @param atom_occupancy The atom occupancy of the current object
              */
            void SetAtomOccupancy(double atom_occupancy);
            /*! \fn
              * A mutator function in order to set the atom temperature factor of the current object
              * Set the atom_temperature_factor_ attribute of the current pdbqt atom
              * @param atom_temperature_factor The atom temperature factor of the current object
              */
            void SetAtomTempretureFactor(double atom_temperature_factor);
            /*! \fn
              * A mutator function in order to set the atom charge of the current object
              * Set the atom_charge_ attribute of the current pdbqt atom
              * @param atom_charge The atom charge of the current object
              */
            void SetAtomCharge(double atom_charge);

            /*! \fn
              * A mutator function in order to set the atom type of the current object
              * Set the atom_type_ attribute of the current pdbqt atom
              * @param atom_type The atom type in a residue set of the current object
              */
            void SetAtomType(std::string atom_type);

            /*! \fn
              * A mutator function in order to set the type of the current object
              * Set the type_ attribute of the current pdbqt atom
              * @param type The type of the current object
              */
            void SetType(std::string type);

            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the pdbqt atom contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cerr);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            std::string card_type_;
            int atom_serial_number_;
            std::string atom_name_;
            char atom_alternate_location_;
            std::string atom_residue_name_;
            char atom_chain_id_;
            int atom_residue_sequence_number_;
            char atom_insertion_code_;
            GeometryTopology::Coordinate atom_orthogonal_coordinate_;
            double atom_occupancy_;
            double atom_temperature_factor_;
            double atom_charge_;
            std::string atom_type_;
            std::string type_;

    };
}

#endif // PDBQTATOM_HPP
