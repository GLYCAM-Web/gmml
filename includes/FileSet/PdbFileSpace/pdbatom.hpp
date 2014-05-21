// Created by: Delaram Rahbarinia
// Modified by: Alireza Khatamian, Delaram Rahbarinia

#ifndef PDBATOM_HPP
#define PDBATOM_HPP

#include <string>
#include <iostream>

#include "../../../includes/Geometry/coordinate.hpp"

namespace PdbFileSpace
{
    class PdbAtom
    {
        public:

            //////////////////////////////////////////////////////////
            //                       Constructor                    //
            //////////////////////////////////////////////////////////            
            /*! \fn
              * Default constructor
              */
            PdbAtom();
            /*! \fn
              * Constructor with required parameters
              * @param line A single line in a pdb file that represents pdb atom in model card
              */
            PdbAtom(std::string& line);
            PdbAtom(int atom_serial_number, std::string atom_name, char atom_alternate_location, std::string residue_name, char chain_id,
                    int residue_sequence_number, char insertion_code, Geometry::Coordinate coordinate, double occupancy, double tempreture_factor,
                    std::string element_symbol, std::string charge);
            PdbAtom(char chain_id, std::string atom_name, std::string residue_name, int residue_sequence_number, char residue_insertion_code, char atom_alternate_location);

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
            /*! \fn
              * An accessor function in order to access to the atom serial number in a pdb atom
              * @return atom_serial_number_ attribute of the current object of this class
              */
            int GetAtomSerialNumber();
            /*! \fn
              * An accessor function in order to access to the atom name in a pdb atom
              * @return atom_name_ attribute of the current object of this class
              */
            std::string GetAtomName();
            /*! \fn
              * An accessor function in order to access to the atom alternate location in a pdb atom
              * @return atom_alternate_location_ attribute of the current object of this class
              */
            char GetAtomAlternateLocation();
            /*! \fn
              * An accessor function in order to access to the atom residue name in a pdb atom
              * @return atom_residue_name_ attribute of the current object of this class
              */
            std::string GetAtomResidueName();
            /*! \fn
              * An accessor function in order to access to the atom chain id in a pdb atom
              * @return atom_chain_id_ attribute of the current object of this class
              */
            char GetAtomChainId();
            /*! \fn
              * An accessor function in order to access to the atom residue sequence number in a pdb atom
              * @return atom_residue_sequence_number_ attribute of the current object of this class
              */
            int GetAtomResidueSequenceNumber();
            /*! \fn
              * An accessor function in order to access to the atom insertion code in a pdb atom
              * @return atom_insertion_code_ attribute of the current object of this class
              */
            char GetAtomInsertionCode();
            /*! \fn
              * An accessor function in order to access to the atom orthogonal coordinate in a pdb atom
              * @return atom_orthogonal_coordinate_ attribute of the current object of this class
              */
            Geometry::Coordinate GetAtomOrthogonalCoordinate();
            /*! \fn
              * An accessor function in order to access to the atom occupancy in a pdb atom
              * @return atom_occupancy_ attribute of the current object of this class
              */
            double GetAtomOccupancy();
            /*! \fn
              * An accessor function in order to access to the atom temperature factor in a pdb atom
              * @return atom_temperature_factor_ attribute of the current object of this class
              */
            double GetAtomTempretureFactor();
            /*! \fn
              * An accessor function in order to access to the atom element symbol in a pdb atom
              * @return atom_element_symbol_ attribute of the current object of this class
              */
            std::string GetAtomElementSymbol();
            /*! \fn
              * An accessor function in order to access to the atom charge in a pdb atom
              * @return atom_charge_ attribute of the current object of this class
              */
            std::string GetAtomCharge();

            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A mutator function in order to set the atom serial number of the current object
              * Set the atom_serial_number_ attribute of the current pdb atom
              * @param atom_serial_number The atom serial number of the current object
              */
            void SetAtomSerialNumber(int atom_serial_number);
            /*! \fn
              * A mutator function in order to set the record name of the current object
              * Set the atom_name_ attribute of the current pdb atom
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
              * Set the atom_residue_name_ attribute of the current pdb atom
              * @param atom_residue_name The atom residue name of the current object
              */
            void SetAtomResidueName(const std::string atom_residue_name);
            /*! \fn
              * A mutator function in order to set the atom chain id of the current object
              * Set the atom_chain_id_ attribute of the current pdb atom
              * @param atom_chain_id The atom chain id of the current object
              */
            void SetAtomChainId(char atom_chain_id);
            /*! \fn
              * A mutator function in order to set the atom residue sequence number of the current object
              * Set the atom_residue_sequence_number_ attribute of the current pdb atom
              * @param atom_residue_sequence_number The atom residue sequence number of the current object
              */
            void SetAtomResidueSequenceNumber(int atom_residue_sequence_number);
            /*! \fn
              * A mutator function in order to set the atom insertion code of the current object
              * Set the atom_insertion_code_ attribute of the current pdb atom
              * @param atom_insertion_code The atom insertion code of the current object
              */
            void SetAtomInsertionCode(char atom_insertion_code);
            /*! \fn
              * A mutator function in order to set the atom orthogonal coordinate of the current object
              * Set the atom_orthogonal_coordinate_ attribute of the current pdb atom
              * @param atom_orthogonal_coordinate The atom orthogonal coordinate of the current object
              */
            void SetAtomOrthogonalCoordinate(Geometry::Coordinate atom_orthogonal_coordinate);
            /*! \fn
              * A mutator function in order to set the atom occupancy of the current object
              * Set the atom_occupancy_ attribute of the current pdb atom
              * @param atom_occupancy The atom occupancy of the current object
              */
            void SetAtomOccupancy(double atom_occupancy);
            /*! \fn
              * A mutator function in order to set the atom temperature factor of the current object
              * Set the atom_temperature_factor_ attribute of the current pdb atom
              * @param atom_temperature_factor The atom temperature factor of the current object
              */
            void SetAtomTempretureFactor(double atom_temperature_factor);
            /*! \fn
              * A mutator function in order to set the atom element symbol of the current object
              * Set the atom_element_symbol_ attribute of the current pdb atom
              * @param atom_element_symbol The atom element symbol of the current object
              */
            void SetAtomElementSymbol(const std::string atom_element_symbol);
            /*! \fn
              * A mutator function in order to set the atom charge of the current object
              * Set the atom_charge_ attribute of the current pdb atom
              * @param atom_charge The atom charge of the current object
              */
            void SetAtomCharge(const std::string atom_charge);

            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            void Print(std::ostream& out = std::cout);



        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            int atom_serial_number_;                    /*!< Atom serial number in a model card of a pdb file */
            std::string atom_name_;                     /*!< Atom name in a single atom record in a model card of a pdb file */
            char atom_alternate_location_;              /*!< Atom residue name in a single atom record in a model card of a pdb file */
            std::string atom_residue_name_;             /*!< Residue name that the atom is assigned to */
            char atom_chain_id_;                        /*!< Chain id that the atom belongs to */
            int atom_residue_sequence_number_;          /*!< Sequence number of the residue that the atom is assigned to */
            char atom_insertion_code_;                  /*!< Insertion code for the atom in it belonging residue */
            Geometry::Coordinate atom_orthogonal_coordinate_;   /*!< Atom coordinate */
            double atom_occupancy_;                             /*!< Atom occupancy */
            double atom_temperature_factor_;                     /*!< Atom temperature factor */
            std::string atom_element_symbol_;                   /*!< Atom element symbol */
            std::string atom_charge_;                           /*!< Atom charge */

    };
}

#endif // PDBATOM_HPP
