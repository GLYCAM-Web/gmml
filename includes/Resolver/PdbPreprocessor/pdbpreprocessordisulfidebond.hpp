#ifndef PDBPREPROCESSORDISULFIDEBOND_HPP
#define PDBPREPROCESSORDISULFIDEBOND_HPP

#include <string>
#include <iostream>

namespace PdbPreprocessorSpace
{
    class PdbPreprocessorDisulfideBond
    {
        public:
            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            PdbPreprocessorDisulfideBond();

            /*! \fn
              * Constructor to initialize attributes of an entity of a disulfide bond
              * @param chain_id_1 Chain id of one of the CYS residues involved in a disulfide bond
              * @param chain_id_2 Chain id of the other CYS residue invoved in a disulfide bond
              * @param residue_sequence_number_1 Sequence number of the first CYS residue involved in a disulfide bond
              * @param residue_sequence_number_2 Sequence number of the second CYS residue involved in a disulfide bond
              * @param distance The distance between the two sulfor atom in each CYS residue
              * @param is_bonded Option to bond/unbond the two CYS residues in a disulfide bond
              * @param residue_insertion_code_1 Insertion code of the first CYS residue invloved in a disulfide bond
              * @param residue_indersion_code_2 Insertion code of the second CYS residue involved in a disulfide bond
              * @param residue_alternate_location_1 Alternate location of the first CYS residue involved in a disulfide bond
              * @param resiude_alternate_location_2 Alternate location of the second CYS residue involved in a disulfide bond
              * @param sulfur_atom_serial_number_1 Sulfur atom serial number of the first CYS residue involved in the disulfide bond
              * @param sulfur_atom_serial_number_2 Sulfur atom serial number of the second CYS residue involved in the disulfide bond
              */
            PdbPreprocessorDisulfideBond(char chain_id_1, char chain_id_2, int residue_sequence_number_1, int residue_sequence_number_2, double distance, bool is_bonded,
                                         char residue_insertion_code_1, char residue_insertion_code_2, char residue_alternate_location_1, char residue_alternate_location_2,
                                         int sulfur_atom_serial_number_1, int sulfur_atom_serial_number_2);

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
/** \addtogroup Molecular_Data_Structure
               * @{
               */
            /*! \fn
              * An accessor function in order to access to the residue chain Id 1
              * @return residue_chain_id_1_ attribute of the current object of this class
              */
            char GetResidueChainId1();
            /*! \fn
              * An accessor function in order to access to the residue chain Id 2
              * @return residue_chain_id_2_ attribute of the current object of this class
              */
            char GetResidueChainId2();
            /*! \fn
              * An accessor function in order to access to the residue sequence number 1
              * @return residue_sequence_Number_1_ attribute of the current object of this class
              */
            int GetResidueSequenceNumber1();
            /*! \fn
              * An accessor function in order to access to the residue sequence number 2
              * @return residue_sequence_Number_2_ attribute of the current object of this class
              */
            int GetResidueSequenceNumber2();
            /*! \fn
              * An accessor function in order to access to the distance
              * @return distance_ attribute of the current object of this class
              */
            double GetDistance();
            /*! \fn
              * An accessor function in order to access to the is bonded attribute
              * @return is_bonded_ attribute of the current object of this class
              */
            bool GetIsBonded();
            /*! \fn
              * An accessor function in order to access to the residue insertion code attribute of the first residue
              * @return residue_insertion_code_1_ attribute of the current object of this class
              */
            char GetResidueInsertionCode1();
            /*! \fn
              * An accessor function in order to access to the residue insertion code attribute of the second residue
              * @return residue_insertion_code_2_ attribute of the current object of this class
              */
            char GetResidueInsertionCode2();
            /*! \fn
              * An accessor function in order to access to the residue alternate location attribute of the first residue
              * @return residue_alternate_location_1_ attribute of the current object of this class
              */
            char GetResidueAlternateLocation1();
            /*! \fn
              * An accessor function in order to access to the residue alternate location attribute of the second residue
              * @return residue_alternate_location_2_ attribute of the current object of this class
              */
            char GetResidueAlternateLocation2();
            /*! \fn
              * An accessor function in order to access to the sulfur atom serial number attribute of the first residue
              * @return sulfur_atom_serial_number_1_ attribute of the current object of this class
              */
            int GetSulfurAtomSerialNumber1();
            /*! \fn
              * An accessor function in order to access to the sulfur atom serial number attribute of the second residue
              * @return sulfur_atom_serial_number_2_ attribute of the current object of this class
              */
            int GetSulfurAtomSerialNumber2();
/** @}*/
            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
/** \addtogroup Manipulators
               * @{
               */
            /*! \fn
              * A mutator function in order to set the residue chain id 1 of the current object
              * Set the residue_chain_id_1_ attribute of the current pdb preprocessor disulfide bond
              * @param residue_chain_id_1 The residue chain id 1 attribute of the current object
              */
            void SetResidueChainId1(char residue_chain_id_1);
            /*! \fn
              * A mutator function in order to set the residue chain id 2 of the current object
              * Set the residue_chain_id_2_ attribute of the current pdb preprocessor disulfide bond
              * @param residue_chain_id_2 The residue chain id 2 attribute of the current object
              */
            void SetResidueChainId2(char residue_chain_id_2);
            /*! \fn
              * A mutator function in order to set the residue sequence number 1 of the current object
              * Set the residue_sequence_number_1_ attribute of the current pdb preprocessor disulfide bond
              * @param residue_sequence_number_1 The residue sequence number 1 attribute of the current object
              */
            void SetResidueSequenceNumber1(int residue_sequence_number_1);
            /*! \fn
              * A mutator function in order to set the residue sequence number 2 of the current object
              * Set the residue_sequence_number_2_ attribute of the current pdb preprocessor disulfide bond
              * @param residue_sequence_number_2 The residue sequence number 2 attribute of the current object
              */
            void SetResidueSequenceNumber2(int residue_sequence_number_2);
            /*! \fn
              * A mutator function in order to set the distance of the current object
              * Set the distance_ attribute of the current pdb preprocessor disulfide bond
              * @param distance The distance attribute of the current object
              */
            void SetDistance(double distance);
            /*! \fn
              * A mutator function in order to set the is bonded of the current object
              * Set the is_bonded_ attribute of the current pdb preprocessor disulfide bond
              * @param is_bonded The is bonded attribute of the current object
              */
            void SetIsBonded(bool is_bonded);
            /*! \fn
              * A mutator function in order to set the first residue insertion code of the current object
              * Set the residue_insertion_code_1_ attribute of the current pdb preprocessor disulfide bond
              * @param residue_insertion_code_1 The insertion code 1 attribute of the current object
              */
            void SetResidueInsertionCode1(char residue_insertion_code_1);
            /*! \fn
              * A mutator function in order to set the second residue insertion code of the current object
              * Set the residue_insertion_code_2_ attribute of the current pdb preprocessor disulfide bond
              * @param residue_insertion_code_2 The insertion code 2 attribute of the current object
              */
            void SetResidueInsertionCode2(char residue_insertion_code_2);
            /*! \fn
              * A mutator function in order to set the first residue alternate location of the current object
              * Set the residue_alternate_location_1_ attribute of the current pdb preprocessor disulfide bond
              * @param residue_alternate_location_1 The residue alternate location 1 attribute of the current object
              */
            void SetResidueAlternateLocation1(char residue_alternate_location_1);
            /*! \fn
              * A mutator function in order to set the second residue alternate location of the current object
              * Set the residue_alternate_location_2_ attribute of the current pdb preprocessor disulfide bond
              * @param residue_alternate_location_2 The residue alternate location 2 attribute of the current object
              */
            void SetResidueAlternateLocation2(char residue_alternate_location_2);
            /*! \fn
              * A mutator function in order to set the first sulfur atom serial number of the current object
              * Set the sulfur_atom_serial_number_1_ attribute of the current pdb preprocessor disulfide bond
              * @param sulfur_atom_serial_number_1 The first sulfur atom serial number attribute of the current object
              */
            void SetSulfurAtomSerialNumber1(int sulfur_atom_serial_number_1);
            /*! \fn
              * A mutator function in order to set the second sulfur atom serial number of the current object
              * Set the sulfur_atom_serial_number_2_ attribute of the current pdb preprocessor disulfide bond
              * @param sulfur_atom_serial_number_2 The second sulfur atom serial number attribute of the current object
              */
            void SetSulfurAtomSerialNumber2(int sulfur_atom_serial_number_2);
/** @}*/
            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the pdb preprocessor disulfide bond contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cerr);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            char residue_chain_id_1_;                   /*!< Chian id of one of the CYS residues involved in the disulfide bond >*/
            char residue_chain_id_2_;                   /*!< Chain id of the other CYS residue involved in the disulfide bond >*/
            int residue_sequence_number_1_;             /*!< Sequence number of the first CYS residue invloved in the disulfide bond >*/
            int residue_sequence_number_2_;             /*!< Sequence number of the second CYS residue involved in the disulfide bond >*/
            double distance_;                           /*!< The distance between the two sulfor atoms in each CYS residue that are involved in the disulfide bond >*/
            bool is_bonded_;                            /*!< Option to bond/unbond the two CYS residues involved in the disulfide bond >*/
            char residue_insertion_code_1_;             /*!< Insertion code of the first CYS residue involved in the disulfide bond >*/
            char residue_insertion_code_2_;             /*!< Insertion code of the second CYS residue involved in the disulfide bond >*/
            char residue_alternate_location_1_;         /*!< Alternate location of the first CYS residue involved in the disulfide bond >*/
            char residue_alternate_location_2_;         /*!< Alternate location of the second CYS residue involved in the disulfide bond >*/
            int sulfur_atom_serial_number_1_;           /*!< Sulfur atom serial number of the first CYS residue involved in the disulfide bond >*/
            int sulfur_atom_serial_number_2_;           /*!< Sulfur atom serial number of the second CYS residue involved in the disulfide bond >*/

    };
}

#endif // PDBPREPROCESSORDISULFIDEBOND_HPP
