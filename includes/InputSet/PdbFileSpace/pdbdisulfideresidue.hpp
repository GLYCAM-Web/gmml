// Created by: Alireza Khatamian
// Modified by: Alireza Khatamian, Delaram Rahbarinia

#ifndef PDBDISULFIDERESIDUE_HPP
#define PDBDISULFIDERESIDUE_HPP

#include <string>
#include <iostream>

namespace PdbFileSpace
{
    class PdbDisulfideResidue
    {
        public:
            //////////////////////////////////////////////////////////
            //                    TYPE DEFINITION                   //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            PdbDisulfideResidue();
            /*! \fn
              * Constructor with required parameters
              * @param residue_name Name of residue involving in a disulfide bond
              * @param residue_chain_id Chain id of the residue involving in a disulfide bond
              * @param residue_sequence_number Sequence number of the residue in a disulfide bond
              * @param residue_insertion_code Insertion code of the residue in a disulfide bond
              * @param symmetry_operator Symmetry operator of a disulfide bond
              */
            PdbDisulfideResidue(const std::string& residue_name, char residue_chain_id, int residue_sequence_number, char residue_insertion_code, int symmetry_operator);

            //////////////////////////////////////////////////////////
            //                         ACCESSOR                     //
            //////////////////////////////////////////////////////////
/** \addtogroup Molecular_Data_Structure
              * @{
              */
            /*! \fn
              * An accessor function in order to access to the residue name in a disulfide residue
              * @return resdiue_name_ attribute of the current object of this class
              */
            std::string GetResidueName();
            /*! \fn
              * An accessor function in order to access to the residue chain id in a disulfide residue
              * @return residue_chain_id_ attribute of the current object of this class
              */
            char GetResidueChainId();
            /*! \fn
              * An accessor function in order to access to the residue insertion code in a disulfide residue
              * @return residue_insertion_code_ attribute of the current object of this class
              */
            char GetResidueInsertionCode();
            /*! \fn
              * An accessor function in order to access to the residue sequence number in a disulfide residue
              * @return residue_sequence_nember_ attribute of the current object of this class
              */
            int GetResidueSequenceNumber();
            /*! \fn
              * An accessor function in order to access to the symmetry operator in a disulfide residue
              * @return symmetry operator_ attribute of the current object of this class
              */
            int GetSymmetryOperator();
/** @}*/
            //////////////////////////////////////////////////////////
            //                          MUTATOR                     //
            //////////////////////////////////////////////////////////
/** \addtogroup Manipulators
              * @{
              */
            /*! \fn
              * A mutator function in order to set the residue name of the current object
              * Set the residue_name_ attribute of the current disulfide residue
              * @param residue_name The residue name of the current object
              */
            void SetResidueName(const std::string residue_name);
            /*! \fn
              * A mutator function in order to set the residue chain id of the current object
              * Set the residue_chain_id_ attribute of the current disulfide residue
              * @param residue_chain_id The residue chain id of the current object
              */
            void SetResidueChainId(char residue_chain_id);
            /*! \fn
              * A mutator function in order to set the residue insertion code of the current object
              * Set the residue_insertion_code_ attribute of the current disulfide residue
              * @param residue_insertion_code The residue insertion code of the current object
              */
            void SetResidueInsertionCode(char residue_insertion_code);
            /*! \fn
              * A mutator function in order to set the residue sequence number of the current object
              * Set the residue_sequence_number_ attribute of the current disulfide residue
              * @param residue_sequence_number The residue sequence number of the current object
              */
            void SetResidueSequenceNumber(int residue_sequence_number);
            /*! \fn
              * A mutator function in order to set the symmetry operator of the current object
              * Set the symmetry_operator_ attribute of the current disulfide residue
              * @param symmetry_operator The symmetry operator of the current object
              */
            void SetSymmetryOperator(int symmetry_operator);
 /** @}*/
            //////////////////////////////////////////////////////////
            //                        FUNCTIONS                     //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                      DISPLAY FUNCTION                //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the disulfide residue bond card contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cerr);

        private:
            //////////////////////////////////////////////////////////
            //                        ATTRIBUTES                    //
            //////////////////////////////////////////////////////////
            std::string residue_name_;          /*!< Residue name involving in a disulfide bond >*/
            char residue_chain_id_;             /*!< Chain identifier of the residue involving in a disulfide bond >*/
            int residue_sequence_number_;       /*!< Sequence number of the residue involving in a disulfide bond >*/
            char residue_insertion_code_;       /*!< Insertion code of the residue involving in a disulfide bond >*/
            int symmetry_operator_;             /*!< Symmetry operator of a disulfide bond >*/
    };
}

#endif // PDBDISULFIDERESIDUE_HPP
