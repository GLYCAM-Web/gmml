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
              * @param residue_name
              * @param residue_chain_identifier
              * @param residue_sequence_number
              * @param residue_insertion_code
              * @param symmetry_operator
              */
            PdbDisulfideResidue(const std::string& residue_name, char residue_chain_identifier, int residue_sequence_number, char residue_insertion_code, int symmetry_operator);

            //////////////////////////////////////////////////////////
            //                         ACCESSOR                     //
            //////////////////////////////////////////////////////////
            /*! \fn
              * An accessor function in order to access to the residue name in a disulfide residue
              * @return resdiue_name_ attribute of the current object of this class
              */
            std::string GetResidueName();
            /*! \fn
              * An accessor function in order to access to the residue chain identifier in a disulfide residue
              * @return residue_chain_identifier_ attribute of the current object of this class
              */
            char GetResidueChainIdentifier();
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

            //////////////////////////////////////////////////////////
            //                          MUTATOR                     //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A mutator function in order to set the residue name of the current object
              * Set the residue_name_ attribute of the current disulfide residue
              * @param residue_name The residue name of the current object
              */
            void SetResidueName(const std::string residue_name);
            /*! \fn
              * A mutator function in order to set the residue chain identifier of the current object
              * Set the residue_chain_identifier_ attribute of the current disulfide residue
              * @param residue_chain_identifier The residue chain identifier of the current object
              */
            void SetResidueChainIdentifier(char residue_chain_identifier);
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
            void Print(std::ostream& out = std::cout);

        private:
            //////////////////////////////////////////////////////////
            //                        ATTRIBUTES                    //
            //////////////////////////////////////////////////////////
            std::string residue_name_;
            char residue_chain_identifier_;
            int residue_sequence_number_;
            char residue_insertion_code_;
            int symmetry_operator_;
    };
}

#endif // PDBDISULFIDERESIDUE_HPP
