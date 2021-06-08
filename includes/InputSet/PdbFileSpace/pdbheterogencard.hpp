// Created by: Alireza Khatamian
// Modified by: Alireza Khatamian, Delaram Rahbarinia

#ifndef PDBHETEROGENCARD_HPP
#define PDBHETEROGENCARD_HPP

#include <string>
#include <iostream>

namespace PdbFileSpace
{
    class PdbHeterogenCard
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
            PdbHeterogenCard();
            /*! \fn
              * Constructor with required parameters
              * @param heterogen_id Heterogen identifier
              * @param chain_identifier Residue chain identifier
              * @param sequence_number Residue sequence number
              * @param insertion_code Residue insertion code
              * @param number_of_heterogen_atoms Number of heterogen atoms in the heterogen
              * @param dscr Short description for a heterogen
              */
            PdbHeterogenCard(const std::string& heterogen_id, char chain_identifier, int sequence_number,
                         char insertion_code, int number_of_heterogen_atoms, const std::string& dscr);
            /*! \fn
              * Constructor with required parameters
              * @param stream_block
              */
            PdbHeterogenCard(std::stringstream& stream_block);

            //////////////////////////////////////////////////////////
            //                         ACCESSOR                     //
            //////////////////////////////////////////////////////////
/** \addtogroup Molecular_Data_Structure
              * @{
              */
            /*! \fn
              * An accessor function in order to access to the heterogen id in a heterogen
              * @return heterogen_id_ attribute of the current object of this class
              */
            std::string GetHeterogenId();
            /*! \fn
              * An accessor function in order to access to the chain identifier in a heterogen
              * @return chain_identifier_ attribute of the current object of this class
              */
            char GetChainId();
            /*! \fn
              * An accessor function in order to access to the sequence number in a heterogen
              * @return sequence number_ attribute of the current object of this class
              */
            int GetSequenceNumber();
            /*! \fn
              * An accessor function in order to access to the insertion code in a heterogen
              * @return insertion_code_ attribute of the current object of this class
              */
            char GetInsertionCode();
            /*! \fn
              * An accessor function in order to access to the number of heterogen atoms in a heterogen
              * @return number_of_heterogen_atoms_ attribute of the current object of this class
              */
            int GetNumberOfHeterogenAtoms();
            /*! \fn
              * An accessor function in order to access to the description in a heterogen
              * @return dscr_ attribute of the current object of this class
              */
            std::string GetDscr();
/** @}*/
            //////////////////////////////////////////////////////////
            //                          MUTATOR                     //
            //////////////////////////////////////////////////////////
/** \addtogroup Manipulators
              * @{
              */
            /*! \fn
              * A mutator function in order to set the heterogen id of the current object
              * Set the heterogen_id_ attribute of the current heterogen
              * @param heterogen_id The heterogen id of the current object
              */
            void SetHeterogenId(const std::string heterogen_id);
            /*! \fn
              * A mutator function in order to set the chain identifier of the current object
              * Set the chain_identifier_ attribute of the current heterogen
              * @param chain_identifier The chain identifier of the current object
              */
            void SetChainId(char chain_identifier);
            /*! \fn
              * A mutator function in order to set the sequence number of the current object
              * Set the sequence_number_ attribute of the current heterogen
              * @param sequence_number The sequence number of the current object
              */
            void SetSequenceNumber(int sequence_number);
            /*! \fn
              * A mutator function in order to set the insertion code of the current object
              * Set the insertion_code_ attribute of the current heterogen
              * @param insertion_code The insertion code of the current object
              */
            void SetInsertionCode(char insertion_code);
            /*! \fn
              * A mutator function in order to set the number of heterogen atoms of the current object
              * Set the number_of_heterogen_atoms_ attribute of the current heterogen
              * @param number_of_heterogen_atoms The number of heterogen atoms of the current object
              */
            void SetNumberOfHeterogenAtoms(int number_of_heterogen_atoms);
            /*! \fn
              * A mutator function in order to set the description of the current object
              * Set the dscr_ attribute of the current heterogen
              * @param dscr The description of the current object
              */
            void SetDscr(const std::string dscr);
/** @}*/
            //////////////////////////////////////////////////////////
            //                        FUNCTIONS                     //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                      DISPLAY FUNCTION                //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the heterogen contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cerr);

        private:
            //////////////////////////////////////////////////////////
            //                        ATTRIBUTES                    //
            //////////////////////////////////////////////////////////
            std::string heterogen_id_;          /*!< Heterogen identifier >*/
            char chain_identifier_;             /*!< Residue chain identifier >*/
            int sequence_number_;               /*!< Residue sequence number >*/
            char insertion_code_;               /*!< Residue insertion code >*/
            int number_of_heterogen_atoms_;     /*!< Number of hetrogen atoms in the current heterogen >*/
            std::string dscr_;                  /*!< Short description for a heterogen >*/
    };
}
#endif // PDBHETEROGENCARD_HPP
