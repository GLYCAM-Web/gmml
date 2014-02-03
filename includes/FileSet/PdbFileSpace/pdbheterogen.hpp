// Created by: Alireza Khatamian
// Modified by: Alireza Khatamian, Delaram Rahbarinia

#ifndef PDBHETEROGEN_HPP
#define PDBHETEROGEN_HPP

#include <string>
#include <iostream>

namespace PdbFileSpace
{
    class PdbHeterogen
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
            PdbHeterogen();
            /*! \fn
              * Constructor with required parameters
              * @param heterogen_id
              * @param chain_identifier
              * @param sequence_number
              * @param insertion_code
              * @param number_of_heterogen_atoms
              * @param dscr
              */
            PdbHeterogen(const std::string& heterogen_id, char chain_identifier, int sequence_number,
                         char insertion_code, int number_of_heterogen_atoms, const std::string& dscr);
            PdbHeterogen(std::stringstream& stream_block);

            //////////////////////////////////////////////////////////
            //                         ACCESSOR                     //
            //////////////////////////////////////////////////////////
            /*! \fn
              * An accessor function in order to access to the heterogen id in a heterogen
              * @return heterogen_id_ attribute of the current object of this class
              */
            std::string GetHeterogenId();
            /*! \fn
              * An accessor function in order to access to the chain identifier in a heterogen
              * @return chain_identifier_ attribute of the current object of this class
              */
            char GetChainIdentifier();
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

            //////////////////////////////////////////////////////////
            //                          MUTATOR                     //
            //////////////////////////////////////////////////////////
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
            void SetChainIdentifier(char chain_identifier);
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

            //////////////////////////////////////////////////////////
            //                        FUNCTIONS                     //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                      DISPLAY FUNCTION                //
            //////////////////////////////////////////////////////////
            void Print(std::ostream& out = std::cout);

        private:
            //////////////////////////////////////////////////////////
            //                        ATTRIBUTES                    //
            //////////////////////////////////////////////////////////
            std::string heterogen_id_;
            char chain_identifier_;
            int sequence_number_;
            char insertion_code_;
            int number_of_heterogen_atoms_;
            std::string dscr_;
    };
}
#endif // PDBHETEROGEN_HPP
