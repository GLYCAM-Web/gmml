// Author: Alireza Khatamian

#ifndef PDBHETEROGEN_HPP
#define PDBHETEROGEN_HPP
// Created by: Alireza Khatamian
// Modified by: Alireza Khatamian, Delaram Rahbarinia

#include <string>

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
            PdbHeterogen();
            PdbHeterogen(const std::string& heterogen_id, char chain_identifier, int sequence_number,
                         char insertion_code, int number_of_heterogen_atoms, const std::string& dscr);
            PdbHeterogen(std::stringstream& stream_block);

            //////////////////////////////////////////////////////////
            //                         ACCESSOR                     //
            //////////////////////////////////////////////////////////
            std::string GetHeterogenId();
            char GetChainIdentifier();
            int GetSequenceNumber();
            char GetInsertionCode();
            int GetNumberOfHeterogenAtoms();
            std::string GetDscr();

            //////////////////////////////////////////////////////////
            //                          MUTATOR                     //
            //////////////////////////////////////////////////////////
            void SetHeterogenId(const std::string heterogen_id);
            void SetChainIdentifier(char chain_identifier);
            void SetSequenceNumber(int sequence_number);
            void SetInsertionCode(char insertion_code);
            void SetNumberOfHeterogenAtoms(int number_of_heterogen_atoms);
            void SetDscr(const std::string dscr);

            //////////////////////////////////////////////////////////
            //                        FUNCTIONS                     //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                      DISPLAY FUNCTION                //
            //////////////////////////////////////////////////////////

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
