// Author: Alireza Khatamian

#ifndef PDBHETEROGEN_HPP
#define PDBHETEROGEN_HPP

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
            PdbHeterogen(const std::string& heterogen_id, const std::string& chain_identifier, int sequence_number,
                         char insertion_code, int number_of_heterogen_atoms, const std::string& dscr);

            //////////////////////////////////////////////////////////
            //                         ACCESSOR                     //
            //////////////////////////////////////////////////////////
            std::string GetHeterogenId();
            std::string GetChainIdentifier();
            int GetSequenceNumber();
            char GetInsertionCode();
            int GetNumberOfHeterogenAtoms();
            std::string GetDscr();

            //////////////////////////////////////////////////////////
            //                          MUTATOR                     //
            //////////////////////////////////////////////////////////
            void SetHeterogenId(const std::string heterogen_id);
            void SetChainIdentifier(const std::string chain_identifier);
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
            std::string chain_identifier_;
            int sequence_number_;
            char insertion_code_;
            int number_of_heterogen_atoms_;
            std::string dscr_;
    };
}
#endif // PDBHETEROGEN_HPP
