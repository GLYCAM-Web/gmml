// Created by: Alireza Khatamian
// Modified by: Alireza Khatamian, Delaram Rahbarinia

#ifndef PDBRESIDUEMODIFICATION_HPP
#define PDBRESIDUEMODIFICATION_HPP

#include <string>

namespace PdbFileSpace
{
    class PdbResidueModification
    {
        public:
            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            PdbResidueModification();
            PdbResidueModification(const std::string& id_code, const std::string& residue_name, char chain_identifier, int sequence_number,
                                   char insertion_code, const std::string& standard_residue_name, const std::string& dscr);
            PdbResidueModification(std::stringstream& stream_block);
            //////////////////////////////////////////////////////////
            //                         ACCESSOR                     //
            //////////////////////////////////////////////////////////
            std::string GetIdCode();
            std::string GetResidueName();
            char GetChainIdentifier();
            int GetSequenceNumber();
            char GetInsertionCode();
            std::string GetStandardResidueName();
            std::string GetDscr();

            //////////////////////////////////////////////////////////
            //                          MUTATOR                     //
            //////////////////////////////////////////////////////////
            void SetIdCode(const std::string id_code);
            void SetResidueName(const std::string residue_name);
            void SetChainIdentifier(char chain_identifier);
            void SetSequenceNumber(int sequence_number);
            void SetInsertionCode(char insertion_code);
            void SetStandardResidueName(const std::string standard_residue_name);
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
            std::string id_code_;
            std::string residue_name_;
            char chain_identifier_;
            int sequence_number_;
            char insertion_code_;
            std::string standard_residue_name_;
            std::string dscr_;
    };
}

#endif // PDBRESIDUEMODIFICATION_HPP
