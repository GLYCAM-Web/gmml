// Author: Alireza Khatamian

#ifndef PDBHELIXRESIDUE_HPP
#define PDBHELIXRESIDUE_HPP

#include <string>

namespace PdbFileSpace
{
    class PdbHelixResidue
    {
        public:
            //////////////////////////////////////////////////////////
            //                    TYPE DEFINITION                   //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            PdbHelixResidue();
            PdbHelixResidue(const std::string& residue_name, char residue_chain_id, int residue_sequence_number, char residue_insertion_code);

            //////////////////////////////////////////////////////////
            //                         ACCESSOR                     //
            //////////////////////////////////////////////////////////
            std::string GetResidueName();
            char GetResidueChainId();
            int GetResidueSequenceNumber();
            char GetResidueInsertionCode();

            //////////////////////////////////////////////////////////
            //                          MUTATOR                     //
            //////////////////////////////////////////////////////////
            void SetResidueName(const std::string residue_name);
            void SetResidueChainId(char residue_chain_id);
            void SetResidueSequenceNumber(int residue_sequence_number);
            void SetResidueInsertionCode(char residue_insertion_code);

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
            std::string residue_name_;
            char residue_chain_id_;
            int residue_sequence_number_;
            char residue_insertion_code_;
    };
}
#endif // PDBHELIXRESIDUE_HPP
