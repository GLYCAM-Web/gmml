// Author: Alireza Khatamian

#ifndef PDBDISULFIDERESIDUE_HPP
#define PDBDISULFIDERESIDUE_HPP

#include <string>

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
            PdbDisulfideResidue();
            PdbDisulfideResidue(const std::string& residue_name, char residue_chain_identifier, int residue_sequence_number, char residue_insertion_code, int symmetry_operator);

            //////////////////////////////////////////////////////////
            //                         ACCESSOR                     //
            //////////////////////////////////////////////////////////
            std::string GetResidueName();
            char GetResidueChainIdentifier();
            char GetResidueInsertionCode();
            int GetResidueSequenceNumber();
            int GetSymmetryOperator();

            //////////////////////////////////////////////////////////
            //                          MUTATOR                     //
            //////////////////////////////////////////////////////////
            void SetResidueName(const std::string residue_name);
            void SetResidueChainIdentifier(char residue_chain_identifier);
            void SetResidueInsertionCode(char residue_insertion_code);
            void SetResidueSequenceNumber(int residue_sequence_number);
            void SetSymmetryOperator(int symmetry_operator);

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
            char residue_chain_identifier_;
            int residue_sequence_number_;
            char residue_insertion_code_;
            int symmetry_operator_;
    };
}

#endif // PDBDISULFIDERESIDUE_HPP
