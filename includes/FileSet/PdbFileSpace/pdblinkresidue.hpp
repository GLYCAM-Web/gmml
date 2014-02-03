// Created by: Alireza Khatamian
// Modified by: Alireza Khatamian, Delaram Rahbarinia

#ifndef PDBLINKRESIDUE_HPP
#define PDBLINKRESIDUE_HPP

#include <string>
#include <iostream>

namespace PdbFileSpace
{
    class PdbLinkResidue
    {
        public:
            //////////////////////////////////////////////////////////
            //                    TYPE DEFINITION                   //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            PdbLinkResidue();
            PdbLinkResidue(const std::string& atom_name, char alternate_location_indicator, const std::string& residue_name, char residue_chain_identifier,
                           int residue_sequence_number, char residue_insertion_code, int symmetry_operator);

            //////////////////////////////////////////////////////////
            //                         ACCESSOR                     //
            //////////////////////////////////////////////////////////
            std::string GetAtomName();
            char GetAlternateLocationIndicator();
            std::string GetResidueName();
            char GetResidueChainIdentifier();
            int GetResidueSequenceNumber();
            char GetResidueInsertionCode();
            int GetSymmetryOperator();

            //////////////////////////////////////////////////////////
            //                          MUTATOR                     //
            //////////////////////////////////////////////////////////
            void SetAtomName(const std::string atom_name);
            void SetAlternateLocationIndicator(char alternate_location_indicator);
            void SetResidueName(const std::string residue_name);
            void SetResidueChainIdentifier(char residue_chain_identifier);
            void SetResidueSequenceNumber(int residue_sequence_number);
            void SetResidueInsertionCode(char residue_insertion_code);
            void SetSymmetryOperator(int symmetry_operator);

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
            std::string atom_name_;
            char alternate_location_indicator_;
            std::string residue_name_;
            char residue_chain_identifier_;
            int residue_sequence_number_;
            char residue_insertion_code_;
            int symmetry_operator_;
    };
}

#endif // PDBLINKRESIDUE_HPP
