// Created by: Alireza Khatamian
// Modified by: Alireza Khatamian, Delaram Rahbarinia

#ifndef PDBSHEETSTRAND_HPP
#define PDBSHEETSTRAND_HPP

#include <string>
#include <vector>
#include <iostream>

namespace PdbFileSpace
{
    enum PdbSheetStrandSense
    {
        FIRST_STRAND = 0,
        PARALLEL = 1,
        ANTI_PARALLEL = -1
    };

    class PdbSheetStrandResidue;
    class PdbSheetStrand
    {
        public:
            //////////////////////////////////////////////////////////
            //                    TYPE DEFINITION                   //
            //////////////////////////////////////////////////////////
            typedef std::vector<PdbSheetStrandResidue*> SheetStrandResidueVector;

            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            PdbSheetStrand();
            PdbSheetStrand(const SheetStrandResidueVector strand_residues, PdbSheetStrandSense sense, const std::string& current_atom, const std::string& previous_atom);
            PdbSheetStrand(std::string& line);

            //////////////////////////////////////////////////////////
            //                         ACCESSOR                     //
            //////////////////////////////////////////////////////////
            SheetStrandResidueVector GetStrandResidues();
            PdbSheetStrandSense GetSense();
            std::string GetCurrentAtom();
            std::string GetPreviousAtom();

            //////////////////////////////////////////////////////////
            //                          MUTATOR                     //
            //////////////////////////////////////////////////////////
            void SetStrandResidues(const SheetStrandResidueVector strand_residues);
            void AddStrandResidue(PdbSheetStrandResidue* strand_residue);
            void SetSense(PdbSheetStrandSense sense);
            void SetCurrentAtom(const std::string current_atom);
            void SetPreviousAtom(const std::string previous_atom);

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
            SheetStrandResidueVector strand_residues_;
            PdbSheetStrandSense sense_;
            std::string current_atom_;
            std::string previous_atom_;
    };
}

#endif // PDBSHEETSTRAND_HPP
