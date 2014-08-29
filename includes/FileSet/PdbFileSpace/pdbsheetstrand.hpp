// Created by: Alireza Khatamian
// Modified by: Alireza Khatamian, Delaram Rahbarinia

#ifndef PDBSHEETSTRAND_HPP
#define PDBSHEETSTRAND_HPP

#include <string>
#include <vector>
#include <iostream>

namespace PdbFileSpace
{
/*! \enum
  * Pdb sheet strand sense enumerator
  */
    enum PdbSheetStrandSense
    {
        FIRST_STRAND = 0,
        PARALLEL = 1,
        ANTI_PARALLEL = -1,
        UnknownStrand
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
            /*! \fn
              * Default constructor
              */
            PdbSheetStrand();
            /*! \fn
              * Constructor with required parameters
              * @param strand_residues
              * @param sense
              * @param current_atom
              * @param previous_atom
              */
            PdbSheetStrand(const SheetStrandResidueVector strand_residues, PdbSheetStrandSense sense, const std::string& current_atom, const std::string& previous_atom);
            /*! \fn
              * Constructor with required parameters
              * @param line
              */
            PdbSheetStrand(std::string& line);

            //////////////////////////////////////////////////////////
            //                         ACCESSOR                     //
            //////////////////////////////////////////////////////////
            /*! \fn
              * An accessor function in order to access to the strand residues in a sheet strand
              * @return strand_residues_ attribute of the current object of this class
              */
            SheetStrandResidueVector GetStrandResidues();
            /*! \fn
              * An accessor function in order to access to the sense in a sheet strand
              * @return sense_ attribute of the current object of this class
              */
            PdbSheetStrandSense GetSense();
            /*! \fn
              * An accessor function in order to access to the current atom in a sheet strand
              * @return current_atom_ attribute of the current object of this class
              */
            std::string GetCurrentAtom();
            /*! \fn
              * An accessor function in order to access to the previous atom in a sheet strand
              * @return previous_atom_ attribute of the current object of this class
              */
            std::string GetPreviousAtom();

            //////////////////////////////////////////////////////////
            //                          MUTATOR                     //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A mutator function in order to set the strand residues of the current object
              * Set the strand_residues_ attribute of the current sheet strand
              * @param strand_residues The strand residues of the current object
              */
            void SetStrandResidues(const SheetStrandResidueVector strand_residues);
            /*! \fn
              * A function in order to add the strand residue to the current object
              * Set the strand_residue_ attribute of the current sheet strand
              * @param strand_residue The strand residue of the current object
              */
            void AddStrandResidue(PdbSheetStrandResidue* strand_residue);
            /*! \fn
              * A mutator function in order to set the sense of the current object
              * Set the sense_ attribute of the current sheet strand
              * @param sense The strand residues of the current object
              */
            void SetSense(PdbSheetStrandSense sense);
            /*! \fn
              * A mutator function in order to set the current atom of the current object
              * Set the current_atom_ attribute of the current sheet strand
              * @param current_atom The current atom of the current object
              */
            void SetCurrentAtom(const std::string current_atom);
            /*! \fn
              * A mutator function in order to set the previous atom of the current object
              * Set the previous_atom_ attribute of the current sheet strand
              * @param previous_atom The previous atom of the current object
              */
            void SetPreviousAtom(const std::string previous_atom);

            //////////////////////////////////////////////////////////
            //                        FUNCTIONS                     //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                      DISPLAY FUNCTION                //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the pdb sheet strand contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
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
