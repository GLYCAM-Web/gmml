#ifndef GMML_INCLUDES_INPUTSET_PDBFILE_RESIDUE_HPP
#define GMML_INCLUDES_INPUTSET_PDBFILE_RESIDUE_HPP

#include <string>
#include <iostream>
#include "includes/InputSet/PdbFile/atomRecord.hpp"

// Oliver Jan 2022
// This class just holds (non-owning) references to atomEntries so that I can more efficiently pass information
// out of the PdbFile class to other parts of GMML in a structure (residue) that is useful to the other parts.
// An equivalent Record to "residue" doesn't directly exist in the PDB format. Rather residues are made from the information in ATOM records. A group of ATOM records with the same "chainID" and "resSeq" (aka residue number) belong in the same "residue".
// This class is owned by CoordinateSection, which also owns the ATOM records. So the refs stored here should never be "dead".
// By ownership I mean responsible for creating and managing lifetime of.
namespace pdb
{
    class Residue
    {
        public:
            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            Residue();
            Residue(std::vector<AtomRecord*> atomRecords);
            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
//            const char& GetResidueChainId() const;
//            const std::string& GetResidueName() const;
//            const int& GetResidueSequenceNumber() const;
//            const char& GetResidueInsertionCode() const;
//            const char& GetResidueAlternateLocation() const;
            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
            //AddAtom;
            //InsertAtomAfter;
            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
           // void Print(std::ostream& out = std::cerr) const;
        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            std::vector<AtomRecord*> atomEntries_; // Residue does not own these. Owned by residues's owner.
    };
}

#endif // GMML_INCLUDES_INPUTSET_PDBFILE_RESIDUE_HPP
