#ifndef GMML_INCLUDES_INPUTSET_PDBFILE_PDBRESIDUE_HPP
#define GMML_INCLUDES_INPUTSET_PDBFILE_PDBRESIDUE_HPP

#include <string>
#include <iostream>
#include <functional>
#include "includes/InputSet/PdbFile/atomRecord.hpp"

// Oliver Jan 2022
// This class just holds (non-owning) references to atomEntries so that I can more efficiently pass information
// out of the PdbFile class to other parts of GMML in a structure (residue) that is useful to the other parts.
// An equivalent Record to "residue" doesn't directly exist in the PDB format. Rather residues are made from the information in ATOM records. A group of ATOM records with the same "chainID" and "resSeq" (aka residue number) belong in the same "residue".
// This class is owned by CoordinateSection, which also owns the ATOM records. So the refs stored here should never be "dead".
// By ownership I mean responsible for creating and managing lifetime of.
namespace pdb
{
    class PdbResidue
    {
        public:
            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            PdbResidue(AtomRecord* atomRecord);
            PdbResidue(std::vector<AtomRecord*> atomRecords);
            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
            std::string GetId() const;
            const std::string& GetChainId() const;
            const std::string& GetName() const;
            const std::string& GetInsertionCode() const;
            const int& GetSequenceNumber() const;
            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
            void AddAtom(AtomRecord* atomRecord);
            void SetName(const std::string name);
            //////////////////////////////////////////////////////////
            //                       FUNCTIONS                      //
            //////////////////////////////////////////////////////////
            AtomRecord* FindAtom(const std::string selector) const;
            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            void Print(std::ostream& out = std::cerr) const;
        private:
            AtomRecord* GetFirstAtom() const;
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            std::vector<AtomRecord*> atomRecords_; // Residue does not own these. Owned by residues's owner.
    };
}
#endif // GMML_INCLUDES_INPUTSET_PDBFILE_PDBRESIDUE_HPP
