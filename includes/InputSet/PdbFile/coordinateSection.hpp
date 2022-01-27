#ifndef GMML_INCLUDES_INPUTSET_PDBFILE_COORDINATESECTION_HPP
#define GMML_INCLUDES_INPUTSET_PDBFILE_COORDINATESECTION_HPP

#include <vector>
#include <iostream>
#include "includes/InputSet/PdbFile/pdbResidue.hpp"
#include "includes/InputSet/PdbFile/atomRecord.hpp"

namespace pdb
{
    class CoordinateSection
    {
        public:
            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            CoordinateSection();
            CoordinateSection(std::stringstream& stream_block);
            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                        FUNCTIONS                     //
            //////////////////////////////////////////////////////////
            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            void Print(std::ostream& out = std::cerr) const;
            std::vector<pdb::PdbResidue> FindResidues(const std::string selector); // Uh oh, private parts are exposed. Wee-ooo wee-ooo.
        private:
            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
            //inline pdb::PdbResidue* GetCurrentResidue() {return &(residues_.front());}
            std::vector<pdb::PdbResidue> GetResidues() const;
            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
            //inline pdb::Residue* CreateNewResidue(AtomRecord* atomRecord) {residues_.emplace_back(atomRecord);}
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            std::vector<pdb::PdbResidue> residues_;              // We organize by residue to mirror the rest of GMML structure.
            std::vector<std::unique_ptr<AtomRecord>> atomRecords_;
    };
}
#endif // GMML_INCLUDES_INPUTSET_PDBFILE_COORDINATESECTION_HPP
