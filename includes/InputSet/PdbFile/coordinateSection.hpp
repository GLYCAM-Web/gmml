#ifndef GMML_INCLUDES_INPUTSET_PDBFILE_COORDINATESECTION_HPP
#define GMML_INCLUDES_INPUTSET_PDBFILE_COORDINATESECTION_HPP

#include <vector>
#include <iostream>
//#include "includes/InputSet/PdbFile/pdbResidue.hpp"
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
        private:
            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
        //    inline pdb::Residue& GetCurrentResidue() {return residues_.front();}
            inline const std::vector<AtomRecord> GetAtomRecords() const {return atomRecords_;}
        //    inline const std::vector<pdb::Residue> GetResidues() const {return residues_;}
            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
       //     inline void CreateNewResidue(AtomRecord* atomRecord) {residues_.emplace_back(atomRecord);}
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
      //      std::vector<pdb::Residue> residues_;              // We organize by residue to mirror the rest of GMML structure.
            std::vector<AtomRecord> atomRecords_;
    };
}
#endif // GMML_INCLUDES_INPUTSET_PDBFILE_COORDINATESECTION_HPP
