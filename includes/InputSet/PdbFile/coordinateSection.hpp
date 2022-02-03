#ifndef GMML_INCLUDES_INPUTSET_PDBFILE_COORDINATESECTION_HPP
#define GMML_INCLUDES_INPUTSET_PDBFILE_COORDINATESECTION_HPP

#include <vector>
#include <iostream>
#include "includes/InputSet/PdbFile/pdbChain.hpp"
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
            //                       FUNCTIONS                      //
            //////////////////////////////////////////////////////////
            std::vector<std::vector<pdb::PdbResidue>> GetProteinChains(); // Exposed for PdbFile. Hmm maybe this shouldn't be a separate class.
            std::vector<pdb::PdbResidue> GetResidues(); // Exposed for PdbFile. Hmm maybe this shouldn't be a separate class.
            std::vector<pdb::PdbResidue> FindResidues(const std::string selector); // Uh oh, private parts are exposed. Wee-ooo wee-ooo.
            AtomRecord* FindAtom(int serialNumber); // Conect records
            void DeleteAtomRecord(AtomRecord* atom);
            void CreateNewAtomRecord(std::string name, GeometryTopology::Coordinate& coord, AtomRecord* sisterAtom);
            AtomRecord* CreateNewAtomRecord(const std::string& atomName, const std::string& residueName, const int& residueSequenceNumber, const GeometryTopology::Coordinate& coord, const std::string& chainId, const int& modelNumber, AtomRecord* previousAtom);
            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            void Print(std::ostream& out = std::cerr) const;
        private:
            std::vector<std::unique_ptr<AtomRecord>>::iterator FindPositionOfAtom(AtomRecord* queryAtom);
            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
            //inline pdb::PdbResidue* GetCurrentResidue() {return &(residues_.front());}
            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
            //inline pdb::Residue* CreateNewResidue(AtomRecord* atomRecord) {residues_.emplace_back(atomRecord);}
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
//            std::vector<pdb::PdbResidue> residues_;              // We organize by residue to mirror the rest of GMML structure.
            std::vector<std::unique_ptr<AtomRecord>> atomRecords_;
    };
}
#endif // GMML_INCLUDES_INPUTSET_PDBFILE_COORDINATESECTION_HPP
