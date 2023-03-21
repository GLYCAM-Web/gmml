#ifndef INCLUDES_CENTRALDATASTRUCTURE_READERS_PDB_PDBSELECTIONS_HPP_
#define INCLUDES_CENTRALDATASTRUCTURE_READERS_PDB_PDBSELECTIONS_HPP_

#include "includes/CentralDataStructure/Readers/Pdb/pdbFile.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbResidue.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbResidueId.hpp"

namespace pdb
{
PdbResidue* residueSelector(const PdbFile& pdbFile, const pdb::ResidueId& residueId, const int modelNumber = 0);
PdbResidue* residueSelector(std::vector<cds::Residue*> residues, const pdb::ResidueId& queryId);
}
#endif /* INCLUDES_CENTRALDATASTRUCTURE_READERS_PDB_PDBSELECTIONS_HPP_ */
