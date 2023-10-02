#ifndef INCLUDES_CENTRALDATASTRUCTURE_EDITORS_AMBERMDPREP_HPP_
#define INCLUDES_CENTRALDATASTRUCTURE_EDITORS_AMBERMDPREP_HPP_
#include "includes/CentralDataStructure/residue.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbPreprocessorInputs.hpp"

// The pdb preprocessing code that lives in pdbFile and pdbModel needs to become free functions that live in here.
// We will want to process things that aren't pdb files. Also the name has changed from preprocessing to AmberMdPrep as
// of ~May 2023
namespace amberMdPrep
{
    bool checkForNonNaturalProteinResidues(std::vector<cds::Residue*> unknownResidues, const cds::Atom* cAtom,
                                           pdb::PreprocessorInformation& ppInfo);
}

#endif /* INCLUDES_CENTRALDATASTRUCTURE_EDITORS_AMBERMDPREP_HPP_ */
