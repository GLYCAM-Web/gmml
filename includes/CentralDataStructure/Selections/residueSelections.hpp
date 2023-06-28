#ifndef INCLUDES_CENTRALDATASTRUCTURE_CDSSELECTIONS_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_CDSSELECTIONS_HPP

#include "includes/CentralDataStructure/residue.hpp"
#include "includes/CentralDataStructure/molecule.hpp"
#include <vector>

namespace cdsSelections
{
    using cds::Residue;
    std::vector<Residue*> selectResiduesByType(std::vector<Residue*> inputResidues, cds::ResidueType queryType);
    std::vector<Residue*> selectResiduesByType(std::vector<Residue*> inputResidues,
                                               std::vector<cds::ResidueType> queryTypes);
    unsigned int findHighestResidueNumber(std::vector<Residue*> residues);
    Residue* FindNeighborResidueConnectedViaSpecificAtom(Residue* queryResidue, const std::string queryAtomName);
} // namespace cdsSelections
#endif
