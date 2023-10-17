#ifndef INCLUDES_CENTRALDATASTRUCTURE_CDSSELECTIONS_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_CDSSELECTIONS_HPP

#include "includes/CentralDataStructure/residue.hpp"
#include "includes/CentralDataStructure/molecule.hpp"
#include <vector>

namespace cdsSelections
{
    using cds::Residue;
    std::vector<Residue*> selectResiduesByType(std::vector<Residue*> inputResidues, cds::ResidueType queryType,
                                               const bool invert = false);
    std::vector<Residue*> selectResiduesByType(std::vector<Residue*> inputResidues,
                                               std::vector<cds::ResidueType> queryTypes, const bool invert = false);
    unsigned int findHighestResidueNumber(std::vector<Residue*> residues);
    Residue* FindNeighborResidueConnectedViaSpecificAtom(Residue* queryResidue, const std::string queryAtomName);
    void FindConnectedResidues(std::vector<Residue*>& visitedList, Residue* current);
    std::vector<Residue*> selectNClosestResidues(std::vector<Residue*> inputResidues, Residue* query,
                                                 unsigned int n = 1);
} // namespace cdsSelections
#endif
