#ifndef INCLUDES_CENTRALDATASTRUCTURE_CDSFUNCTIONS_CDSFUNCTIONS_HPP_
#define INCLUDES_CENTRALDATASTRUCTURE_CDSFUNCTIONS_CDSFUNCTIONS_HPP_
#include "includes/CentralDataStructure/atom.hpp"
#include "includes/CentralDataStructure/residue.hpp"
#include <vector>

namespace cds
{
    void bondAtomsByDistance(std::vector<cds::Atom*> atoms);
    void bondAtomsByDistanceSerial(std::vector<cds::Atom*> atoms);
    void bondAtomsAndResiduesByDistance(cds::Residue* residueA, cds::Residue* residueB);
    void bondAtomsAndResiduesByDistance(std::vector<cds::Residue*> residues);
} // namespace cds
#endif /* INCLUDES_CENTRALDATASTRUCTURE_CDSFUNCTIONS_CDSFUNCTIONS_HPP_ */
