#ifndef INCLUDES_CENTRALDATASTRUCTURE_CDSFUNCTIONS_ATOMICCONNECTIVITY_HPP_
#define INCLUDES_CENTRALDATASTRUCTURE_CDSFUNCTIONS_ATOMICCONNECTIVITY_HPP_

#include "includes/CentralDataStructure/residue.hpp"

namespace cds
{
    void setBondingForResidue(cds::Residue* proteinRes);
    bool autoConnectSuccessiveResidues(cds::Residue* cTermRes, cds::Residue* nTermRes);
    void setProteinConnectivity(std::vector<cds::Residue*> proteinResidues);
    void setAtomicConnectivity(std::vector<cds::Residue*> residues);
} // namespace cds

#endif /* INCLUDES_CENTRALDATASTRUCTURE_CDSFUNCTIONS_ATOMICCONNECTIVITY_HPP_ */
