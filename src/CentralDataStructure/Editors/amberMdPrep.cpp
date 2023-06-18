#include "includes/CentralDataStructure/Editors/amberMdPrep.hpp"

bool amberMdPrep::checkForNonNaturalProteinResidues(std::vector<cds::Residue*> unknownResidues, const cds::Atom* cAtom, pdb::PreprocessorInformation &ppInfo)
{
    for(auto & unknownResidue : unknownResidues)
    {
        if(unknownResidue->FindAtom("N") != nullptr && unknownResidue->FindAtom("N")->isWithinBondingDistance(cAtom))
        {
            ppInfo.nonNaturalProteinResidues_.emplace_back(unknownResidue->getId());
            return true;
        }
    }
    return false;
}


