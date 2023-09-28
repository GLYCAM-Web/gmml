#include "includes/CentralDataStructure/cdsFunctions/atomicConnectivity.hpp"
#include "includes/CentralDataStructure/cdsFunctions/bondByDistance.hpp"
#include "includes/CentralDataStructure/Selections/residueSelections.hpp"
#include "includes/MolecularMetadata/proteinBonding.hpp"
#include "includes/CodeUtils/logging.hpp"

void cds::setBondingForResidue(cds::Residue* proteinRes)
{
    for (auto& bondPair : biology::getBackboneBonding())
    {
        cds::Atom* firstAtom  = proteinRes->FindAtom(bondPair.first);
        cds::Atom* secondAtom = proteinRes->FindAtom(bondPair.second);
        if (firstAtom != nullptr && secondAtom != nullptr)
        {
            firstAtom->addBond(secondAtom);
        }
    }
    for (auto& bondPair : biology::getSidechainBonding(proteinRes->getName()))
    {
        cds::Atom* firstAtom  = proteinRes->FindAtom(bondPair.first);
        cds::Atom* secondAtom = proteinRes->FindAtom(bondPair.second);
        if (firstAtom != nullptr && secondAtom != nullptr)
        {
            firstAtom->addBond(secondAtom);
        }
    }
    return;
}

bool cds::autoConnectSuccessiveResidues(cds::Residue* cTermRes, cds::Residue* nTermRes)
{
    cds::Atom* cAtom = cTermRes->FindAtom("C");
    cds::Atom* nAtom = nTermRes->FindAtom("N");
    if (cAtom != nullptr && nAtom != nullptr && cAtom->isWithinBondingDistance(nAtom))
    {
        cAtom->addBond(nAtom);
        return true;
    }
    return false;
}

void cds::setProteinConnectivity(std::vector<cds::Residue*> proteinResidues)
{
    if (proteinResidues.empty())
    {
        return;
    }
    if (proteinResidues.size() == 1)
    {
        cds::setBondingForResidue(proteinResidues.front());
        return;
    }
    cds::Residue* previousRes = proteinResidues.front();
    for (std::vector<cds::Residue*>::iterator it = proteinResidues.begin() + 1; it != proteinResidues.end(); ++it)
    {
        cds::setBondingForResidue(*it); // internal bonding
        if (!cds::autoConnectSuccessiveResidues(previousRes,
                                                *it)) // Automatically bond the N and C atoms of successive residues
        {
            gmml::log(__LINE__, __FILE__, gmml::WAR,
                      "Gap detected between " + previousRes->getStringId() + (*it)->getStringId());
        }
        previousRes = *it;
    }
    return;
}

void cds::setAtomicConnectivity(std::vector<cds::Residue*> residues)
{
    cds::setProteinConnectivity(cdsSelections::selectResiduesByType(residues, ResidueType::Protein));
    bool invertSelection = true;
    cds::bondAtomsAndResiduesByDistance(
        cdsSelections::selectResiduesByType(residues, ResidueType::Protein, invertSelection));
}
