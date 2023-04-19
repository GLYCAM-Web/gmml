#include "includes/CentralDataStructure/cdsFunctions.hpp"

void cds::bondAtomsByDistanceSerial(std::vector<cds::Atom*> atoms)
{
    for (std::vector<cds::Atom*>::iterator it1 = atoms.begin(); it1 != atoms.end(); ++it1)
    {
        for (std::vector<cds::Atom*>::iterator it2 = std::next(it1); it2 != atoms.end(); ++it2)
        {
            atomicBonds::bondAtomsIfClose(*it1, *it2);
        }
    }
}

void cds::bondAtomsAndResiduesByDistance(cds::Residue* residueA, cds::Residue* residueB)
{
    bool residuesAreConnected = false;
    for (auto & atomA : residueA->getAtoms())
    {
        for (auto & atomB : residueB->getAtoms())
        {
            if (atomicBonds::bondAtomsIfClose(atomA, atomB))
            {
                residuesAreConnected = true; // only needs to be true once to connect residues.
            }
        }
    }
    if (residuesAreConnected)
    {
        std::string edgeName = residueA->getStringId() + "--" + residueB->getStringId();
        residueA->addNeighbor(edgeName, residueB);
    }
}

void cds::bondAtomsAndResiduesByDistance(std::vector<cds::Residue*> residues)
{
    for (std::vector<cds::Residue*>::iterator it1 = residues.begin(); it1 != residues.end(); ++it1)
    {   // First bond by distance for atoms within each residue
        cds::Residue* res1 = *it1;
        std::vector<cds::Atom*> res1Atoms = res1->getAtoms();
        cds::bondAtomsByDistanceSerial(res1Atoms);
        // Then for each residue, find other residues within reasonable residue distance.
        for (std::vector<cds::Residue*>::iterator it2 = std::next(it1); it2 != residues.end(); ++it2)
        {
            cds::Residue* res2 = *it2;
            std::vector<cds::Atom*> res2Atoms = res2->getAtoms();
            double residueDistance = res1Atoms.at(0)->calculateDistance(res2Atoms.at(0));
            if (residueDistance < constants::residueDistanceOverlapCutoff * 1.2)
            {
                cds::bondAtomsAndResiduesByDistance(res1, res2);
            }
        }
    }
}





