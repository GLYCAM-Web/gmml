#ifndef INCLUDES_CENTRALDATASTRUCTURE_SELECTIONS_ATOMSELECTIONS_HPP_
#define INCLUDES_CENTRALDATASTRUCTURE_SELECTIONS_ATOMSELECTIONS_HPP_

#include "includes/CentralDataStructure/coordinate.hpp"
#include "includes/CentralDataStructure/atom.hpp"
#include "includes/CentralDataStructure/residue.hpp"
#include <vector>
#include <string>
#include <algorithm> //find
#include "templatedSelections.hpp"

using cds::Atom;
using cds::Residue;

namespace cdsSelections
{
    Atom* getNonCarbonHeavyAtomNumbered(std::vector<Atom*> atoms, const std::string queryNumber);
    void FindConnectedAtoms(std::vector<Atom*>& visitedAtoms, Atom* currentAtom);
    Atom* getNeighborNamed(const Atom* queryAtom, const std::string neighborName);
    Atom* selectNeighborNotInAtomVector(const Atom* atomWithNeighbors, std::vector<Atom*> queryAtoms);
    std::vector<Atom*> findCycleAtoms(cds::Atom* const starterAtom);
    Atom* guessAnomericAtomByForeignNeighbor(const Residue* queryResidue);
    Atom* guessAnomericAtomByInternalNeighbors(const std::vector<cds::Atom*> atoms);
    std::vector<Coordinate*> getCoordinates(std::vector<Atom*> queryAtoms);
    unsigned long int CountInternalHeavyAtomBonds(std::vector<Atom*> queryAtoms);
    std::vector<Atom*> FindHeavyAtoms(std::vector<Atom*> queryAtoms);
    std::vector<std::string> FindNamesOfAtoms(std::vector<Atom*> queryAtoms);
    unsigned long int CountAtomsWithinBondingDistance(const Atom* queryAtom, std::vector<Atom*> otherAtoms);
    std::vector<Atom*> FindAtomsWithinDistance(const Atom* queryAtom, std::vector<Atom*> otherAtoms,
                                               double distance = 1.0);
} // namespace cdsSelections
#endif
