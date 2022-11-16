#ifndef INCLUDES_CENTRALDATASTRUCTURE_SELECTIONS_ATOMSELECTIONS_HPP_
#define INCLUDES_CENTRALDATASTRUCTURE_SELECTIONS_ATOMSELECTIONS_HPP_

#include "includes/CentralDataStructure/atom.hpp"
#include "includes/CentralDataStructure/residue.hpp"
#include "includes/CodeUtils/templatedSelections.hpp"
#include <vector>
#include <string>
#include <algorithm> //find

namespace cdsSelections
{
cds::Atom* getNonCarbonHeavyAtomNumbered(std::vector<cds::Atom*> atoms, const std::string queryNumber);
void FindConnectedAtoms(std::vector<cds::Atom*> &visitedAtoms, cds::Atom* currentAtom);
cds::Atom* getNeighborNamed(const cds::Atom* queryAtom, const std::string neighborName);
} // namespace

#endif

