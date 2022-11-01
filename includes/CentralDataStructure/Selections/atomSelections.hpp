#ifndef INCLUDES_CENTRALDATASTRUCTURE_SELECTIONS_ATOMSELECTIONS_HPP_
#define INCLUDES_CENTRALDATASTRUCTURE_SELECTIONS_ATOMSELECTIONS_HPP_

#include "includes/CodeUtils/templatedSelections.hpp"
#include <vector>
#include <string>
#include <algorithm> //find

namespace TemplatedSelections
{

template <class atomT>
atomT* getNonCarbonHeavyAtomNumbered(std::vector<atomT*> atoms, const std::string queryNumber)
{
    for (auto &atom : atoms)
    {    // Assumes atom names like C2, O2 or N2. Nothing else should match.
        if (atom->getName().size() > 1)
        {
            const std::string number = atom->getName().substr(1);
            const std::string element = atom->getName().substr(0,1); // Only the first character
            if ( (number == queryNumber) && (element != "C") && (element != "H"))
            {
                return atom;
            }
        }
    }
    return nullptr;
}

template <class atomT>
void FindConnectedAtoms(std::vector<atomT*> &visitedAtoms, atomT* currentAtom)
{
    visitedAtoms.push_back(currentAtom);
    for(auto &neighbor : currentAtom->getNeighbors())
    {
        if(!codeUtils::isElementPresent(visitedAtoms.begin(), visitedAtoms.end(), neighbor))
        {
            TemplatedSelections::FindConnectedAtoms(visitedAtoms, neighbor); // recursive function call
        }
    }
    return;
}

} // namespace




#endif

