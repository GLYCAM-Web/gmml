#include "includes/CentralDataStructure/Selections/atomSelections.hpp"

cds::Atom* cdsSelections::getNonCarbonHeavyAtomNumbered(std::vector<cds::Atom*> atoms, const std::string queryNumber)
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

void cdsSelections::FindConnectedAtoms(std::vector<cds::Atom*> &visitedAtoms, cds::Atom* currentAtom)
{
    visitedAtoms.push_back(currentAtom);
    for(auto &neighbor : currentAtom->getNeighbors())
    {
        if( std::find(visitedAtoms.begin(), visitedAtoms.end(), neighbor) == visitedAtoms.end())
        { // Keep looking if neighbor wasn't yet visited.
            cdsSelections::FindConnectedAtoms(visitedAtoms, neighbor); // recursive function call
        }
    }
    return;
}

// ToDo ensure all nodes have names, do this up there.

cds::Atom* cdsSelections::getNeighborNamed(const cds::Atom* queryAtom, const std::string neighborName)
{
    return codeUtils::findElementWithName(queryAtom->getNeighbors(), neighborName);
}

