#include "includes/CentralDataStructure/Selections/atomSelections.hpp"
#include "includes/CodeUtils/logging.hpp"

Atom* cdsSelections::getNonCarbonHeavyAtomNumbered(std::vector<Atom*> atoms, const std::string queryNumber)
{
    for (auto& atom : atoms)
    { // Assumes atom names like C2, O2 or N2. Nothing else should match.
        if (atom->getName().size() > 1)
        {
            const std::string number  = atom->getName().substr(1);
            const std::string element = atom->getName().substr(0, 1); // Only the first character
            if ((number == queryNumber) && (element != "C") && (element != "H"))
            {
                return atom;
            }
        }
    }
    return nullptr;
}

void cdsSelections::FindConnectedAtoms(std::vector<Atom*>& visitedAtoms, Atom* currentAtom)
{
    visitedAtoms.push_back(currentAtom);
    for (auto& neighbor : currentAtom->getNeighbors())
    {
        if (std::find(visitedAtoms.begin(), visitedAtoms.end(), neighbor) == visitedAtoms.end())
        {                                                              // Keep looking if neighbor wasn't yet visited.
            cdsSelections::FindConnectedAtoms(visitedAtoms, neighbor); // recursive function call
        }
    }
    return;
}

// ToDo ensure all nodes have names, do this up there.
Atom* cdsSelections::getNeighborNamed(const Atom* queryAtom, const std::string neighborName)
{
    return codeUtils::findElementWithName(queryAtom->getNeighbors(), neighborName);
}

// Step 1. If the C1 atom has a neighbor that isn't in the queryResidue, return C1.
// Step 2. If the C2 atom has a neighbor that isn't in the queryResidue, return C2.
// Step 3. Panic.
Atom* cdsSelections::guessAnomericAtom(cds::Residue* queryResidue)
{
    std::vector<std::string> usualSuspects = {"C1", "C2"};
    for (auto& suspectName : usualSuspects)
    {
        Atom* potentialAnomer = codeUtils::findElementWithName(queryResidue->getAtoms(), suspectName);
        if (cdsSelections::selectNeighborNotInAtomVector(potentialAnomer, queryResidue->getAtoms()) != nullptr)
        { // If atom has a foreign neighbor.
            return potentialAnomer;
        }
    }
    std::string message = "Did not find a C1 or C2 with a foreign neighbor in residue: " + queryResidue->getStringId() +
                          ", thus no anomeric atom was found.";
    gmml::log(__LINE__, __FILE__, gmml::ERR, message);
    throw std::runtime_error(message);
    return nullptr;
}

Atom* cdsSelections::selectNeighborNotInAtomVector(const Atom* atomWithNeighbors, std::vector<Atom*> queryAtoms)
{
    for (auto& neighbor : atomWithNeighbors->getNeighbors())
    {
        if (std::find(queryAtoms.begin(), queryAtoms.end(), neighbor) == queryAtoms.end())
        {
            return neighbor;
        }
    }
    return nullptr;
}

std::vector<Coordinate*> cdsSelections::getCoordinates(std::vector<Atom*> queryAtoms)
{
    std::vector<Coordinate*> coords;
    for (auto& atom : queryAtoms)
    {
        coords.push_back(atom->getCoordinate());
    }
    return coords;
}

unsigned long int cdsSelections::CountInternalHeavyAtomBonds(std::vector<Atom*> queryAtoms)
{
    // This belong in a Molecule class. Algo:
    // Get heavy atoms in molecule. Check number of neighbors. Check if there are more atoms within distance d than
    // neighbors. Scream if so.
    unsigned long int count       = 0;
    std::vector<Atom*> heavyAtoms = cdsSelections::FindHeavyAtoms(queryAtoms);
    for (auto& atom : heavyAtoms)
    {
        count += cdsSelections::CountAtomsWithinBondingDistance(atom, heavyAtoms);
    }
    return count;
}

std::vector<Atom*> cdsSelections::FindHeavyAtoms(std::vector<Atom*> queryAtoms)
{
    std::vector<Atom*> foundAtoms;
    foundAtoms.reserve(queryAtoms.size());
    std::vector<std::string> heavyList = {"C", "O", "N", "S", "P"};
    for (auto& atom : queryAtoms)
    {
        if (std::find(heavyList.begin(), heavyList.end(), atom->getElement()) != heavyList.end())
        {
            foundAtoms.push_back(atom);
        }
    }
    return foundAtoms;
}

std::vector<std::string> cdsSelections::FindNamesOfAtoms(std::vector<Atom*> queryAtoms)
{
    std::vector<std::string> foundNames;
    foundNames.reserve(queryAtoms.size());
    for (auto& atom : queryAtoms)
    {
        foundNames.push_back(atom->getName());
    }
    return foundNames;
}

unsigned long int cdsSelections::CountAtomsWithinBondingDistance(const Atom* queryAtom, std::vector<Atom*> otherAtoms)
{
    unsigned long int count = 0;
    for (auto& otherAtom : otherAtoms)
    {
        if (queryAtom->isWithinBondingDistance(otherAtom))
        {
            ++count;
        }
    }
    return count;
}

std::vector<Atom*> cdsSelections::FindAtomsWithinDistance(const Atom* queryAtom, std::vector<Atom*> otherAtoms,
                                                          double distance)
{
    std::vector<Atom*> foundAtoms;
    for (auto& otherAtom : otherAtoms)
    {
        if (queryAtom->getCoordinate()->Distance(otherAtom->getCoordinate()) < distance)
        {
            foundAtoms.push_back(otherAtom);
        }
    }
    return foundAtoms;
}
