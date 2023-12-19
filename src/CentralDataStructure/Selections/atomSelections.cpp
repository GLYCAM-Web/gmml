#include "includes/CentralDataStructure/Selections/atomSelections.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/MolecularModeling/TemplateGraph/Algorithms/include/TotalCycleDecomposition.hpp"

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
Atom* cdsSelections::guessAnomericAtomByForeignNeighbor(const cds::Residue* queryResidue)
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

std::vector<Atom*> cdsSelections::findCycleAtoms(cds::Atom* const starterAtom)
{
    glygraph::Graph<cds::Atom> atomGraph(starterAtom);
    std::vector<std::pair<std::unordered_set<glygraph::Node<Atom>*>, std::unordered_set<glygraph::Edge<Atom>*>>>
        g1Cycles = cycle_decomp::totalCycleDetect(atomGraph);
    std::vector<Atom*> cycleAtoms;
    for (std::pair<std::unordered_set<glygraph::Node<Atom>*>, std::unordered_set<glygraph::Edge<Atom>*>> currCyclePair :
         g1Cycles)
    {
        // std::cout << "Cycle:""\n\tNodes: ";
        for (glygraph::Node<Atom>* currAtom : currCyclePair.first)
        {
            // std::cout << currAtom->getName() + ", ";
            cycleAtoms.push_back(currAtom->getDeriviedClass());
        }
    }
    return cycleAtoms;
}

// For now it's the lowest numbered (e.g. C1 lower than C6) ring atom connected to the ring oxygen.
Atom* cdsSelections::guessAnomericAtomByInternalNeighbors(const std::vector<cds::Atom*> atoms)
{
    std::vector<Atom*> cycleAtoms = cdsSelections::findCycleAtoms(atoms.at(0));
    cds::Atom* cycleOxygen        = nullptr;
    for (auto& cycleAtom : cycleAtoms)
    {
        if (cycleAtom->getDeriviedClass()->getElement() == "O")
        {
            cycleOxygen = cycleAtom->getDeriviedClass();
        }
    }
    if (cycleOxygen == nullptr)
    {
        throw std::runtime_error("Did not find a ring oxygen when trying to guess what the anomeric carbon is. Your "
                                 "atom names, are likely, strange.");
    }
    // std::cout << "\nCycle oxygen is " << cycleOxygen->getName() << std::endl;
    std::vector<Atom*> anomerCandidates = cycleOxygen->getNeighbors();
    // So deciding this isn't straightforward. For me use case of prep files this is fine. For reading form the PDB I
    // want a meeting first to figure out what we really need. Using a lambda to sort the candidates so that C1 appears
    // before C2 etc.
    std::sort(anomerCandidates.begin(), anomerCandidates.end(),
              [](Atom*& a, Atom*& b)
              {
                  return (a->getNumberFromName() < b->getNumberFromName());
              });
    Atom* bestOne = anomerCandidates.front();
    // std::cout << "Anomer is " << bestOne->getName() << "\n";
    return bestOne;
    //    for(auto & anomerCandidate : anomerCandidates)
    //    {
    //        for (auto & candidateNeighbor: anomerCandidate->getNeighbors())
    //        {
    //            if (cdsSelections::selectNeighborNotInAtomVector(candidateNeighbor, cycleAtoms) != nullptr)
    //            { // candidateNeighbor is not in the ring and has type O N or S
    //                std::string element = candidateNeighbor->getElement();
    //                if (element == "S" || element == "O" || element == "N")
    //                {
    //
    //                }
    //            }
    //        }
    //    }
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
