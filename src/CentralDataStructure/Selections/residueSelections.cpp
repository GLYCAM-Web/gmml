#include "includes/CentralDataStructure/Selections/residueSelections.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CentralDataStructure/Selections/atomSelections.hpp"

using cds::Residue;

std::vector<Residue*> cdsSelections::selectResiduesByType(std::vector<Residue*> inputResidues,
                                                          cds::ResidueType queryType, const bool invert)
{ // Quality of life wrapper: calls the below function with one queryType.
    return selectResiduesByType(inputResidues, std::vector<cds::ResidueType> {queryType}, invert);
}

std::vector<Residue*> cdsSelections::selectResiduesByType(std::vector<Residue*> inputResidues,
                                                          std::vector<cds::ResidueType> queryTypes, const bool invert)
{
    std::vector<Residue*> selectedResidues;
    for (auto& residue : inputResidues)
    {
        auto findResult = std::find(queryTypes.begin(), queryTypes.end(), residue->GetType());
        if ((findResult != queryTypes.end() && !invert) || (findResult == queryTypes.end() && invert))
        {
            selectedResidues.push_back(residue);
        }
    }
    return selectedResidues;
}

unsigned int cdsSelections::findHighestResidueNumber(std::vector<Residue*> residues)
{
    unsigned int highest = residues.back()->getNumber(); // Good start.
    for (auto& residue : residues)
    {
        unsigned int resNumInt = residue->getNumber();
        if (resNumInt > highest)
        {
            highest = resNumInt;
        }
    }
    return highest;
}

// Not having Atom know which Residue it is in makes this funky. Make a decision about whether that happens or not.
Residue* cdsSelections::FindNeighborResidueConnectedViaSpecificAtom(Residue* queryResidue,
                                                                    const std::string queryAtomName)
{
    cds::Atom* queryAtom = queryResidue->FindAtom(queryAtomName);
    if (queryAtom == nullptr)
    {
        gmml::log(__LINE__, __FILE__, gmml::WAR,
                  "Warning: An atom named " + queryAtomName + " not found in residue: " + queryResidue->getStringId());
        return nullptr;
    }
    cds::Atom* foreignAtomNeighbor = cdsSelections::selectNeighborNotInAtomVector(queryAtom, queryResidue->getAtoms());
    if (foreignAtomNeighbor == nullptr)
    {
        gmml::log(__LINE__, __FILE__, gmml::WAR,
                  "Warning: Did not find foreign neighbors for an atom named " + queryAtomName +
                      " in residue: " + queryResidue->getStringId());
        return nullptr;
    }
    for (auto& residueNeighbor : queryResidue->getNeighbors())
    {
        if (residueNeighbor->contains(foreignAtomNeighbor))
        {
            return residueNeighbor; // happy path.
        }
    }
    gmml::log(__LINE__, __FILE__, gmml::WAR,
              "Warning: Did not find a neighbor residue connected via " + queryAtomName +
                  " to residue: " + queryResidue->getStringId());
    return nullptr;
}
