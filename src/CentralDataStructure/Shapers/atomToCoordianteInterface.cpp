#include "includes/CentralDataStructure/Shapers/atomToCoordinateInterface.hpp"
#include "includes/CentralDataStructure/Shapers/shapers.hpp"
#include "includes/CentralDataStructure/Measurements/measurements.hpp"
#include "includes/CentralDataStructure/Selections/atomSelections.hpp"
#include "includes/MolecularMetadata/GLYCAM/bondlengthbytypepair.hpp"
#include "includes/CodeUtils/logging.hpp"

std::vector<Coordinate*> cds::getCoordinatesFromAtoms(std::vector<cds::Atom*> atoms)
{
    std::vector<Coordinate*> coordinates;
    for(auto & atom : atoms)
    {
        coordinates.push_back(atom->getCoordinate());
    }
    return coordinates;
}

cds::Coordinate cds::CreateMissingCoordinateForTetrahedralAtom(cds::Atom* centralAtom, const double distance)
{
    if(centralAtom->getNeighbors().size() < 1)
    {
        std::stringstream ss;
        ss << "Error in CreateMissingCoordinateForTetrahedralAtom. centralAtom neighbors is " <<
                centralAtom->getNeighbors().size() << " for " << centralAtom->getId();
        gmml::log(__LINE__,__FILE__,gmml::ERR, ss.str());
        throw std::runtime_error(ss.str());
    }
    std::vector<Coordinate*> threeNeighborCoords;
    std::vector<Atom*> neighbors = centralAtom->getNeighbors();
    for (auto &neighbor : neighbors )
    {
//        std::cout << neighbor->getName() << ", " << neighbor->getCoordinate()->ToString() << "\n";
        threeNeighborCoords.push_back(neighbor->getCoordinate());
    }
    return cds::CreateMissingCoordinateForTetrahedralAtom(centralAtom->getCoordinate(), threeNeighborCoords, distance);
}

void cds::FindAtomsToMoveAndSetAngle(cds::Atom* a, cds::Atom* b, cds::Atom* c, const double angle)
{
    std::vector<cds::Atom*> atomsToMove;
    atomsToMove.push_back(b);
    cdsSelections::FindConnectedAtoms(atomsToMove, c);
    atomsToMove.erase(atomsToMove.begin()); // this is expensive
    std::vector<Coordinate*> coordsToMove = cds::getCoordinatesFromAtoms(atomsToMove);
    cds::SetAngle(a->getCoordinate(), b->getCoordinate(), c->getCoordinate(), angle, coordsToMove);
    return;
}

// parentAtom (e.g. O of OME), childAtom (e.g. C1 of Gal1-, S1 of SO3)
void cds::FindAtomsToMoveSetDistance(cds::Atom* parentAtom, cds::Atom* childAtom)
{ // Figure out distance
//    std::cout << "parent is " << parentAtom->getName() << " " << parentAtom->getCoordinate()->ToString() << "\n";
//    std::cout << "child is " << childAtom->getName() << " " << childAtom->getCoordinate()->ToString() << "\n";
    gmml::MolecularMetadata::GLYCAM::BondLengthByTypePairContainer bondLengthByTypePairContainer;
    //std::cout << "Types for parent is " << parentAtom->getType() << " and child is " << childAtom->getType() << "\n";
    double distance = bondLengthByTypePairContainer.GetBondLengthForAtomTypes(parentAtom->getType(), childAtom->getType());
    //std::cout << "distance to new atom sill be: " << distance << "\n";
    // Create an atom c that is will superimpose onto the a atom, bringing b atom with it.
    Coordinate c = cds::CreateMissingCoordinateForTetrahedralAtom(childAtom, distance);
    //std::cout << "New tetraAtom for child is: " << c.ToString() << "\n";
    Coordinate cToParent(parentAtom->getCoordinate()->GetX() - c.GetX(), parentAtom->getCoordinate()->GetY() - c.GetY(), parentAtom->getCoordinate()->GetZ() - c.GetZ());
//    std::cout << "cToParent is " << cToParent.ToString() << "\n";
    // Figure out which atoms will move
    std::vector<cds::Atom*> atomsToRotate;
    atomsToRotate.push_back(parentAtom); // add Parent atom so search doesn't go through it.
    cdsSelections::FindConnectedAtoms(atomsToRotate, childAtom);
    atomsToRotate.erase(atomsToRotate.begin()); // delete the parentAtom
    for(auto & atom : atomsToRotate)
    {
        atom->getCoordinate()->Translate(cToParent.GetX(), cToParent.GetY(), cToParent.GetZ());
//       std::cout << "Moved " << atom->getName() << " to new position:\n";
//        atom->getCoordinate()->Print(std::cout);
//        std::cout << "\n";
    }
    return;
}
