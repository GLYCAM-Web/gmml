#include "includes/CentralDataStructure/Shapers/geometryTopologyInterface.hpp"
#include "includes/GeometryTopology/geometrytopology.hpp"
#include "includes/CentralDataStructure/Selections/atomSelections.hpp"
#include "includes/MolecularMetadata/GLYCAM/bondlengthbytypepair.hpp"
#include "includes/CodeUtils/logging.hpp"

using GeometryTopology::Coordinate;

std::vector<Coordinate*> GeometryTopology::getCoordinatesFromAtoms(std::vector<cds::Atom*> atoms)
{
    std::vector<Coordinate*> coordinates;
    for(auto & atom : atoms)
    {
        coordinates.push_back(atom->getCoordinate());
    }
    return coordinates;
}

Coordinate CreateMissingCoordinateForTetrahedralAtom(cds::Atom* centralAtom, const double distance)
{
    if(centralAtom->getNeighbors().size() != 4)
    {
        std::stringstream ss;
        ss << "Error in CreateMissingCoordinateForTetrahedralAtom. centralAtom neighbors = " <<
                centralAtom->getNeighbors().size() << " for " << centralAtom->getId();
        gmml::log(__LINE__,__FILE__,gmml::ERR, ss.str());
        throw std::runtime_error(ss.str());
    }
    std::vector<Coordinate*> threeNeighborCoords;
    for (auto &neighbor : centralAtom->getNeighbors())
    {
        threeNeighborCoords.push_back(neighbor->getCoordinate());
    }
    return GeometryTopology::CreateMissingCoordinateForTetrahedralAtom(centralAtom->getCoordinate(), threeNeighborCoords, distance);
}

void GeometryTopology::FindAtomsToMoveAndSetAngle(cds::Atom* a, cds::Atom* b, cds::Atom* c, const double angle)
{
    std::vector<cds::Atom*> atomsToMove;
    atomsToMove.push_back(b);
    TemplatedSelections::FindConnectedAtoms(atomsToMove, c);
    atomsToMove.erase(atomsToMove.begin()); // this is expensive
    std::vector<Coordinate*> coordsToMove = GeometryTopology::getCoordinatesFromAtoms(atomsToMove);
    GeometryTopology::SetAngle(a->getCoordinate(), b->getCoordinate(), c->getCoordinate(), angle, coordsToMove);
    return;
}

void GeometryTopology::FindAtomsToMoveSetDistance(cds::Atom* a, cds::Atom* b)
{ // Figure out distance
    gmml::MolecularMetadata::GLYCAM::BondLengthByTypePairContainer bondLengthByTypePairContainer;
    double distance = bondLengthByTypePairContainer.GetBondLengthForAtomTypes(a->getType(), b->getType());
    // Figure out position of where the b atom should end up relative to a
    Coordinate c = GeometryTopology::CreateMissingCoordinateForTetrahedralAtom(a, distance);
    // Figure out which atoms will move
    std::vector<cds::Atom*> atomsToRotate;
    atomsToRotate.push_back(a);
    TemplatedSelections::FindConnectedAtoms(atomsToRotate, b);
    atomsToRotate.erase(atomsToRotate.begin());
    for(auto & atom : atomsToRotate)
    {
        atom->getCoordinate()->Translate(c.GetX(), c.GetY(), c.GetZ());
    }
    return;
}
