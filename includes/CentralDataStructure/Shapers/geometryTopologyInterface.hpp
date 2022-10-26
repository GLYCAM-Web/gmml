#ifndef INCLUDES_CENTRALDATASTRUCTURE_SHAPERS_GEOMETRYTOPOLOGYINTERFACE_HPP_
#define INCLUDES_CENTRALDATASTRUCTURE_SHAPERS_GEOMETRYTOPOLOGYINTERFACE_HPP_

#include "includes/GeometryTopology/coordinate.hpp"
#include "includes/GeometryTopology/geometrytopology.hpp"
#include "includes/CentralDataStructure/Selections/atomSelections.hpp"

#include <vector>

using GeometryTopology::Coordinate;
namespace GeometryTopology
{
template <class atomT>
void GeometryTopology::SetAngle(atomT *a, atomT *b, atomT *c, const double angle)
{
    std::vector<atomT*> atomsToMove;
    atomsToMove.push_back(b);

    TemplatedSelections::FindConnectedAtoms(atomsToMove, c);
    std::vector<Coordinate*> coordsToMove;
    atomsToMove.erase(atomsToMove.begin()); //feck
    for(auto & atom : atomsToMove)
    {
        coordsToMove.push_back(atom->getCoordinate());
    }
    GeometryTopology::SetAngle(a->getCoordinate(), b->getCoordinate(), c->getCoordinate(), angle, coordsToMove);
    return;
}

} // namespace
#endif



