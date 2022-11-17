#ifndef INCLUDES_CENTRALDATASTRUCTURE_SHAPERS_GEOMETRYTOPOLOGYINTERFACE_HPP_
#define INCLUDES_CENTRALDATASTRUCTURE_SHAPERS_GEOMETRYTOPOLOGYINTERFACE_HPP_

// This file replaces some of the functions that are in geometrytopology, the ones that take MoleculeModeling classes, and replaces them
// It still calls the Coordinate accepting classes in geometrytopology.
// The idea is to separate geometrytopology into the parts that use coordinates only.
#include "includes/CentralDataStructure/atom.hpp"
#include "includes/GeometryTopology/coordinate.hpp"
#include "includes/GeometryTopology/geometrytopology.hpp"
#include "includes/CentralDataStructure/Selections/atomSelections.hpp"
#include "includes/MolecularMetadata/GLYCAM/bondlengthbytypepair.hpp"
#include <vector>

using GeometryTopology::Coordinate;
namespace GeometryTopology
{
Coordinate CreateMissingCoordinateForTetrahedralAtom(cds::Atom* centralAtom, const double distance);
void FindAtomsToMoveAndSetAngle(cds::Atom* a, cds::Atom* b, cds::Atom* c, const double angle);
void FindAtomsToMoveSetDistance(cds::Atom* a, cds::Atom* b);
} // namespace
#endif
