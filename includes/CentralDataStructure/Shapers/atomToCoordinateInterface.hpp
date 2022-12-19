#ifndef INCLUDES_CENTRALDATASTRUCTURE_SHAPERS_ATOMTOCOORDINATEINTERFACE_HPP_
#define INCLUDES_CENTRALDATASTRUCTURE_SHAPERS_ATOMTOCOORDINATEINTERFACE_HPP_

// This file replaces some of the functions that are in geometrytopology, the ones that take MoleculeModeling classes, and replaces them
// It still calls the Coordinate accepting classes in geometrytopology.
// The idea is to separate geometrytopology into the parts that use coordinates only.
#include "includes/CentralDataStructure/atom.hpp"
#include "includes/CentralDataStructure/coordinate.hpp"
#include "includes/CentralDataStructure/Selections/atomSelections.hpp"
#include "includes/MolecularMetadata/GLYCAM/bondlengthbytypepair.hpp"
#include <vector>
#include "shapers.hpp"

using cds::Coordinate;
namespace cds
{
std::vector<Coordinate*> getCoordinatesFromAtoms(std::vector<cds::Atom*> atoms);
Coordinate CreateMissingCoordinateForTetrahedralAtom(cds::Atom* centralAtom, const double distance = 1.0);
void FindAtomsToMoveAndSetAngle(cds::Atom* a, cds::Atom* b, cds::Atom* c, const double angle);
void FindAtomsToMoveSetDistance(cds::Atom* a, cds::Atom* b);
} // namespace
#endif
