#ifndef INCLUDES_CENTRALDATASTRUCTURE_MEASUREMENTS_MEASUREMENTS_HPP_
#define INCLUDES_CENTRALDATASTRUCTURE_MEASUREMENTS_MEASUREMENTS_HPP_

#include "includes/GeometryTopology/coordinate.hpp"
#include "includes/CentralDataStructure/atom.hpp"

using GeometryTopology::Coordinate;

namespace cds
{
Coordinate calculateGeometricCenter(const std::vector<Coordinate*> coords);
std::vector<Coordinate*> getCoordinatesFromAtoms(std::vector<cds::Atom*> atoms);
} // namespace
#endif /* INCLUDES_CENTRALDATASTRUCTURE_MEASUREMENTS_MEASUREMENTS_HPP_ */
