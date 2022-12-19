#ifndef INCLUDES_CENTRALDATASTRUCTURE_MEASUREMENTS_MEASUREMENTS_HPP_
#define INCLUDES_CENTRALDATASTRUCTURE_MEASUREMENTS_MEASUREMENTS_HPP_

#include "includes/CentralDataStructure/coordinate.hpp"

namespace cds
{

Coordinate calculateGeometricCenter(const std::vector<Coordinate*> coords);
double CalculateDihedralAngle(Coordinate* a1, Coordinate* a2, Coordinate* a3, Coordinate* a4, const bool returnRadians = false);
Coordinate CreateMissingCoordinateForTetrahedralAtom(const Coordinate *centralCoord, std::vector<Coordinate*> threeNeighbors, const double distance = 1.0);
Coordinate get_cartesian_point_from_internal_coords(const Coordinate& a, const Coordinate& b, const Coordinate& c, double angle_Degrees, double dihedral_Degrees, double distance_Angstrom);
Coordinate subtract_coordinates(const Coordinate& minuaend, const Coordinate& subtrahend);

} // namespace
#endif /* INCLUDES_CENTRALDATASTRUCTURE_MEASUREMENTS_MEASUREMENTS_HPP_ */
