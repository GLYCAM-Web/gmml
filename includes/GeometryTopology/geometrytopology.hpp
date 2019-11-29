#ifndef GEOMETRYTOPOLOGY_HPP
#define GEOMETRYTOPOLOGY_HPP

#include "../MolecularModeling/atom.hpp"
#include "coordinate.hpp"

namespace GeometryTopology
{

Coordinate get_cartesian_point_from_internal_coords(MolecularModeling::Atom *a, MolecularModeling::Atom *b, MolecularModeling::Atom *c, double theta_Degrees, double phi_Degrees, double distance_Angstrom);

Coordinate get_cartesian_point_from_internal_coords(Coordinate a, Coordinate b, Coordinate c, double theta_Degrees, double phi_Degrees, double distance_Angstrom);

Coordinate subtract_coordinates(Coordinate minuaend, Coordinate subtrahend);

double calculateDistanceFromPointToLineBetweenTwoPoints(Coordinate queryPoint, Coordinate linePointA, Coordinate linePointB);

Coordinate CreateMissingCoordinateForTetrahedralAtom(Coordinate *centralCoord, CoordinateVector threeNeighbors);

}

#endif // GEOMETRYTOPOLOGY_HPP

//Eigen:
//Vector3d v(1,2,3);
//Vector3d w(0,1,2);
//cout << "Dot product: " << v.dot(w) << endl;
//cout << "Cross product:\n" << v.cross(w) << endl;
// Eigen::ParametrizedLine calculateParametrizedLineFrom3DPoints(Coordinate originCoord, Coordinate otherCoord);

