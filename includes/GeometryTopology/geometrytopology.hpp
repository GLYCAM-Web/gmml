#ifndef GEOMETRYTOPOLOGY_HPP
#define GEOMETRYTOPOLOGY_HPP

#include "../MolecularModeling/atom.hpp"
#include "coordinate.hpp"

using MolecularModeling::Atom;
namespace GeometryTopology
{

Coordinate get_cartesian_point_from_internal_coords(Atom *a, Atom *b, Atom *c, double theta_Degrees, double phi_Degrees, double distance_Angstrom);

Coordinate get_cartesian_point_from_internal_coords(Coordinate a, Coordinate b, Coordinate c, double theta_Degrees, double phi_Degrees, double distance_Angstrom);

Coordinate subtract_coordinates(Coordinate minuaend, Coordinate subtrahend);

double calculateDistanceFromPointToLineBetweenTwoPoints(Coordinate queryPoint, Coordinate linePointA, Coordinate linePointB);

Coordinate CreateMissingCoordinateForTetrahedralAtom(Atom *atom, const double distance = 1.0);

Coordinate CreateMissingCoordinateForTetrahedralAtom(Coordinate *centralCoord, CoordinateVector threeNeighbors, const double distance);

double CalculateDihedralAngle(Coordinate* a1, Coordinate* a2, Coordinate* a3, Coordinate* a4, bool returnRadians = false);

void SetDihedralAngle(Coordinate* a1, Coordinate* a2, Coordinate* a3, Coordinate* a4, const double dihedral_angle, std::vector<Coordinate*> movingCoords);

void SetAngle(Atom *a, Atom *b, Atom *c, const double angle);

void SetAngle(Coordinate* a1, Coordinate* a2, Coordinate* a3, const double angle, std::vector<Coordinate*> coordinatesToMove);

void SetDistance(Atom *a, Atom *b);

}

#endif // GEOMETRYTOPOLOGY_HPP

//Eigen:
//Vector3d v(1,2,3);
//Vector3d w(0,1,2);
//cout << "Dot product: " << v.dot(w) << endl;
//cout << "Cross product:\n" << v.cross(w) << endl;
// Eigen::ParametrizedLine calculateParametrizedLineFrom3DPoints(Coordinate originCoord, Coordinate otherCoord);

