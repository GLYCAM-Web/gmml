#ifndef GEOMETRYTOPOLOGY_HPP
#define GEOMETRYTOPOLOGY_HPP



#include "coordinate.hpp"
#include "../utils.hpp"


namespace GeometryTopology
{

Coordinate get_cartesian_point_from_internal_coords(MolecularModeling::Atom *a, MolecularModeling::Atom *b, MolecularModeling::Atom *c, double theta_Degrees, double phi_Degrees, double distance_Angstrom);

Coordinate get_cartesian_point_from_internal_coords(Coordinate a, Coordinate b, Coordinate c, double theta_Degrees, double phi_Degrees, double distance_Angstrom);

Coordinate subtract_coordinates(Coordinate minuaend, Coordinate subtrahend);

}

#endif // GEOMETRYTOPOLOGY_HPP
