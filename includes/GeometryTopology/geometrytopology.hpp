#ifndef GEOMETRY_TOPOLOGY_HPP
#define GEOMETRY_TOPOLOGY_HPP

#include "angle.hpp"
#include "cell.hpp"
#include "coordinate.hpp"
#include "dihedral.hpp"
#include "distance.hpp"
#include "grid.hpp"
#include "plane.hpp"

namespace GeometryTopology
{
Coordinate get_cartesian_point_from_internal_coords(MolecularModeling::Atom *a, MolecularModeling::Atom *b, MolecularModeling::Atom *c,
     double theta_Degrees, double phi_Degrees, double distance_Angstrom);

}

#endif // GEOMETRY_TOPOLOGY_HPP
