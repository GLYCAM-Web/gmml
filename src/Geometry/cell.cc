#include "../../includes/Geometry/cell.hpp"
#include "../../includes/Geometry/coordinate.hpp"

using namespace Geometry;

//////////////////////////////////////////////////////////
//                       Constructor                    //
//////////////////////////////////////////////////////////
Cell::Cell()
{
    min_corner_ = Coordinate();
    max_corner_ = Coordinate();
    cell_charge_ = 0.0;
    cell_potential_energy_ = 0.0;
}

Cell::Cell(Coordinate min, Coordinate max)
{
    min_corner_ = Coordinate(min.GetX(), min.GetY(), min.GetZ());
    max_corner_ = Coordinate(max.GetX(), max.GetY(), max.GetZ());
    cell_charge_ = 0.0;
    cell_potential_energy_ = 0.0;
}

Cell::Cell(Coordinate min, Coordinate max, double charge, double potential_energy)
{
    min_corner_ = Coordinate(min.GetX(), min.GetY(), min.GetZ());
    max_corner_ = Coordinate(max.GetX(), max.GetY(), max.GetZ());
    cell_charge_ = charge;
    cell_potential_energy_ = potential_energy;
}
