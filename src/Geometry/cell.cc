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

//////////////////////////////////////////////////////////
//                           ACCESSOR                   //
//////////////////////////////////////////////////////////
Geometry::Coordinate Cell::GetMinCorner()
{
    return min_corner_;
}

Geometry::Coordinate Cell::GetMaxCorner()
{
    return max_corner_;
}

double Cell::GetCellCharge()
{
    return cell_charge_;
}

double Cell::GetCellPotentialEnergy()
{
    return cell_potential_energy_;
}

double Cell::GetCellLength()
{
    return max_corner_.GetX() - min_corner_.GetX();
}

double Cell::GetCellWidth()
{
    max_corner_.GetY() - min_corner_.GetY();
}

double Cell::GetCellHeight()
{
    max_corner_.GetZ() - min_corner_.GetZ();
}

//////////////////////////////////////////////////////////
//                           MUTATOR                    //
//////////////////////////////////////////////////////////
void Cell::SetMinCorner(Geometry::Coordinate min)
{
    min_corner_.SetX(min.GetX());
    min_corner_.SetY(min.GetY());
    min_corner_.SetZ(min.GetZ());
}

void Cell::SetMaxCorner(Geometry::Coordinate max)
{
    max_corner_.SetX(max.GetX());
    max_corner_.SetY(max.GetY());
    max_corner_.SetZ(max.GetZ());
}

void Cell::SetCellCharge(double charge)
{
    cell_charge_ = charge;
}

void Cell::SetCellPotentialEnergy(double potential_energy)
{
    cell_potential_energy_ = potential_energy;
}

//////////////////////////////////////////////////////////+
//                         FUNCTIONS                    //
//////////////////////////////////////////////////////////
Geometry::Coordinate Cell::GetCellCenter()
{
    return Coordinate((max_corner_.GetX() - min_corner_.GetX()) / 2, (max_corner_.GetY() - min_corner_.GetY()) / 2,
                      (max_corner_.GetZ() - min_corner_.GetZ()) / 2);
}

//////////////////////////////////////////////////////////
//                     DISPLAY FUNCTIONS                //
//////////////////////////////////////////////////////////
