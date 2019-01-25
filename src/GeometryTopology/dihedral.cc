#include "../../includes/GeometryTopology/dihedral.hpp"

using GeometryTopology::Dihedral;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
Dihedral::Dihedral() {}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
GeometryTopology::Coordinate::CoordinateVector Dihedral::GetCoordinates()
{
    return coordinates_;
}
double Dihedral::GetTorsion()
{
    return torsion_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void Dihedral::SetCoordinates(GeometryTopology::Coordinate::CoordinateVector coordinates)
{
    coordinates_.clear();
    for(GeometryTopology::Coordinate::CoordinateVector::iterator it = coordinates.begin(); it != coordinates.end(); it++)
    {
        coordinates_.push_back(*it);
    }
}
void Dihedral::AddCoordinate(Coordinate *coordinate)
{
    coordinates_.push_back(coordinate);
}
void Dihedral::SetTorsion(double torsion)
{
    torsion_ = torsion;
}

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void Dihedral::Print(std::ostream &out)
{
    out << "";
}
