#include "../../includes/GeometryTopology/angle.hpp"

using GeometryTopology::Angle;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
Angle::Angle() {}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
GeometryTopology::Coordinate::CoordinateVector Angle::GetCoordinates()
{
    return coordinates_;
}
double Angle::GetAngle()
{
    return angle_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void Angle::SetCoordinates(GeometryTopology::Coordinate::CoordinateVector coordinates)
{
    coordinates_.clear();
    for(GeometryTopology::Coordinate::CoordinateVector::iterator it = coordinates.begin(); it != coordinates.end(); it++)
    {
        coordinates_.push_back(*it);
    }
}
void Angle::AddCoordinate(Coordinate *coordinate)
{
    coordinates_.push_back(coordinate);
}
void Angle::SetAngle(double angle)
{
    angle_ = angle;
}

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void Angle::Print(std::ostream &out)
{
    out << "";
}
