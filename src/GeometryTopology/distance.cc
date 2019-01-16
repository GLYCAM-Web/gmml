#include "../../includes/GeometryTopology/distance.hpp"

using GeometryTopology::Distance;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
Distance::Distance() {}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
GeometryTopology::Coordinate::CoordinateVector Distance::GetCoordinates()
{
    return coordinates_;
}
double Distance::GetLength()
{
    return length_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void Distance::SetCoordinates(GeometryTopology::Coordinate::CoordinateVector coordinates)
{
    coordinates_.clear();
    for(GeometryTopology::Coordinate::CoordinateVector::iterator it = coordinates.begin(); it != coordinates.end(); it++)
    {
        coordinates_.push_back(*it);
    }
}
void Distance::AddCoordinate(Coordinate *coordinate)
{
    coordinates_.push_back(coordinate);
}
void Distance::SetLenght(double length)
{
    length_ = length;
}

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void Distance::Print(std::ostream &out)
{
    out << "";
}
