#include "../../../includes/Geometry/InternalCoordinate/distance.hpp"

using namespace std;
using namespace Geometry;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
Distance::Distance() {}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
Distance::CoordinateVector Distance::GetCoordinates()
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
void Distance::SetCoordinates(CoordinateVector coordinates)
{
    coordinates_ = coordinates;
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
void Distance::Print(ostream &out)
{
}



