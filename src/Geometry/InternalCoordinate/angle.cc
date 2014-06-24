#include "../../../includes/Geometry/InternalCoordinate/angle.hpp"

using namespace std;
using namespace Geometry;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
Angle::Angle() {}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
Angle::CoordinateVector Angle::GetCoordinates()
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
void Angle::SetCoordinates(CoordinateVector coordinates)
{
    coordinates_ = coordinates;
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
void Angle::Print(ostream &out)
{
}




