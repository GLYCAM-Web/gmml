#include "../../includes/GeometryTopology/plane.hpp"

using GeometryTopology::Plane;

//////////////////////////////////////////////////////////
//                       Constructor                    //
//////////////////////////////////////////////////////////
Plane::Plane() {}

//////////////////////////////////////////////////////////
//                           ACCESSOR                   //
//////////////////////////////////////////////////////////
GeometryTopology::Coordinate Plane::GetV1()
{
    return v1_;
}

GeometryTopology::Coordinate Plane::GetV2()
{
    return v2_;
}

GeometryTopology::Coordinate Plane::GetUnitNormalVector()
{
    GeometryTopology::Coordinate v1 =  GeometryTopology::Coordinate(v1_.GetX(), v1_.GetY(), v1_.GetZ());
    GeometryTopology::Coordinate v2 =  GeometryTopology::Coordinate(v2_.GetX(), v2_.GetY(), v2_.GetZ());
    v1.CrossProduct(v2);
    v1.Normalize();
    return GeometryTopology::Coordinate(v1.GetX(), v1.GetY(), v1.GetZ());
}
//////////////////////////////////////////////////////////
//                           MUTATOR                    //
//////////////////////////////////////////////////////////
void Plane::SetV1(GeometryTopology::Coordinate v1)
{
    v1_.SetX(v1.GetX());
    v1_.SetY(v1.GetY());
    v1_.SetZ(v1.GetZ());
}
void Plane::SetV1(double x, double y, double z)
{
    v1_.SetX(x);
    v1_.SetY(y);
    v1_.SetZ(z);
}
void Plane::SetV2(GeometryTopology::Coordinate v2)
{
    v2_.SetX(v2.GetX());
    v2_.SetY(v2.GetY());
    v2_.SetZ(v2.GetZ());
}
void Plane::SetV2(double x, double y, double z)
{
    v2_.SetX(x);
    v2_.SetY(y);
    v2_.SetZ(z);
}

//////////////////////////////////////////////////////////
//                         FUNCTIONS                    //
//////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////
//                     DISPLAY FUNCTIONS                //
//////////////////////////////////////////////////////////
void Plane::Print(std::ostream& out)
{
    out << "";
}
