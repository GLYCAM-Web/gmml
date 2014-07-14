
#include "../../../includes/Geometry/Graph/edgeattribute.hpp"

using namespace std;
using namespace Geometry;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
EdgeAttribute::EdgeAttribute() {}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
double EdgeAttribute::GetForce()
{
    return force_;
}
double EdgeAttribute::GetWeight()
{
    return weight_;
}
double EdgeAttribute::GetLength()
{
    return length_;
}
gmml::EdgeTypeDescriptor EdgeAttribute::GetEdgeType()
{
    return edge_type_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void EdgeAttribute::SetForce(double force)
{
    force_ = force;
}
void EdgeAttribute::SetWeight(double weight)
{
    weight_ = weight;
}
void EdgeAttribute::SetLength(double length)
{
    length_ = length;
}
void EdgeAttribute::SetEdgeType(gmml::EdgeTypeDescriptor edge_type)
{
    edge_type_ = edge_type;
}

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void EdgeAttribute::Print(ostream &out)
{
}








