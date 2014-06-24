#include "../../../includes/Geometry/InternalCoordinate/dihedral.hpp"

using namespace std;
using namespace Geometry;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
Dihedral::Dihedral() {}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
Dihedral::CoordinateVector Dihedral::GetCoordinates()
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
void Dihedral::SetCoordinates(CoordinateVector coordinates)
{
    coordinates_ = coordinates;
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
void Dihedral::Print(ostream &out)
{
}





