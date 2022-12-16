#include "includes/CentralDataStructure/Shapers/cdsGeometryTopology.hpp"
#include "includes/CentralDataStructure/Shapers/rotationMatrix.hpp"
#include "includes/utils.hpp"  //gmml::DIST_EPSILON gmml::ConvertDegree2Radian

void cds::SetDihedralAngle(Coordinate* a1, Coordinate* a2, Coordinate* a3, Coordinate* a4, const double dihedral_angle, std::vector<Coordinate*>& movingCoords )
{
    Coordinate b1 = *a2;
    b1.operator -(*a1);
    Coordinate b2 = *a3;
    b2.operator -(*a2);
    Coordinate b3 = *a4;
    b3.operator -(*a3);
    Coordinate b4 = b2;
    b4.operator *(-1);
    Coordinate b2xb3 = b2;
    b2xb3.CrossProduct(b3);
    Coordinate b1_m_b2n = b1;
    b1_m_b2n.operator *(b2.length());
    Coordinate b1xb2 = b1;
    b1xb2.CrossProduct(b2);
    double current_dihedral = atan2(b1_m_b2n.DotProduct(b2xb3), b1xb2.DotProduct(b2xb3));
    RotationMatrix rotationMatrix(&b4, a2, current_dihedral - constants::degree2Radian(dihedral_angle));
    rotationMatrix.rotateCoordinates(movingCoords);
}

void cds::SetAngle(Coordinate* a1, Coordinate* a2, Coordinate* a3, const double angle, std::vector<Coordinate*> coordinatesToMove)
{
    double current_angle = 0.0;
    Coordinate b1 = *a1;
    b1.operator -(*a2);
    Coordinate b2 = *a3;
    b2.operator -(*a2);
    current_angle = acos((b1.DotProduct(b2)) / (b1.length() * b2.length() + gmml::DIST_EPSILON));
    double rotation_angle = gmml::ConvertDegree2Radian(angle) - current_angle;
    Coordinate direction = b1;
    direction.CrossProduct(b2);
    direction.Normalize();
    RotationMatrix rotationMatrix(&direction, a2, rotation_angle);
    rotationMatrix.rotateCoordinates(coordinatesToMove);
}
