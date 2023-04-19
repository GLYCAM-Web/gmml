#include "includes/CentralDataStructure/Shapers/shapers.hpp"
#include "includes/CentralDataStructure/Shapers/rotationMatrix.hpp"
#include "includes/CentralDataStructure/Measurements/measurements.hpp" // CalculateDihedralAngle
#include "includes/CodeUtils/constants.hpp"

void cds::SetDihedralAngle(Coordinate* a1, Coordinate* a2, Coordinate* a3, Coordinate* a4, const double dihedral_angle, std::vector<Coordinate*>& movingCoords )
{
    double current_dihedral = cds::CalculateDihedralAngle(a1, a2, a3, a4, true);
    Coordinate direction_a2Toa3 = *a3;
    direction_a2Toa3.operator -(*a2);
    direction_a2Toa3.operator *(-1); // No idea.
    RotationMatrix rotationMatrix(&direction_a2Toa3, a2, current_dihedral - constants::degree2Radian(dihedral_angle));
    rotationMatrix.rotateCoordinates(movingCoords);
}

void cds::SetAngle(Coordinate* a1, Coordinate* a2, Coordinate* a3, const double angle, std::vector<Coordinate*> coordinatesToMove)
{
    double current_angle = cds::CalculateAngle(a1, a2, a3, true);
    double rotation_angle = constants::degree2Radian(angle) - current_angle;
    Coordinate direction = *a1;
    direction.operator -(*a2);
    Coordinate b2 = *a3;
    b2.operator -(*a2);
    direction.CrossProduct(b2);
    direction.Normalize();
    RotationMatrix rotationMatrix(&direction, a2, rotation_angle);
    rotationMatrix.rotateCoordinates(coordinatesToMove);
}
