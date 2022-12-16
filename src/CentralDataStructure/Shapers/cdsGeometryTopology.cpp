#include "includes/CentralDataStructure/Shapers/cdsGeometryTopology.hpp"
#include "includes/CentralDataStructure/Shapers/rotationMatrix.hpp"
#include "includes/utils.hpp"  //gmml::DIST_EPSILON gmml::ConvertDegree2Radian

// Memory fountain:
double** cds::GenerateRotationMatrix(Coordinate* direction, Coordinate* parent, double angle)
{
    double** rotation_matrix = new double*[3];
    for(int i = 0; i < 3; i++)
        rotation_matrix[i] = new double[4];

    direction->Normalize();
    double u = direction->GetX();
    double v = direction->GetY();
    double w = direction->GetZ();

    double a = parent->GetX();
    double b = parent->GetY();
    double c = parent->GetZ();

    double u2 = u*u;
    double v2 = v*v;
    double w2 = w*w;
    double cos_rotation_angle = cos(angle);
    double sin_rotation_angle = sin(angle);

    rotation_matrix[0][3] = a * (v2 + w2) - u * (b * v + c * w) + (u * (b * v + c * w) - a * (v2 + w2)) * cos_rotation_angle + (b * w - c * v) * sin_rotation_angle;
    rotation_matrix[1][3] = b * (u2 + w2) - v * (a * u + c * w) + (v * (a * u + c * w) - b * (u2 + w2)) * cos_rotation_angle + (c * u - a * w) * sin_rotation_angle;
    rotation_matrix[2][3] = c * (u2 + v2) - w * (a * u + b * v) + (w * (a * u + b * v) - c * (u2 + v2)) * cos_rotation_angle + (a * v - b * u) * sin_rotation_angle;

    rotation_matrix[0][0] = u2 + (v2 + w2) * cos_rotation_angle;
    rotation_matrix[0][1] = u * v * (1 - cos_rotation_angle) - w * sin_rotation_angle;
    rotation_matrix[0][2] = u * w * (1 - cos_rotation_angle) + v * sin_rotation_angle;

    rotation_matrix[1][0] = u * v * (1 - cos_rotation_angle) + w * sin_rotation_angle;
    rotation_matrix[1][1] = v2 + (u2 + w2) * cos_rotation_angle;
    rotation_matrix[1][2] = v * w * (1 - cos_rotation_angle) - u * sin_rotation_angle;

    rotation_matrix[2][0] = u * w * (1 - cos_rotation_angle) - v * sin_rotation_angle;
    rotation_matrix[2][1] = v * w * (1 - cos_rotation_angle) + u * sin_rotation_angle;
    rotation_matrix[2][2] = w2 + (u2 + v2) * cos_rotation_angle;


    return rotation_matrix;
}

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
//    double** dihedral_angle_matrix = cds::GenerateRotationMatrix(&b4, a2, current_dihedral - constants::degree2Radian(dihedral_angle));
//    for(auto &coord : movingCoords)
//    {
//        double x =dihedral_angle_matrix[0][0] * coord->GetX() + dihedral_angle_matrix[0][1] * coord->GetY() +
//                dihedral_angle_matrix[0][2] * coord->GetZ() + dihedral_angle_matrix[0][3];
//        double y = dihedral_angle_matrix[1][0] * coord->GetX() + dihedral_angle_matrix[1][1] * coord->GetY() +
//                dihedral_angle_matrix[1][2] * coord->GetZ() + dihedral_angle_matrix[1][3];
//        double z = dihedral_angle_matrix[2][0] * coord->GetX() + dihedral_angle_matrix[2][1] * coord->GetY() +
//                dihedral_angle_matrix[2][2] * coord->GetZ() + dihedral_angle_matrix[2][3];
//        coord->SetX(x);
//        coord->SetY(y);
//        coord->SetZ(z);
//    }
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
    double** rotation_matrix = cds::GenerateRotationMatrix(&direction, a2, rotation_angle);
    for(auto &atom_coordinate : coordinatesToMove)
    {
//        std::cout << "Old coordinate Before setting angle is: ";
//        atom_coordinate->Print(std::cout);
//        std::cout << "\n";
        double x = rotation_matrix[0][0] * atom_coordinate->GetX() + rotation_matrix[0][1] * atom_coordinate->GetY() +
                rotation_matrix[0][2] * atom_coordinate->GetZ() + rotation_matrix[0][3];
        double y = rotation_matrix[1][0] * atom_coordinate->GetX() + rotation_matrix[1][1] * atom_coordinate->GetY() +
                rotation_matrix[1][2] * atom_coordinate->GetZ() + rotation_matrix[1][3];
        double z = rotation_matrix[2][0] * atom_coordinate->GetX() + rotation_matrix[2][1] * atom_coordinate->GetY() +
                rotation_matrix[2][2] * atom_coordinate->GetZ() + rotation_matrix[2][3];
        atom_coordinate->SetX(x);
        atom_coordinate->SetY(y);
        atom_coordinate->SetZ(z);
//        std::cout << "New coordinate after setting angle is: ";
//        atom_coordinate->Print(std::cout);
//        std::cout << "\n";
    }
}
