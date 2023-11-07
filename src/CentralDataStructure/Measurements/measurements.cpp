#include "includes/CentralDataStructure/Measurements/measurements.hpp"
#include <eigen3/Eigen/Geometry>

using cds::Coordinate;

Coordinate cds::calculateGeometricCenter(const std::vector<Coordinate*> coords)
{
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
    for (auto& coord : coords)
    {
        x += coord->GetX();
        y += coord->GetY();
        z += coord->GetZ();
    }
    if (x == 0.0 || y == 0.0 || x == 0.0 || coords.size() == 0)
    {
        throw std::runtime_error("Oliver what were you thinking in cds::calculateGeometricCenter?");
    }
    x = x / coords.size();
    y = y / coords.size();
    z = z / coords.size();
    return Coordinate(x, y, z);
}

double cds::CalculateAngle(Coordinate* a1, Coordinate* a2, Coordinate* a3, const bool returnRadians)
{ // returns Degrees by default, must set bool returnRadians to true for radians.
    double current_angle = 0.0;
    Coordinate b1        = *a1;
    b1.operator-(*a2);
    Coordinate b2 = *a3;
    b2.operator-(*a2);
    current_angle = acos((b1.DotProduct(b2)) / (b1.length() * b2.length() + constants::DIST_EPSILON));
    if (returnRadians)
    {
        return current_angle;
    }
    return (current_angle * (180 / constants::PI_RADIAN)); // Convert to degrees
}

double cds::CalculateDihedralAngle(Coordinate* a1, Coordinate* a2, Coordinate* a3, Coordinate* a4,
                                   const bool returnRadians)
{ // returns Degrees by default, must set bool returnRadians to true for radians.
    Coordinate b1 = *a2;
    b1.operator-(*a1);
    Coordinate b2 = *a3;
    b2.operator-(*a2);
    Coordinate b3 = *a4;
    b3.operator-(*a3);
    Coordinate b4 = b2;
    b4.operator*(-1);

    Coordinate b2xb3 = b2;
    b2xb3.CrossProduct(b3);

    Coordinate b1_m_b2n = b1;
    b1_m_b2n.operator*(b2.length());

    Coordinate b1xb2 = b1;
    b1xb2.CrossProduct(b2);

    double current_dihedral_angle = atan2(b1_m_b2n.DotProduct(b2xb3), b1xb2.DotProduct(b2xb3));
    if (returnRadians)
    {
        return current_dihedral_angle;
    }
    return (current_dihedral_angle * (180 / constants::PI_RADIAN)); // Convert to degrees
}

Coordinate cds::CreateCoordinateForCenterAwayFromNeighbors(const Coordinate* centralCoord,
                                                           std::vector<Coordinate*> threeNeighbors,
                                                           const double distance)
{
    Coordinate combinedVs(0.0, 0.0, 0.0);
    //    double xValue = 0.0;
    //    double yValue = 0.0;
    //    double zValue = 0.0;
    for (auto& neighbor : threeNeighbors)
    {

        //        xValue += centralCoord->GetX() - neighbor->GetX();
        //        yValue += centralCoord->GetY() - neighbor->GetY();
        //        zValue += centralCoord->GetZ() - neighbor->GetZ();
        Coordinate temp(centralCoord->GetX() - neighbor->GetX(), centralCoord->GetY() - neighbor->GetY(),
                        centralCoord->GetZ() - neighbor->GetZ());
        combinedVs += temp;
    }
    //    Coordinate combinedVs(xValue, yValue, zValue);
    combinedVs.Normalize();
    return Coordinate(centralCoord->GetX() + (combinedVs.GetX() * distance),
                      centralCoord->GetY() + (combinedVs.GetY() * distance),
                      centralCoord->GetZ() + (combinedVs.GetZ() * distance));
}

Coordinate cds::calculateCoordinateFromInternalCoords(const Coordinate& a, const Coordinate& b, const Coordinate& c,
                                                      double angle_Degrees, double dihedral_Degrees,
                                                      double distance_Angstrom)
{
    //  std::cout << "Distance: " << distance_Angstrom << std::endl;
    //  std::cout << "Angle: " << angle_Degrees << std::endl;
    //  std::cout << "Dihedral: " << dihedral_Degrees << std::endl;
    double theta_Radians = constants::degree2Radian(angle_Degrees);
    double phi_Radians   = constants::degree2Radian(dihedral_Degrees);
    Coordinate lmn_x, lmn_y, lmn_z;
    double x_p, y_p, z_p;

    // ToDo no. Overload the - operator properly in Coordinate.
    Coordinate cb = cds::subtractCoordinates(b, c); // original
    Coordinate ba = cds::subtractCoordinates(a, b); // original

    lmn_y = ba;
    lmn_y.CrossProduct(cb);
    lmn_y.Normalize();

    lmn_z = cb;
    lmn_z.Normalize();

    lmn_x = lmn_z;
    lmn_x.CrossProduct(lmn_y);

    x_p = distance_Angstrom * sin(theta_Radians) * cos(phi_Radians);
    y_p = distance_Angstrom * sin(theta_Radians) * sin(phi_Radians);
    z_p = distance_Angstrom * cos(theta_Radians);

    double new_x = lmn_x.GetX() * x_p + lmn_y.GetX() * y_p + lmn_z.GetX() * z_p + c.GetX();
    double new_y = lmn_x.GetY() * x_p + lmn_y.GetY() * y_p + lmn_z.GetY() * z_p + c.GetY();
    double new_z = lmn_x.GetZ() * x_p + lmn_y.GetZ() * y_p + lmn_z.GetZ() * z_p + c.GetZ();

    Coordinate new_coordinate(new_x, new_y, new_z);
    return new_coordinate;
}

// This is only here until i SORT OUT THE STUOPID OVERLOADS IN COORDINATE CLASS
Coordinate cds::subtractCoordinates(const Coordinate& minuaend, const Coordinate& subtrahend)
{
    Coordinate new_coordinate((minuaend.GetX() - subtrahend.GetX()), (minuaend.GetY() - subtrahend.GetY()),
                              (minuaend.GetZ() - subtrahend.GetZ()));
    return new_coordinate;
}

double cds::CalculateMaxDistanceBetweenCoordinates(std::vector<Coordinate*> coords)
{
    double maxDistance = 0.0;
    for (std::vector<Coordinate*>::iterator it1 = coords.begin(); it1 != coords.end(); ++it1)
    {
        Coordinate* coord1 = (*it1);
        for (std::vector<Coordinate*>::iterator it2 = it1; it2 != coords.end(); ++it2)
        {
            Coordinate* coord2 = (*it2);
            if (coord1->Distance(coord2) > maxDistance)
            {
                maxDistance = coord1->Distance(coord2);
            }
        }
    }
    return maxDistance;
}

// I meant to time which is faster
bool cds::CheckIfOtherCoordinateIsWithinDistance(const Coordinate* a, const Coordinate* b, const double distance)
{
    double xDiff = a->GetX() - b->GetX();
    double yDiff = a->GetY() - b->GetY();
    double zDiff = a->GetZ() - b->GetZ();
    if ((xDiff * xDiff + yDiff * yDiff + zDiff * zDiff) < distance * distance)
    {
        return true;
    }
    return false;
}

// bool cds::CheckIfOtherCoordinateIsWithinDistanceA(const Coordinate* a, const Coordinate* b, const double distance)
//{
//     double xDiff = std::abs(a->GetX() - b->GetX());
//     double yDiff = std::abs(a->GetY() - b->GetY());
//     double zDiff = std::abs(a->GetZ() - b->GetZ());
//     if ((xDiff < distance) && (yDiff < distance) && (zDiff < distance))
//     {
//         return cds::CheckIfOtherCoordinateIsWithinDistance(a, b, distance);
//     }
//     return false;
// }
