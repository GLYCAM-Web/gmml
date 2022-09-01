#include <eigen3/Eigen/Geometry>

#include "../../includes/GeometryTopology/geometrytopology.hpp"
#include "../../includes/utils.hpp"

using GeometryTopology::Coordinate;

Coordinate GeometryTopology::subtract_coordinates(const Coordinate& minuaend, const Coordinate& subtrahend)
{
    Coordinate new_coordinate( (minuaend.GetX()-subtrahend.GetX()), (minuaend.GetY()-subtrahend.GetY()), (minuaend.GetZ()-subtrahend.GetZ()) );
    return new_coordinate;
}

Coordinate GeometryTopology::get_cartesian_point_from_internal_coords(const Coordinate* a, const Coordinate* b, const Coordinate* c, double angle_Degrees, double dihedral_Degrees, double distance_Angstrom)
{
	return GeometryTopology::get_cartesian_point_from_internal_coords(*a, *b, *c, angle_Degrees, dihedral_Degrees, distance_Angstrom);
}

Coordinate GeometryTopology::get_cartesian_point_from_internal_coords(const Coordinate& a, const Coordinate& b, const Coordinate& c, double angle_Degrees, double dihedral_Degrees, double distance_Angstrom)
{
//	std::cout << "Distance: " << distance_Angstrom << std::endl;
//	std::cout << "Angle: " << angle_Degrees << std::endl;
//	std::cout << "Dihedral: " << dihedral_Degrees << std::endl;
    double theta_Radians = gmml::ConvertDegree2Radian(angle_Degrees);
    double phi_Radians = gmml::ConvertDegree2Radian(dihedral_Degrees);
    Coordinate lmn_x, lmn_y, lmn_z;
    double x_p, y_p, z_p;

    Coordinate cb = GeometryTopology::subtract_coordinates(b, c); // original
    Coordinate ba = GeometryTopology::subtract_coordinates(a, b); // original

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

    double new_x = lmn_x.GetX()*x_p + lmn_y.GetX()*y_p + lmn_z.GetX()*z_p + c.GetX();
    double new_y = lmn_x.GetY()*x_p + lmn_y.GetY()*y_p + lmn_z.GetY()*z_p + c.GetY();
    double new_z = lmn_x.GetZ()*x_p + lmn_y.GetZ()*y_p + lmn_z.GetZ()*z_p + c.GetZ();

    Coordinate new_coordinate ( new_x, new_y, new_z );
    return new_coordinate;
}

Coordinate GeometryTopology::get_cartesian_point_from_internal_coords(MolecularModeling::Atom *a, MolecularModeling::Atom *b, MolecularModeling::Atom *c,
                                                                                        double angle_Degrees, double dihedral_Degrees,  double distance_Angstrom)
{
    return GeometryTopology::get_cartesian_point_from_internal_coords(a->GetCoordinate(), b->GetCoordinate(), c->GetCoordinate(), angle_Degrees, dihedral_Degrees, distance_Angstrom);
}

double GeometryTopology::calculateDistanceFromPointToLineBetweenTwoPoints(const Coordinate& queryPoint, const Coordinate& linePointA, const Coordinate& linePointB)
{
//                           (queryPoint)
//                           q ------------------ The "parallelogram".
//                          /|
//                         / | height (distance between queryPoint and line)
//                        /  |
//            linePointA a-------------------b linePointB
    // Get vector (aq) from linePointA to queryPoint.
    // Get vector (ab) from linePointA to linePointB. It is a direction vector for the line
    // Calculate magnitude of vector product of aq ab. This gives area of the parallelogram that they form.
    // Area of parallelogram is also base x height. We want to know height. We know base is distance between linePointA and linePointB.
    // h = area / base. h = | aq x ab | / | ab |.
    // || indicates magnitude

    // Get vector (aq) from linePointA to queryPoint.
    Eigen::Vector3d aq(queryPoint.GetX() - linePointA.GetX(), queryPoint.GetY() - linePointA.GetY(), queryPoint.GetZ() - linePointA.GetZ());
    // Get vector (ab) from linePointA to linePointB. It is a direction vector for the line.
    Eigen::Vector3d ab(linePointB.GetX() - linePointA.GetX(), linePointB.GetY() - linePointA.GetY(), linePointB.GetZ() - linePointA.GetZ());
    // Calculate magnitude of vector product of aq ab. This gives area of the parallelogram that they form.
    Eigen::Vector3d aq_x_ab = aq.cross(ab);
    double area =  aq_x_ab.norm(); // Area Of Parallelogram // Norm is magnitude aka length.

    double base = linePointA.Distance(linePointB); // ab.norm() should work too.
    double height = area / base;

    //std::cout << "Distance: " << height << " base is " << ab.norm() << " or " << base << "\n";
    return height;
}

Coordinate GeometryTopology::CreateMissingCoordinateForTetrahedralAtom(Coordinate *centralCoord, CoordinateVector threeNeighbors)
{
    Eigen::Vector3d combinedVs;
    for(auto &neighbor : threeNeighbors)
    {
        Eigen::Vector3d temp(centralCoord->GetX() - neighbor->GetX(),
                             centralCoord->GetY() - neighbor->GetY(),
                             centralCoord->GetZ() - neighbor->GetZ());
        combinedVs += temp;
    }
    combinedVs.normalize();
    Coordinate newCoord(centralCoord->GetX() + combinedVs(0),
                        centralCoord->GetY() + combinedVs(1),
                        centralCoord->GetZ() + combinedVs(2));
    return newCoord;
}

double GeometryTopology::CalculateDihedralAngle(GeometryTopology::Coordinate* a1, GeometryTopology::Coordinate* a2, GeometryTopology::Coordinate* a3, GeometryTopology::Coordinate* a4, bool returnRadians)
{   // returns Degrees by default, must set bool returnRadians to true for radians.
    GeometryTopology::Coordinate b1 = a2;
    b1.operator -(*a1);
    GeometryTopology::Coordinate b2 = a3;
    b2.operator -(*a2);
    GeometryTopology::Coordinate b3 = a4;
    b3.operator -(*a3);
    GeometryTopology::Coordinate b4 = b2;
    b4.operator *(-1);

    GeometryTopology::Coordinate b2xb3 = b2;
    b2xb3.CrossProduct(b3);

    GeometryTopology::Coordinate b1_m_b2n = b1;
    b1_m_b2n.operator *(b2.length());

    GeometryTopology::Coordinate b1xb2 = b1;
    b1xb2.CrossProduct(b2);

    double current_dihedral_angle = atan2(b1_m_b2n.DotProduct(b2xb3), b1xb2.DotProduct(b2xb3));

    if (returnRadians)
    {
        return current_dihedral_angle;
    }
    return (current_dihedral_angle * (180 / gmml::PI_RADIAN) ); // Convert to degrees
}
