#include "../../includes/GeometryTopology/geometrytopology.hpp"
#include "../../includes/External_Libraries/Eigen_Algebra_Template_Library/Geometry"
#include "../../includes/utils.hpp"
//void whywontthislink::test(){
//    return;
//}

using GeometryTopology::Coordinate;

Coordinate GeometryTopology::subtract_coordinates(Coordinate minuaend, Coordinate subtrahend)
{
    Coordinate new_coordinate( (minuaend.GetX()-subtrahend.GetX()), (minuaend.GetY()-subtrahend.GetY()), (minuaend.GetZ()-subtrahend.GetZ()) );
    return new_coordinate;
}

Coordinate GeometryTopology::get_cartesian_point_from_internal_coords(Coordinate a, Coordinate b, Coordinate c, double theta_Degrees, double phi_Degrees, double distance_Angstrom)
{     // theta is the angle between 3 atoms. Phi is the torsion between 4 atoms.
    //Convert from Degrees to Radians
    if ( theta_Degrees < 0.0 ) {theta_Degrees += 360.0;}
    double theta_Radians = gmml::ConvertDegree2Radian(theta_Degrees);
    double phi_Radians = gmml::ConvertDegree2Radian(phi_Degrees);
    //! \todo Add these to the debugging mechanism once the DebugLevel class (or whatever) is implemented.
    // std::cout << "   The three coords are:  " << std::endl;
    // std::cout << "      a:  " << std::endl;
    // std::cout << "         X  :  " << a.GetX() << std::endl;
    // std::cout << "         Y  :  " << a.GetY() << std::endl;
    // std::cout << "         Z  :  " << a.GetZ() << std::endl;
    // std::cout << "      b:  " << std::endl;
    // std::cout << "         X  :  " << b.GetX() << std::endl;
    // std::cout << "         Y  :  " << b.GetY() << std::endl;
    // std::cout << "         Z  :  " << b.GetZ() << std::endl;
    // std::cout << "      c:  " << std::endl;
    // std::cout << "         X  :  " << c.GetX() << std::endl;
    // std::cout << "         Y  :  " << c.GetY() << std::endl;
    // std::cout << "         Z  :  " << c.GetZ() << std::endl;

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
    //! \todo Add these to the debugging mechanism once the DebugLevel class (or whatever) is implemented.
    // std::cout << "   The NEW coords are:  " << std::endl;
    // std::cout << "         X  :  " << new_coordinate.GetX() << std::endl;
    // std::cout << "         Y  :  " << new_coordinate.GetY() << std::endl;
    // std::cout << "         Z  :  " << new_coordinate.GetZ() << std::endl;
    return new_coordinate;
}

Coordinate GeometryTopology::get_cartesian_point_from_internal_coords(MolecularModeling::Atom *a, MolecularModeling::Atom *b, MolecularModeling::Atom *c,
                                                                                        double theta_Degrees, double phi_Degrees,  double distance_Angstrom)
{
    return GeometryTopology::get_cartesian_point_from_internal_coords(a->GetCoordinate(), b->GetCoordinate(), c->GetCoordinate(), theta_Degrees, phi_Degrees, distance_Angstrom);
}

double GeometryTopology::calculateDistanceFromPointToLineBetweenTwoPoints(Coordinate queryPoint, Coordinate linePointA, Coordinate linePointB)
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

