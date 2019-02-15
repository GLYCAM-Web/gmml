#include "../../includes/GeometryTopology/geometrytopology.hpp"

GeometryTopology::Coordinate GeometryTopology::get_cartesian_point_from_internal_coords(MolecularModeling::Atom *a, MolecularModeling::Atom *b, MolecularModeling::Atom *c,
                                                                                        double theta_Degrees, double phi_Degrees, double distance_Angstrom)
{
    GeometryTopology::Coordinate x = (a->GetCoordinates().at(0));
    GeometryTopology::Coordinate y = (b->GetCoordinates().at(0));
    GeometryTopology::Coordinate z = (c->GetCoordinates().at(0));
    GeometryTopology::Coordinate new_coordinate;

    new_coordinate = new_coordinate.get_cartesian_point_from_internal_coords(
                x, y, z, theta_Degrees, phi_Degrees, distance_Angstrom);
//    new_coordinate.Print();
//    std::cout << "\n";

    return new_coordinate;



}

GeometryTopology::Coordinate GeometryTopology::get_cartesian_point_from_internal_coords(Coordinate a, Coordinate b, Coordinate c,
            double theta_Degrees, double phi_Degrees, double distance_Angstrom)
{     // theta is the angle between 3 atoms. Phi is the torsion between 4 atoms.

       //Convert from Degrees to Radians
       if ( theta_Degrees < 0.0 ) theta_Degrees += 360.0;
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

       Coordinate cb = subtract_coordinates(b, c); // original
       Coordinate ba = subtract_coordinates(a, b); // original

       lmn_y = ba;
       lmn_y.CrossProduct(cb);
       lmn_y.Normalize();

       lmn_z = cb;
       lmn_z.Normalize();

       lmn_x = lmn_z;
       lmn_x.CrossProduct(lmn_y);

       x_p = distance_Angstrom *  sin(theta_Radians) * cos(phi_Radians);
       y_p = distance_Angstrom * sin(theta_Radians) * sin(phi_Radians);
       z_p = distance_Angstrom * cos(theta_Radians);

       double new_x = lmn_x.GetX()*x_p + lmn_y.GetX()*y_p + lmn_z.GetX()*z_p + c.GetX();
       double new_y = lmn_x.GetY()*x_p + lmn_y.GetY()*y_p + lmn_z.GetY()*z_p + c.GetY();
       double new_z = lmn_x.GetZ()*x_p + lmn_y.GetZ()*y_p + lmn_z.GetZ()*z_p + c.GetZ();

       Coordinate new_coordinate ( new_x, new_y, new_z );
//! \todo Add these to the debugging mechanism once the DebugLevel class (or whatever) is implemented.
 std::cout << "   The NEW coords are:  " << std::endl;
 std::cout << "         X  :  " << new_coordinate.GetX() << std::endl;
 std::cout << "         Y  :  " << new_coordinate.GetY() << std::endl;
 std::cout << "         Z  :  " << new_coordinate.GetZ() << std::endl;

       return new_coordinate;
}
