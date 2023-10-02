#include <eigen3/Eigen/Geometry>

#include "includes/GeometryTopology/geometrytopology.hpp"
#include "includes/utils.hpp"
#include "includes/MolecularMetadata/GLYCAM/bondlengthbytypepair.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/MolecularModeling/atomnode.hpp"

using GeometryTopology::Coordinate;

Coordinate GeometryTopology::subtract_coordinates(const Coordinate& minuaend, const Coordinate& subtrahend)
{
    Coordinate new_coordinate((minuaend.GetX() - subtrahend.GetX()), (minuaend.GetY() - subtrahend.GetY()),
                              (minuaend.GetZ() - subtrahend.GetZ()));
    return new_coordinate;
}

Coordinate GeometryTopology::get_cartesian_point_from_internal_coords(const Coordinate* a, const Coordinate* b,
                                                                      const Coordinate* c, double angle_Degrees,
                                                                      double dihedral_Degrees, double distance_Angstrom)
{
    return GeometryTopology::get_cartesian_point_from_internal_coords(*a, *b, *c, angle_Degrees, dihedral_Degrees,
                                                                      distance_Angstrom);
}

Coordinate GeometryTopology::get_cartesian_point_from_internal_coords(const Coordinate& a, const Coordinate& b,
                                                                      const Coordinate& c, double angle_Degrees,
                                                                      double dihedral_Degrees, double distance_Angstrom)
{
    //	std::cout << "Distance: " << distance_Angstrom << std::endl;
    //	std::cout << "Angle: " << angle_Degrees << std::endl;
    //	std::cout << "Dihedral: " << dihedral_Degrees << std::endl;
    double theta_Radians = gmml::ConvertDegree2Radian(angle_Degrees);
    double phi_Radians   = gmml::ConvertDegree2Radian(dihedral_Degrees);
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

    double new_x = lmn_x.GetX() * x_p + lmn_y.GetX() * y_p + lmn_z.GetX() * z_p + c.GetX();
    double new_y = lmn_x.GetY() * x_p + lmn_y.GetY() * y_p + lmn_z.GetY() * z_p + c.GetY();
    double new_z = lmn_x.GetZ() * x_p + lmn_y.GetZ() * y_p + lmn_z.GetZ() * z_p + c.GetZ();

    Coordinate new_coordinate(new_x, new_y, new_z);
    return new_coordinate;
}

Coordinate GeometryTopology::get_cartesian_point_from_internal_coords(MolecularModeling::Atom* a,
                                                                      MolecularModeling::Atom* b,
                                                                      MolecularModeling::Atom* c, double angle_Degrees,
                                                                      double dihedral_Degrees, double distance_Angstrom)
{
    return GeometryTopology::get_cartesian_point_from_internal_coords(
        a->GetCoordinate(), b->GetCoordinate(), c->GetCoordinate(), angle_Degrees, dihedral_Degrees, distance_Angstrom);
}

double GeometryTopology::calculateDistanceFromPointToLineBetweenTwoPoints(const Coordinate& queryPoint,
                                                                          const Coordinate& linePointA,
                                                                          const Coordinate& linePointB)
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
    // Area of parallelogram is also base x height. We want to know height. We know base is distance between linePointA
    // and linePointB. h = area / base. h = | aq x ab | / | ab |.
    // || indicates magnitude

    // Get vector (aq) from linePointA to queryPoint.
    Eigen::Vector3d aq(queryPoint.GetX() - linePointA.GetX(), queryPoint.GetY() - linePointA.GetY(),
                       queryPoint.GetZ() - linePointA.GetZ());
    // Get vector (ab) from linePointA to linePointB. It is a direction vector for the line.
    Eigen::Vector3d ab(linePointB.GetX() - linePointA.GetX(), linePointB.GetY() - linePointA.GetY(),
                       linePointB.GetZ() - linePointA.GetZ());
    // Calculate magnitude of vector product of aq ab. This gives area of the parallelogram that they form.
    Eigen::Vector3d aq_x_ab = aq.cross(ab);
    double area             = aq_x_ab.norm(); // Area Of Parallelogram // Norm is magnitude aka length.

    double base   = linePointA.Distance(linePointB); // ab.norm() should work too.
    double height = area / base;

    // std::cout << "Distance: " << height << " base is " << ab.norm() << " or " << base << "\n";
    return height;
}

Coordinate GeometryTopology::CreateMissingCoordinateForTetrahedralAtom(MolecularModeling::Atom* centralAtom,
                                                                       const double distance)
{
    if (centralAtom->GetNode()->GetNodeNeighbors().size() != 4)
    {
        std::stringstream ss;
        ss << "Error in CreateMissingCoordinateForTetrahedralAtom. centralAtom neighbors is "
           << centralAtom->GetNode()->GetNodeNeighbors().size() << " for " << centralAtom->GetId();
        gmml::log(__LINE__, __FILE__, gmml::ERR, ss.str());
        throw std::runtime_error(ss.str());
    }
    std::vector<Coordinate*> threeNeighborCoords;
    for (auto& neighbor : centralAtom->GetNode()->GetNodeNeighbors())
    {
        threeNeighborCoords.push_back(neighbor->GetCoordinate());
    }
    return GeometryTopology::CreateMissingCoordinateForTetrahedralAtom(centralAtom->GetCoordinate(),
                                                                       threeNeighborCoords, distance);
}

Coordinate GeometryTopology::CreateMissingCoordinateForTetrahedralAtom(Coordinate* centralCoord,
                                                                       CoordinateVector threeNeighbors,
                                                                       const double distance)
{
    Eigen::Vector3d combinedVs;
    for (auto& neighbor : threeNeighbors)
    {
        Eigen::Vector3d temp(centralCoord->GetX() - neighbor->GetX(), centralCoord->GetY() - neighbor->GetY(),
                             centralCoord->GetZ() - neighbor->GetZ());
        combinedVs += temp;
    }
    combinedVs.normalize();
    Coordinate newCoord(centralCoord->GetX() + (combinedVs(0) * distance),
                        centralCoord->GetY() + (combinedVs(1) * distance),
                        centralCoord->GetZ() + (combinedVs(2)) * distance);
    return newCoord;
}

double GeometryTopology::CalculateDihedralAngle(Coordinate* a1, Coordinate* a2, Coordinate* a3, Coordinate* a4,
                                                bool returnRadians)
{ // returns Degrees by default, must set bool returnRadians to true for radians.
    Coordinate b1 = a2;
    b1.operator-(*a1);
    Coordinate b2 = a3;
    b2.operator-(*a2);
    Coordinate b3 = a4;
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
    return (current_dihedral_angle * (180 / gmml::PI_RADIAN)); // Convert to degrees
}

void GeometryTopology::SetDihedralAngle(Coordinate* a1, Coordinate* a2, Coordinate* a3, Coordinate* a4,
                                        const double dihedral_angle, std::vector<Coordinate*>& movingCoords)
{
    Coordinate b1 = a2;
    b1.operator-(*a1);
    Coordinate b2 = a3;
    b2.operator-(*a2);
    Coordinate b3 = a4;
    b3.operator-(*a3);
    Coordinate b4 = b2;
    b4.operator*(-1);
    Coordinate b2xb3 = b2;
    b2xb3.CrossProduct(b3);
    Coordinate b1_m_b2n = b1;
    b1_m_b2n.operator*(b2.length());
    Coordinate b1xb2 = b1;
    b1xb2.CrossProduct(b2);
    double current_dihedral = atan2(b1_m_b2n.DotProduct(b2xb3), b1xb2.DotProduct(b2xb3));
    double** dihedral_angle_matrix =
        gmml::GenerateRotationMatrix(&b4, a2, current_dihedral - gmml::ConvertDegree2Radian(dihedral_angle));
    for (auto& coord : movingCoords)
    {
        double x = dihedral_angle_matrix[0][0] * coord->GetX() + dihedral_angle_matrix[0][1] * coord->GetY() +
                   dihedral_angle_matrix[0][2] * coord->GetZ() + dihedral_angle_matrix[0][3];
        double y = dihedral_angle_matrix[1][0] * coord->GetX() + dihedral_angle_matrix[1][1] * coord->GetY() +
                   dihedral_angle_matrix[1][2] * coord->GetZ() + dihedral_angle_matrix[1][3];
        double z = dihedral_angle_matrix[2][0] * coord->GetX() + dihedral_angle_matrix[2][1] * coord->GetY() +
                   dihedral_angle_matrix[2][2] * coord->GetZ() + dihedral_angle_matrix[2][3];
        coord->SetX(x);
        coord->SetY(y);
        coord->SetZ(z);
    }
}

void GeometryTopology::SetAngle(MolecularModeling::Atom* a, MolecularModeling::Atom* b, MolecularModeling::Atom* c,
                                const double angle)
{
    std::vector<MolecularModeling::Atom*> atomsToRotate;
    atomsToRotate.push_back(b);
    c->FindConnectedAtoms(atomsToRotate); // this is too slow, just get the coords here.
    std::vector<Coordinate*> coords;
    atomsToRotate.erase(atomsToRotate.begin()); // feck
    for (auto& atom : atomsToRotate)
    {
        coords.push_back(atom->GetCoordinate());
    }
    GeometryTopology::SetAngle(a->GetCoordinate(), b->GetCoordinate(), c->GetCoordinate(), angle, coords);
    return;
}

void GeometryTopology::SetAngle(Coordinate* a1, Coordinate* a2, Coordinate* a3, const double angle,
                                std::vector<Coordinate*> coordinatesToMove)
{
    double current_angle = 0.0;
    Coordinate b1        = *a1;
    b1.operator-(*a2);
    Coordinate b2 = *a3;
    b2.operator-(*a2);
    current_angle         = acos((b1.DotProduct(b2)) / (b1.length() * b2.length() + gmml::DIST_EPSILON));
    double rotation_angle = gmml::ConvertDegree2Radian(angle) - current_angle;
    Coordinate direction  = b1;
    direction.CrossProduct(b2);
    direction.Normalize();
    double** rotation_matrix = gmml::GenerateRotationMatrix(&direction, a2, rotation_angle);
    for (auto& atom_coordinate : coordinatesToMove)
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

void GeometryTopology::SetDistance(MolecularModeling::Atom* a, MolecularModeling::Atom* b)
{ // Figure out distance
    gmml::MolecularMetadata::GLYCAM::BondLengthByTypePairContainer bondLengthByTypePairContainer;
    double distance = bondLengthByTypePairContainer.GetBondLengthForAtomTypes(a->MolecularDynamicAtom::GetAtomType(),
                                                                              b->MolecularDynamicAtom::GetAtomType());
    // Figure out position of where the b atom should end up relative to a
    Coordinate c    = GeometryTopology::CreateMissingCoordinateForTetrahedralAtom(a, distance);
    // Figure out which atoms will move
    std::vector<MolecularModeling::Atom*> atomsToRotate;
    atomsToRotate.push_back(a);
    b->FindConnectedAtoms(atomsToRotate);       // this is too slow, just get the coords here.
    atomsToRotate.erase(atomsToRotate.begin()); // feck
    for (auto& atom : atomsToRotate)
    {
        atom->GetCoordinate()->Translate(c.GetX(), c.GetY(), c.GetZ());
    }
    return;
}
