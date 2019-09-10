#include "../../includes/MolecularModeling/superimposition.hpp"

using MolecularModeling::Atom;
using MolecularModeling::AtomVector;
using MolecularModeling::Assembly;
using MolecularModeling::AssemblyVector;
using GeometryTopology::CoordinateVector;

CoordinateVector gmml::GetCoordinatesInAtomVector(AtomVector *atoms)
{
    CoordinateVector coordinates;
    for(AtomVector::iterator it = atoms->begin(); it != atoms->end(); ++it)
    {
        Atom *atom = *it;
        coordinates.push_back(atom->GetCoordinate());
    }
    return coordinates;
}

CoordinateVector gmml::GetCoordinatesInAssemblyVector(AssemblyVector *assemblies)
{
    CoordinateVector coordinates;
    for(AssemblyVector::iterator it1 = assemblies->begin(); it1 != assemblies->end(); ++it1)
    {
        Assembly *assembly = *it1;
        CoordinateVector extra_coordinates = assembly->GetAllCoordinates();
        coordinates.insert( coordinates.end(), extra_coordinates.begin(), extra_coordinates.end() ); // "Concatentate"
    }
    return coordinates;
}

//void gmml::GenerateMatrixFromAtomVectorCoordinates(AtomVector *atoms, Eigen::Matrix3Xd *matrix)
//{
//    int col = 0; // Column index for matrix
//    for(AtomVector::iterator it = atoms->begin(); it != atoms->end(); it++)
//    {
//        (*matrix)(0, col) = (*it)->GetCoordinates().at(0)->GetX();
//        (*matrix)(1, col) = (*it)->GetCoordinates().at(0)->GetY();
//        (*matrix)(2, col) = (*it)->GetCoordinates().at(0)->GetZ();
//        col++;
//    }
//}

void gmml::GenerateMatrixFromCoordinates(CoordinateVector *coordinates, Eigen::Matrix3Xd *matrix)
{
    int col = 0; // Column index for matrix
    for(CoordinateVector::iterator coordinate = coordinates->begin(); coordinate != coordinates->end(); ++coordinate)
    {
        (*matrix)(0, col) = (*coordinate)->GetX();
        (*matrix)(1, col) = (*coordinate)->GetY();
        (*matrix)(2, col) = (*coordinate)->GetZ();
        ++col;
    }
}

//void gmml::ReplaceAtomVectorCoordinatesFromMatrix(AtomVector *atoms, Eigen::Matrix3Xd *matrix)
//{
//    int col = 0; // Column index for matrix
//    for(AtomVector::iterator it = atoms->begin(); it != atoms->end(); it++)
//    {
//        (*it)->GetCoordinates().at(0)->SetX( (*matrix)(0, col) );
//        (*it)->GetCoordinates().at(0)->SetY( (*matrix)(1, col) );
//        (*it)->GetCoordinates().at(0)->SetZ( (*matrix)(2, col) );
//        col++;
//    }
//}

void gmml::ReplaceCoordinatesFromMatrix(CoordinateVector *coordinates, Eigen::Matrix3Xd *matrix)
{
    int col = 0; // Column index for matrix
    for(CoordinateVector::iterator coordinate = coordinates->begin(); coordinate != coordinates->end(); ++coordinate)
    {
        (*coordinate)->SetX( (*matrix)(0, col) );
        (*coordinate)->SetY( (*matrix)(1, col) );
        (*coordinate)->SetZ( (*matrix)(2, col) );
        ++col;
    }
}

Eigen::Affine3d gmml::Find3DAffineTransform(Eigen::Matrix3Xd in, Eigen::Matrix3Xd out)
{
    // Default output
    Eigen::Affine3d A;
    A.linear() = Eigen::Matrix3d::Identity(3, 3);
    A.translation() = Eigen::Vector3d::Zero();

    if (in.cols() != out.cols())
        throw "Find3DAffineTransform(): input data mis-match";

    // First find the scale, by finding the ratio of sums of some distances,
    // then bring the datasets to the same scale.
    double dist_in = 0, dist_out = 0;
    for (int col = 0; col < in.cols()-1; col++)
    {
        dist_in  += (in.col(col+1) - in.col(col)).norm();
        dist_out += (out.col(col+1) - out.col(col)).norm();
    }
    if (dist_in <= 0 || dist_out <= 0)
        return A;
    double scale = dist_out/dist_in;
    out /= scale;

    // Find the centroids then shift to the origin
    Eigen::Vector3d in_ctr = Eigen::Vector3d::Zero();
    Eigen::Vector3d out_ctr = Eigen::Vector3d::Zero();
    for (int col = 0; col < in.cols(); col++)
    {
        in_ctr  += in.col(col);
        out_ctr += out.col(col);
    }
    in_ctr /= in.cols();
    out_ctr /= out.cols();
    for (int col = 0; col < in.cols(); col++)
    {
        in.col(col)  -= in_ctr;
        out.col(col) -= out_ctr;
    }

    // SVD
    Eigen::MatrixXd Cov = in * out.transpose();
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(Cov, Eigen::ComputeThinU | Eigen::ComputeThinV);

    // Find the rotation
    double d = (svd.matrixV() * svd.matrixU().transpose()).determinant();
    if (d > 0)
        d = 1.0;
    else
        d = -1.0;
    Eigen::Matrix3d I = Eigen::Matrix3d::Identity(3, 3);
    I(2, 2) = d;
    Eigen::Matrix3d R = svd.matrixV() * I * svd.matrixU().transpose();

    // The final transform
    A.linear() = scale * R;
    A.translation() = scale*(out_ctr - R*in_ctr);

    return A;
}

void gmml::Superimpose(CoordinateVector moving, CoordinateVector target)
{
    Eigen::Matrix3Xd movingMatrix(3, moving.size()), targetMatrix(3, target.size()) ;

    //Create Matrices containing co-ordinates of moving and target
    GenerateMatrixFromCoordinates(&moving, &movingMatrix);
    GenerateMatrixFromCoordinates(&target, &targetMatrix);

    // Figure out how to move assembly moving onto target
    Eigen::Affine3d Affine = Find3DAffineTransform(movingMatrix, targetMatrix);

    // Create a matirx containing the moved co-ordinates of assembly moving.
    Eigen::Matrix3Xd movedMatrix = (Affine * movingMatrix);

    // Replace co-ordinates of moving with moved co-ordinates
    ReplaceCoordinatesFromMatrix(&moving, &movedMatrix);
}

void gmml::Superimpose(CoordinateVector moving, CoordinateVector target, CoordinateVector alsoMoving)
{
    Eigen::Matrix3Xd movingMatrix(3, moving.size()), targetMatrix(3, target.size());
    Eigen::Matrix3Xd alsoMovingMatrix(3, alsoMoving.size()); // separate from above line for clarity

    // Create a matrices containing co-ordinates of assembly moving and target
    GenerateMatrixFromCoordinates(&moving, &movingMatrix);
    GenerateMatrixFromCoordinates(&target, &targetMatrix);
    GenerateMatrixFromCoordinates(&alsoMoving, &alsoMovingMatrix);

    // Figure out how to move assembly moving onto target
    Eigen::Affine3d Affine = Find3DAffineTransform(movingMatrix, targetMatrix);

    // Create a matrix containing the moved co-ordinates of assembly moving and alsoMoving.
    Eigen::Matrix3Xd movedMatrix = (Affine * movingMatrix);
    Eigen::Matrix3Xd alsoMovedMatrix = (Affine * alsoMovingMatrix);

    // Replace co-ordinates of moving with moved co-ordinates
    ReplaceCoordinatesFromMatrix(&moving, &movedMatrix);
    ReplaceCoordinatesFromMatrix(&alsoMoving, &alsoMovedMatrix);
}

//void gmml::Superimpose(AtomVector moving, AtomVector target)
//{
//    Eigen::Matrix3Xd movingMatrix(3, moving.size()), targetMatrix(3, target.size()) ;

//    //Create Matrices containing co-ordinates of moving and target
//    GenerateMatrixFromAtomVectorCoordinates(&moving, &movingMatrix);
//    GenerateMatrixFromAtomVectorCoordinates(&target, &targetMatrix);

//    // Figure out how to move assembly moving onto target
//    Eigen::Affine3d Affine = Find3DAffineTransform(movingMatrix, targetMatrix);

//    // Create a matirx containing the moved co-ordinates of assembly moving.
//    Eigen::Matrix3Xd movedMatrix = (Affine * movingMatrix);

//    // Replace co-ordinates of moving with moved co-ordinates
//    ReplaceAtomVectorCoordinatesFromMatrix(&moving, &movedMatrix);
//}

//void gmml::Superimpose(AtomVector moving, AtomVector target, AtomVector alsoMoving)
//{
//    Eigen::Matrix3Xd movingMatrix(3, moving.size()), targetMatrix(3, target.size());
//    Eigen::Matrix3Xd alsoMovingMatrix(3, alsoMoving.size()); // separate from above line for clarity

//    // Create a matrices containing co-ordinates of assembly moving and target
//    GenerateMatrixFromAtomVectorCoordinates(&moving, &movingMatrix);
//    GenerateMatrixFromAtomVectorCoordinates(&target, &targetMatrix);
//    GenerateMatrixFromAtomVectorCoordinates(&alsoMoving, &alsoMovingMatrix);

//    // Figure out how to move assembly moving onto target
//    Eigen::Affine3d Affine = Find3DAffineTransform(movingMatrix, targetMatrix);

//    // Create a matrix containing the moved co-ordinates of assembly moving and alsoMoving.
//    Eigen::Matrix3Xd movedMatrix = (Affine * movingMatrix);
//    Eigen::Matrix3Xd alsoMovedMatrix = (Affine * alsoMovingMatrix);

//    // Replace co-ordinates of moving with moved co-ordinates
//    ReplaceAtomVectorCoordinatesFromMatrix(&moving, &movedMatrix);
//    ReplaceAtomVectorCoordinatesFromMatrix(&alsoMoving, &alsoMovedMatrix);
//}

void gmml::Superimpose(AtomVector moving, AtomVector target)
{
    CoordinateVector moving_coordinates = gmml::GetCoordinatesInAtomVector(&moving);
    CoordinateVector target_coordinates = gmml::GetCoordinatesInAtomVector(&target);

    gmml::Superimpose(moving_coordinates, target_coordinates);
}

void gmml::Superimpose(AtomVector moving, AtomVector target, AtomVector alsoMoving)
{
    CoordinateVector moving_coordinates = gmml::GetCoordinatesInAtomVector(&moving);
    CoordinateVector target_coordinates = gmml::GetCoordinatesInAtomVector(&target);
    CoordinateVector alsoMoving_coordinates = gmml::GetCoordinatesInAtomVector(&alsoMoving);

    gmml::Superimpose(moving_coordinates, target_coordinates, alsoMoving_coordinates);
}

void gmml::Superimpose(Assembly *moving, Assembly *target)
{
    CoordinateVector moving_coordinates = moving->GetAllCoordinates();
    CoordinateVector target_coordinates = target->GetAllCoordinates();

    gmml::Superimpose(moving_coordinates, target_coordinates);
}

void gmml::Superimpose(Assembly *moving, Assembly *target, Assembly *alsoMoving)
{
    CoordinateVector moving_coordinates = moving->GetAllCoordinates();
    CoordinateVector target_coordinates = target->GetAllCoordinates();
    CoordinateVector alsoMoving_coordinates = alsoMoving->GetAllCoordinates();

    gmml::Superimpose(moving_coordinates, target_coordinates, alsoMoving_coordinates);

}

void gmml::Superimpose(Assembly *moving, Assembly *target, AssemblyVector *alsoMoving)
{

    CoordinateVector moving_coordinates = moving->GetAllCoordinates();
    CoordinateVector target_coordinates = target->GetAllCoordinates();
    CoordinateVector alsoMoving_coordinates = gmml::GetCoordinatesInAssemblyVector(alsoMoving);

    gmml::Superimpose(moving_coordinates, target_coordinates, alsoMoving_coordinates);

}

//void gmml::Superimpose(Assembly *moving, Assembly *target, AssemblyVector *alsoMoving)
//{

//    AtomVector moving_atoms = moving->GetAllAtomsOfAssembly();
//    AtomVector target_atoms = target->GetAllAtomsOfAssembly();

//    Eigen::Matrix3Xd movingMatrix(3, moving_atoms.size()), targetMatrix(3, target_atoms.size());

//        // Create a matrices containing co-ordinates of assembly moving and target
//    GenerateMatrixFromAtomVectorCoordinates(&moving_atoms, &movingMatrix);
//    GenerateMatrixFromAtomVectorCoordinates(&target_atoms, &targetMatrix);

//    // Figure out how to move assembly moving onto target
//    Eigen::Affine3d Affine = Find3DAffineTransform(movingMatrix, targetMatrix);

//    // Create a matrix containing the moved co-ordinates of assembly moving and alsoMoving.
//    Eigen::Matrix3Xd movedMatrix = (Affine * movingMatrix);

//    // Replace co-ordinates of moving with moved co-ordinates
//    ReplaceAtomVectorCoordinatesFromMatrix(&moving_atoms, &movedMatrix);

//    // Also move every assembly in also moving
//    for(AssemblyVector::iterator it = alsoMoving->begin(); it != alsoMoving->end(); it++)
//    {
//        AtomVector alsoMoving_atoms = (*it)->GetAllAtomsOfAssembly();
//        Eigen::Matrix3Xd alsoMovingMatrix(3, alsoMoving_atoms.size());
//        GenerateMatrixFromAtomVectorCoordinates(&alsoMoving_atoms, &alsoMovingMatrix);
//        Eigen::Matrix3Xd alsoMovedMatrix = (Affine * alsoMovingMatrix);
//        ReplaceAtomVectorCoordinatesFromMatrix(&alsoMoving_atoms, &alsoMovedMatrix);
//    }
//}

// A function to test Find3DAffineTransform()
/* void gmml::TestFind3DAffineTransform()
        {
        // Create datasets with known transform
        Eigen::Matrix3Xd in(3, 100), out(3, 100);
        Eigen::Quaternion<double> Q(1, 3, 5, 2);
        Q.normalize();
        Eigen::Matrix3d R = Q.toRotationMatrix();
        double scale = 2.0;
        for (int row = 0; row < in.rows(); row++) {
            for (int col = 0; col < in.cols(); col++) {
                in(row, col) = Eigen::log( (2*row + 10.0)/sqrt(1.0*col + 4.0) + sqrt(col*1.0)/(row + 1.0) );
            }
        }
        Eigen::Vector3d S;
        S << -5, 6, -27;
        for (int col = 0; col < in.cols(); col++)
            out.col(col) = scale*R*in.col(col) + S;

        Eigen::Affine3d A = Find3DAffineTransform(in, out);

        // See if we got the transform we expected
        if ( (scale*R-A.linear()).cwiseAbs().maxCoeff() > 1e-13 ||
             (S-A.translation()).cwiseAbs().maxCoeff() > 1e-13)
            throw "Could not determine the affine transform accurately enough";
    }*/
