#include "../includes/MolecularModeling/superimposition.hpp"

// OG: Have something else I'm changing, so can't update this now too. But there is lots of repeated code here.
// OG: Should just have an AtomVector based function and have everything else call that.

using namespace MolecularModeling;

/*
void gmml::GenerateMatrixFromAssembyCoordinates(Assembly *assembly, Eigen::Matrix3Xd *matrix)
{
    int col = 0; // Column index for matrix
    CoordinateVector assemblyCoordinates = assembly->GetAllCoordinates();
    for(CoordinateVector::iterator it = assemblyCoordinates.begin(); it != assemblyCoordinates.end(); it++)
    {
        (*matrix)(0, col) = (*it)->GetX();
        (*matrix)(1, col) = (*it)->GetY();
        (*matrix)(2, col) = (*it)->GetZ();
        col++;
    }
}

void gmml::ReplaceAssemblyCoordinatesFromMatrix(Assembly *assembly, Eigen::Matrix3Xd *matrix)
{
    int col = 0;
    CoordinateVector assemblyCoordinates = assembly->GetAllCoordinates();
    for(CoordinateVector::iterator it = assemblyCoordinates.begin(); it != assemblyCoordinates.end(); it++)
    {
        (*it)->SetX( (*matrix)(0, col) );
        (*it)->SetY( (*matrix)(1, col) );
        (*it)->SetZ( (*matrix)(2, col) );
        col++;
    }
}
*/

void gmml::GenerateMatrixFromAtomVectorCoordinates(AtomVector *atoms, Eigen::Matrix3Xd *matrix)
{
    int col = 0; // Column index for matrix
    for(AtomVector::iterator it = atoms->begin(); it != atoms->end(); it++)
    {
        (*matrix)(0, col) = (*it)->GetCoordinates().at(0)->GetX();
        (*matrix)(1, col) = (*it)->GetCoordinates().at(0)->GetY();
        (*matrix)(2, col) = (*it)->GetCoordinates().at(0)->GetZ();
        col++;
    }
}


void gmml::ReplaceAtomVectorCoordinatesFromMatrix(AtomVector *atoms, Eigen::Matrix3Xd *matrix)
{
    int col = 0; // Column index for matrix
    for(AtomVector::iterator it = atoms->begin(); it != atoms->end(); it++)
    {
        (*it)->GetCoordinates().at(0)->SetX( (*matrix)(0, col) );
        (*it)->GetCoordinates().at(0)->SetY( (*matrix)(1, col) );
        (*it)->GetCoordinates().at(0)->SetZ( (*matrix)(2, col) );
        col++;
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

void gmml::Superimpose(AtomVector moving, AtomVector target)
{
    Eigen::Matrix3Xd movingMatrix(3, moving.size()), targetMatrix(3, target.size()) ;

    //Create Matrices containing co-ordinates of moving and target
    GenerateMatrixFromAtomVectorCoordinates(&moving, &movingMatrix);
    GenerateMatrixFromAtomVectorCoordinates(&target, &targetMatrix);

    // Figure out how to move assembly moving onto target
    Eigen::Affine3d Affine = Find3DAffineTransform(movingMatrix, targetMatrix);

    // Create a matirx containing the moved co-ordinates of assembly moving.
    Eigen::Matrix3Xd movedMatrix = (Affine * movingMatrix);

    // Replace co-ordinates of moving with moved co-ordinates
    ReplaceAtomVectorCoordinatesFromMatrix(&moving, &movedMatrix);
}

void gmml::Superimpose(AtomVector moving, AtomVector target, AtomVector alsoMoving)
{
    Eigen::Matrix3Xd movingMatrix(3, moving.size()), targetMatrix(3, target.size());
    Eigen::Matrix3Xd alsoMovingMatrix(3, alsoMoving.size()); // separate from above line for clarity

    // Create a matrices containing co-ordinates of assembly moving and target
    GenerateMatrixFromAtomVectorCoordinates(&moving, &movingMatrix);
    GenerateMatrixFromAtomVectorCoordinates(&target, &targetMatrix);
    GenerateMatrixFromAtomVectorCoordinates(&alsoMoving, &alsoMovingMatrix);

    // Figure out how to move assembly moving onto target
    Eigen::Affine3d Affine = Find3DAffineTransform(movingMatrix, targetMatrix);

    // Create a matrix containing the moved co-ordinates of assembly moving and alsoMoving.
    Eigen::Matrix3Xd movedMatrix = (Affine * movingMatrix);
    Eigen::Matrix3Xd alsoMovedMatrix = (Affine * alsoMovingMatrix);

    // Replace co-ordinates of moving with moved co-ordinates
    ReplaceAtomVectorCoordinatesFromMatrix(&moving, &movedMatrix);
    ReplaceAtomVectorCoordinatesFromMatrix(&alsoMoving, &alsoMovedMatrix);
}

void gmml::Superimpose(MolecularModeling::Assembly *moving, MolecularModeling::Assembly *target)
{
    AtomVector moving_atoms = moving->GetAllAtomsOfAssembly();
    AtomVector target_atoms = target->GetAllAtomsOfAssembly();

    gmml::Superimpose(moving_atoms, target_atoms);

    /*
    Eigen::Matrix3Xd movingMatrix(3, moving->GetAllCoordinates().size()), targetMatrix(3, target->GetAllCoordinates().size());

    // Create a matrices containing co-ordinates of assembly moving and target
    GenerateMatrixFromAssembyCoordinates(moving, &movingMatrix);
    GenerateMatrixFromAssembyCoordinates(target, &targetMatrix);

    // Figure out how to move assembly moving onto target
    Eigen::Affine3d Affine = Find3DAffineTransform(movingMatrix, targetMatrix);

    // Create a matirx containing the moved co-ordinates of assembly moving.
    Eigen::Matrix3Xd movedMatrix = (Affine * movingMatrix);

    // Replace co-ordinates of moving with moved co-ordinates
    ReplaceAssemblyCoordinatesFromMatrix(moving, &movedMatrix);
    */
}

void gmml::Superimpose(Assembly *moving, Assembly *target, Assembly *alsoMoving)
{
    AtomVector moving_atoms = moving->GetAllAtomsOfAssembly();
    AtomVector target_atoms = target->GetAllAtomsOfAssembly();
    AtomVector alsoMoving_atoms = alsoMoving->GetAllAtomsOfAssembly();

    gmml::Superimpose(moving_atoms, target_atoms, alsoMoving_atoms);

    /*

    Eigen::Matrix3Xd movingMatrix(3, moving->GetAllCoordinates().size()), targetMatrix(3, target->GetAllCoordinates().size());
    Eigen::Matrix3Xd alsoMovingMatrix(3, alsoMoving->GetAllCoordinates().size()); // separate from above line for clarity

    // Create a matrices containing co-ordinates of assembly moving and target
    GenerateMatrixFromAssembyCoordinates(moving, &movingMatrix);
    GenerateMatrixFromAssembyCoordinates(target, &targetMatrix);
    GenerateMatrixFromAssembyCoordinates(alsoMoving, &alsoMovingMatrix);

    // Figure out how to move assembly moving onto target
    Eigen::Affine3d Affine = Find3DAffineTransform(movingMatrix, targetMatrix);

    // Create a matrix containing the moved co-ordinates of assembly moving and alsoMoving.
    Eigen::Matrix3Xd movedMatrix = (Affine * movingMatrix);
    Eigen::Matrix3Xd alsoMovedMatrix = (Affine * alsoMovingMatrix);


    // Replace co-ordinates of moving with moved co-ordinates
    ReplaceAssemblyCoordinatesFromMatrix(moving, &movedMatrix);
    ReplaceAssemblyCoordinatesFromMatrix(alsoMoving, &alsoMovedMatrix);
    */
}

void gmml::Superimpose(Assembly *moving, Assembly *target, AssemblyVector *alsoMoving)
{

    AtomVector moving_atoms = moving->GetAllAtomsOfAssembly();
    AtomVector target_atoms = target->GetAllAtomsOfAssembly();

    Eigen::Matrix3Xd movingMatrix(3, moving_atoms.size()), targetMatrix(3, target_atoms.size());

        // Create a matrices containing co-ordinates of assembly moving and target
    GenerateMatrixFromAtomVectorCoordinates(&moving_atoms, &movingMatrix);
    GenerateMatrixFromAtomVectorCoordinates(&target_atoms, &targetMatrix);

    // Figure out how to move assembly moving onto target
    Eigen::Affine3d Affine = Find3DAffineTransform(movingMatrix, targetMatrix);

    // Create a matrix containing the moved co-ordinates of assembly moving and alsoMoving.
    Eigen::Matrix3Xd movedMatrix = (Affine * movingMatrix);

    // Replace co-ordinates of moving with moved co-ordinates
    ReplaceAtomVectorCoordinatesFromMatrix(&moving_atoms, &movedMatrix);

    // Also move every assembly in also moving
    for(AssemblyVector::iterator it = alsoMoving->begin(); it != alsoMoving->end(); it++)
    {
        AtomVector alsoMoving_atoms = (*it)->GetAllAtomsOfAssembly();
        Eigen::Matrix3Xd alsoMovingMatrix(3, alsoMoving_atoms.size());
        GenerateMatrixFromAtomVectorCoordinates(&alsoMoving_atoms, &alsoMovingMatrix);
        Eigen::Matrix3Xd alsoMovedMatrix = (Affine * alsoMovingMatrix);
        ReplaceAtomVectorCoordinatesFromMatrix(&alsoMoving_atoms, &alsoMovedMatrix);
    }
}

/*
void gmml::Superimpose(Assembly *moving, Assembly *target, AssemblyVector *alsoMoving)
{

    Eigen::Matrix3Xd movingMatrix(3, moving->GetAllCoordinates().size()), targetMatrix(3, target->GetAllCoordinates().size());

    // Create a matrices containing co-ordinates of assembly moving and target
    GenerateMatrixFromAssembyCoordinates(moving, &movingMatrix);
    GenerateMatrixFromAssembyCoordinates(target, &targetMatrix);


    // Figure out how to move assembly moving onto target
    Eigen::Affine3d Affine = Find3DAffineTransform(movingMatrix, targetMatrix);

    // Create a matrix containing the moved co-ordinates of assembly moving and alsoMoving.
    Eigen::Matrix3Xd movedMatrix = (Affine * movingMatrix);

    // Replace co-ordinates of moving with moved co-ordinates
    ReplaceAssemblyCoordinatesFromMatrix(moving, &movedMatrix);

    // Also move every assembly in also moving
    for(AssemblyVector::iterator it = alsoMoving->begin(); it != alsoMoving->end(); it++)
    {
        Assembly *assembly = (*it); // I hate (*it)->
        Eigen::Matrix3Xd alsoMovingMatrix(3, assembly->GetAllCoordinates().size());
        GenerateMatrixFromAssembyCoordinates(assembly, &alsoMovingMatrix);
        Eigen::Matrix3Xd alsoMovedMatrix = (Affine * alsoMovingMatrix);
        ReplaceAssemblyCoordinatesFromMatrix(assembly, &alsoMovedMatrix);
    }
}
*/

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

