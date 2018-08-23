#ifndef SUPERIMPOSITION_HPP
#define SUPERIMPOSITION_HPP

#include "../Eigen_Algebra_Template_Library/Geometry"
#include "../GeometryTopology/coordinate.hpp"
#include "assembly.hpp"
#include "atom.hpp"

//*******************************************

typedef std::vector<MolecularModeling::Atom*> AtomVector;
typedef std::vector<GeometryTopology::Coordinate*> CoordinateVector;
typedef std::vector<MolecularModeling::Assembly*> AssemblyVector;

//*******************************************

namespace gmml
{
    CoordinateVector GetCoordinatesInAtomVector(AtomVector *atoms);

    void GenerateMatrixFromAssembyCoordinates(MolecularModeling::Assembly *assembly, Eigen::Matrix3Xd *matrix);

    void ReplaceAssemblyCoordinatesFromMatrix(MolecularModeling::Assembly *assembly, Eigen::Matrix3Xd *matrix);

    void GenerateMatrixFromAtomVectorCoordinates(AtomVector *atoms, Eigen::Matrix3Xd *matrix);

    void ReplaceAtomVectorCoordinatesFromMatrix(AtomVector *atoms, Eigen::Matrix3Xd *matrix);

    void GenerateMatrixFromCoordinates(CoordinateVector *coordinates, Eigen::Matrix3Xd *matrix);

    void ReplaceCoordinatesFromMatrix(CoordinateVector *coordinates, Eigen::Matrix3Xd *matrix);

    Eigen::Affine3d Find3DAffineTransform(Eigen::Matrix3Xd in, Eigen::Matrix3Xd out);

    void Superimpose(CoordinateVector moving, CoordinateVector target);

    void Superimpose(CoordinateVector moving, CoordinateVector target, CoordinateVector alsoMoving);

    void Superimpose(AtomVector moving, AtomVector target);

    void Superimpose(AtomVector moving, AtomVector target, AtomVector alsoMoving);

    void Superimpose(MolecularModeling::Assembly *moving, MolecularModeling::Assembly *target);

    void Superimpose(MolecularModeling::Assembly *moving, MolecularModeling::Assembly *target, MolecularModeling::Assembly *alsoMoving);

    void Superimpose(MolecularModeling::Assembly *moving, MolecularModeling::Assembly *target, AssemblyVector *alsoMoving);
    // A function to test Find3DAffineTransform()
    // void TestFind3DAffineTransform();
}


#endif // SUPERIMPOSITION_HPP
