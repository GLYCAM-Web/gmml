#ifndef SUPERIMPOSITION_HPP
#define SUPERIMPOSITION_HPP

#include "../Eigen_Algebra_Template_Library/Geometry"
#include "../GeometryTopology/coordinate.hpp"
#include "assembly.hpp"
#include "atom.hpp"

//*******************************************

typedef std::vector<MolecularModeling::Atom*> AtomVector;
typedef std::vector<MolecularModeling::Assembly*> AssemblyVector;

//*******************************************

namespace gmml
{
    GeometryTopology::Coordinate::CoordinateVector GetCoordinatesInAtomVector(AtomVector *atoms);

    GeometryTopology::Coordinate::CoordinateVector GetCoordinatesInAssemblyVector(AssemblyVector *assemblies);

    void GenerateMatrixFromAssembyCoordinates(MolecularModeling::Assembly *assembly, Eigen::Matrix3Xd *matrix);

    void ReplaceAssemblyCoordinatesFromMatrix(MolecularModeling::Assembly *assembly, Eigen::Matrix3Xd *matrix);

    //void GenerateMatrixFromAtomVectorCoordinates(AtomVector *atoms, Eigen::Matrix3Xd *matrix);

    //void ReplaceAtomVectorCoordinatesFromMatrix(AtomVector *atoms, Eigen::Matrix3Xd *matrix);

    void GenerateMatrixFromCoordinates(GeometryTopology::Coordinate::CoordinateVector *coordinates, Eigen::Matrix3Xd *matrix);

    void ReplaceCoordinatesFromMatrix(GeometryTopology::Coordinate::CoordinateVector *coordinates, Eigen::Matrix3Xd *matrix);

    Eigen::Affine3d Find3DAffineTransform(Eigen::Matrix3Xd in, Eigen::Matrix3Xd out);

    void Superimpose(GeometryTopology::Coordinate::CoordinateVector moving, GeometryTopology::Coordinate::CoordinateVector target);

    void Superimpose(GeometryTopology::Coordinate::CoordinateVector moving, GeometryTopology::Coordinate::CoordinateVector target, GeometryTopology::Coordinate::CoordinateVector alsoMoving);

    void Superimpose(AtomVector moving, AtomVector target);

    void Superimpose(AtomVector moving, AtomVector target, AtomVector alsoMoving);

    void Superimpose(MolecularModeling::Assembly *moving, MolecularModeling::Assembly *target);

    void Superimpose(MolecularModeling::Assembly *moving, MolecularModeling::Assembly *target, MolecularModeling::Assembly *alsoMoving);

    void Superimpose(MolecularModeling::Assembly *moving, MolecularModeling::Assembly *target, AssemblyVector *alsoMoving);
    // A function to test Find3DAffineTransform()
    // void TestFind3DAffineTransform();
}


#endif // SUPERIMPOSITION_HPP
