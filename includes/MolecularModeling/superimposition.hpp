#ifndef SUPERIMPOSITION_HPP
#define SUPERIMPOSITION_HPP

#include "../External_Libraries/Eigen_Algebra_Template_Library/Geometry"
#include "../GeometryTopology/coordinate.hpp"
#include "assembly.hpp"
#include "atom.hpp"

namespace gmml
{

using MolecularModeling::AtomVector;
using MolecularModeling::AssemblyVector;
using MolecularModeling::Assembly;

GeometryTopology::CoordinateVector GetCoordinatesInAtomVector(AtomVector *atoms);

GeometryTopology::CoordinateVector GetCoordinatesInAssemblyVector(AssemblyVector *assemblies);

void GenerateMatrixFromAssembyCoordinates(Assembly *assembly, Eigen::Matrix3Xd *matrix);

void ReplaceAssemblyCoordinatesFromMatrix(Assembly *assembly, Eigen::Matrix3Xd *matrix);

//void GenerateMatrixFromAtomVectorCoordinates(AtomVector *atoms, Eigen::Matrix3Xd *matrix);

//void ReplaceAtomVectorCoordinatesFromMatrix(AtomVector *atoms, Eigen::Matrix3Xd *matrix);

void GenerateMatrixFromCoordinates(GeometryTopology::CoordinateVector *coordinates, Eigen::Matrix3Xd *matrix);

void ReplaceCoordinatesFromMatrix(GeometryTopology::CoordinateVector *coordinates, Eigen::Matrix3Xd *matrix);

Eigen::Affine3d Find3DAffineTransform(Eigen::Matrix3Xd in, Eigen::Matrix3Xd out);

void Superimpose(GeometryTopology::CoordinateVector moving, GeometryTopology::CoordinateVector target);

void Superimpose(GeometryTopology::CoordinateVector moving, GeometryTopology::CoordinateVector target, GeometryTopology::CoordinateVector alsoMoving);

void Superimpose(AtomVector moving, AtomVector target);

void Superimpose(AtomVector moving, AtomVector target, AtomVector alsoMoving);

void Superimpose(Assembly *moving, Assembly *target);

void Superimpose(Assembly *moving, Assembly *target, Assembly *alsoMoving);

void Superimpose(Assembly *moving, Assembly *target, AssemblyVector *alsoMoving);
// A function to test Find3DAffineTransform()
// void TestFind3DAffineTransform();
}


#endif // SUPERIMPOSITION_HPP
