#ifndef INCLUDES_CENTRALDATASTRUCTURE_EDITORS_SUPERIMPOSITION_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_EDITORS_SUPERIMPOSITION_HPP

#include "includes/CentralDataStructure/coordinate.hpp"
#include <eigen3/Eigen/Geometry>

using cds::Coordinate;
namespace cds
{
void GenerateMatrixFromCoordinates(std::vector<Coordinate*> *coordinates, Eigen::Matrix3Xd *matrix);
void ReplaceCoordinatesFromMatrix(std::vector<Coordinate*> *coordinates, Eigen::Matrix3Xd *matrix);
Eigen::Affine3d Find3DAffineTransform(Eigen::Matrix3Xd in, Eigen::Matrix3Xd out);
void Superimpose(std::vector<Coordinate*> moving, std::vector<Coordinate*> target);
void Superimpose(std::vector<Coordinate*> moving, std::vector<Coordinate*> target, std::vector<Coordinate*> alsoMoving);
// A function to test Find3DAffineTransform()
// void TestFind3DAffineTransform();
}


#endif // SUPERIMPOSITION_HPP
