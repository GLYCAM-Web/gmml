#ifndef GMML_INCLUDES_CENTRALDATASTRUCTURE_OVERLAPS_OVERLAPS_HPP
#define GMML_INCLUDES_CENTRALDATASTRUCTURE_OVERLAPS_OVERLAPS_HPP

#include "includes/CodeUtils/logging.hpp"
#include "includes/CodeUtils/constants.hpp"
#include "includes/GeometryTopology/coordinate.hpp"
#include "includes/CentralDataStructure/atom.hpp"

#include <sstream>
#include <vector>
#include <iostream>

namespace cds
{
using GeometryTopology::Coordinate;
using cds::Atom;
// ToDo this should be in Coordinate.
bool CheckIfOtherCoordinateIsWithinDistance(const Coordinate* a, const Coordinate* b, const double distance);
bool CheckIfOtherCoordinateIsWithinDistanceA(const Coordinate* a, const Coordinate* b, const double distance);
double CalculateAtomicOverlaps(Atom* atomA, Atom* atomB, double radiusA = 0.0, double radiusB = 0.0);
double CalculateAtomicOverlaps(std::vector<Atom*> atomsA, std::vector<Atom*> atomsB, bool print = false);
double CalculateAtomicOverlapsBetweenNonBondedAtoms(std::vector<Atom*>& atomsA, std::vector<Atom*>& atomsB);
unsigned int CountOverlappingAtoms(std::vector<Atom*>& atomsA, std::vector<Atom*>& atomsB);
} // namespace
#endif
