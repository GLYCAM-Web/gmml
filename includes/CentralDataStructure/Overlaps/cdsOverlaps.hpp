#ifndef GMML_INCLUDES_CENTRALDATASTRUCTURE_OVERLAPS_OVERLAPS_HPP
#define GMML_INCLUDES_CENTRALDATASTRUCTURE_OVERLAPS_OVERLAPS_HPP

#include "includes/CodeUtils/logging.hpp"
#include "includes/CodeUtils/constants.hpp"
#include "includes/CentralDataStructure/coordinate.hpp"
#include "includes/CentralDataStructure/atom.hpp"
#include "includes/CentralDataStructure/residue.hpp"

#include <sstream>
#include <vector>
#include <iostream>

namespace cds
{
using cds::Coordinate;
using cds::Atom;
// ToDo this should be in Coordinate.
bool CheckIfOtherCoordinateIsWithinDistance(const Coordinate* a, const Coordinate* b, const double distance);
bool CheckIfOtherCoordinateIsWithinDistanceA(const Coordinate* a, const Coordinate* b, const double distance);
double CalculateAtomicOverlaps(Atom* atomA, Atom* atomB, double radiusA = 0.0, double radiusB = 0.0);
double CalculateAtomicOverlaps(std::vector<Atom*> atomsA, std::vector<Atom*> atomsB, bool print = false);
double CalculateAtomicOverlapsBetweenNonBondedAtoms(std::vector<Atom*>& atomsA, std::vector<Atom*>& atomsB);
unsigned int CountOverlappingResidues(const std::vector<cds::Residue*>& residuesA, const std::vector<cds::Residue*>& residuesB);
unsigned int CountOverlappingAtoms(const std::vector<Atom*>& atomsA, const std::vector<Atom*>& atomsB);
unsigned int CountOverlappingAtoms(const std::vector<cds::Residue*>& residuesA, const std::vector<cds::Residue*>& residuesB);

} // namespace
#endif
