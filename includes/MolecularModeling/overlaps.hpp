#ifndef OVERLAPS_HPP
#define OVERLAPS_HPP

#include "atom.hpp"
#include "common.hpp"

namespace gmml
{

    double CalculateAtomicOverlaps(MolecularModeling::AtomVector atomsA, MolecularModeling::AtomVector atomsB);

    double CalculateAtomicOverlaps(MolecularModeling::Atom *atomA, MolecularModeling::Atom *atomB, double radiusA = -0.1, double radiusB = -0.1);

}


#endif // OVERLAPS_HPP
