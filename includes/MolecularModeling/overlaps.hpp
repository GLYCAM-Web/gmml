#ifndef OVERLAPS_HPP
#define OVERLAPS_HPP

//#include "assembly.hpp"
#include "atom.hpp"
#include "common.hpp"
//*******************************************

//typedef std::vector<MolecularModeling::Atom*> AtomVector;
//typedef std::vector<MolecularModeling::Assembly*> AssemblyVector;

//*******************************************

namespace gmml
{

    double CalculateAtomicOverlaps(MolecularModeling::AtomVector atomsA, MolecularModeling::AtomVector atomsB);

    double CalculateAtomicOverlaps(MolecularModeling::Atom *atomA, MolecularModeling::Atom *atomB, double radiusA = -0.1, double radiusB = -0.1);

}


#endif // SUPERIMPOSITION_HPP
