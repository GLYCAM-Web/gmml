#ifndef SELECTIONS_H
#define SELECTIONS_H

//#include "../../../includes/gmml.hpp"
#include "gmml.hpp"
#include "../atom.hpp"
#include "../residue.hpp"
#include "../common.hpp" // For the typedefs...
#include <regex>

using namespace MolecularModeling;
namespace selection
{
AtomVector AtomsWithinDistanceOf(Atom *query_atom, double distance, AtomVector atoms);
void FindAtomsConnectingResidues(Atom *current_atom, Residue *second_residue, AtomVector *connecting_atoms, bool *found_neighbor);
//void FindAtomsInPathToCycle(Atom *previous_atom, Atom *current_atom, AtomVector *atom_path, bool *found_cycle_point, Atom *&cycle_point);
//void FindAtomsInPathToBackboneNAtom(Atom *previous_atom, Atom *current_atom, AtomVector *atom_path, bool *found_N_atom);
//Pass pointer by reference in *&cycle_point, as I need to modify the actual pointer and not a copy of it.
bool FindCyclePoint(Atom *previous_atom, Atom *current_atom, AtomVector *atom_path, bool *found_cycle_point, Atom *&cycle_point);
bool FindPathBetweenTwoAtoms(Atom *current_atom, Atom *target_atom, AtomVector *atom_path, bool *found);
//bool CheckIfCycle(Atom *previous_atom, Atom *current_atom, AtomVector *atom_path);
void ClearAtomDescriptions(Residue *residue);
AtomVector FindCyclePoints(Atom *atom);
bool FindRotationPointsForNonCycles(Atom *previous_atom, Atom *current_atom, AtomVector *rotation_points);
Atom* FindCyclePointNeighbor(const AtomVector atom_path, Atom *cycle_point);
Atom* FindAtomNeighborThatMatchesQuery(Atom *atom, std::string query);
}

#endif // SELECTIONS_H
