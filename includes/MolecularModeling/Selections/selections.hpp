#ifndef SELECTIONS_H
#define SELECTIONS_H

#include "../atom.hpp"
#include "../atomnode.hpp"
#include "../residue.hpp"

namespace selection
{
MolecularModeling::AtomVector AtomsWithinDistanceOf(MolecularModeling::Atom *query_atom, double distance, MolecularModeling::AtomVector atoms);
void FindAtomsConnectingResidues(MolecularModeling::Atom *current_atom, MolecularModeling::Residue *second_residue, MolecularModeling::AtomVector *connecting_atoms, bool *found_neighbor);
//void FindAtomsInPathToCycle(MolecularModeling::Atom *previous_atom, MolecularModeling::Atom *current_atom, MolecularModeling::AtomVector *atom_path, bool *found_cycle_point, MolecularModeling::Atom *&cycle_point);
//void FindAtomsInPathToBackboneNAtom(MolecularModeling::Atom *previous_atom, MolecularModeling::Atom *current_atom, MolecularModeling::AtomVector *atom_path, bool *found_N_atom);
//Pass pointer by reference in *&cycle_point, as I need to modify the actual pointer and not a copy of it.
bool FindCyclePoint(MolecularModeling::Atom *previous_atom, MolecularModeling::Atom *current_atom, MolecularModeling::AtomVector *atom_path, bool *found_cycle_point, MolecularModeling::Atom *&cycle_point);
bool FindPathBetweenTwoAtoms(MolecularModeling::Atom *current_atom, MolecularModeling::Atom *target_atom, MolecularModeling::AtomVector *atom_path, bool *found);
//bool CheckIfCycle(MolecularModeling::Atom *previous_atom, MolecularModeling::Atom *current_atom, MolecularModeling::AtomVector *atom_path);
void ClearAtomDescriptions(MolecularModeling::Residue *residue);
MolecularModeling::AtomVector FindCyclePoints(MolecularModeling::Atom *atom);
bool FindRotationPointsForNonCycles(MolecularModeling::Atom *previous_atom, MolecularModeling::Atom *current_atom, MolecularModeling::AtomVector *rotation_points);
MolecularModeling::Atom* FindCyclePointNeighbor(const MolecularModeling::AtomVector atom_path, MolecularModeling::Atom *cycle_point);
//commented out the below because it is commented out in the .cc file, and breaks GEMS -DM
// MolecularModeling::Atom* FindAtomNeighborThatMatchesQuery(MolecularModeling::Atom *atom, std::string query);
MolecularModeling::Residue* FindResidue(MolecularModeling::Assembly &assembly, const std::string query);
double GetMaxDistanceBetweenAtoms(MolecularModeling::AtomVector atoms);
MolecularModeling::AtomVector GetAtomsCommonToBothAtomVectors(MolecularModeling::AtomVector a, MolecularModeling::AtomVector b);
MolecularModeling::AtomVector GetAtomsin_a_Notin_b_AtomVectors(MolecularModeling::AtomVector a, MolecularModeling::AtomVector b);
MolecularModeling::AtomVector FindOtherAtomsWithinMolecule(MolecularModeling::Atom *queryAtom);
}

#endif // SELECTIONS_H
