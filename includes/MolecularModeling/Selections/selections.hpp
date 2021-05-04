#ifndef SELECTIONS_H
#define SELECTIONS_H

#include "../atom.hpp"
#include "../atomnode.hpp"
#include "../residue.hpp"

namespace selection
{
using MolecularModeling::Atom; 
using MolecularModeling::AtomVector; 
using MolecularModeling::Residue; 
using MolecularModeling::ResidueVector; 
using MolecularModeling::Assembly;
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
//commented out the below because it is commented out in the .cc file, and breaks GEMS -DM
// Atom* FindAtomNeighborThatMatchesQuery(Atom *atom, std::string query);
Residue* FindResidue(Assembly &assembly, const std::string query);
double GetMaxDistanceBetweenAtoms(AtomVector atoms);
AtomVector GetAtomsCommonToBothAtomVectors(AtomVector a, AtomVector b);
AtomVector GetAtomsin_a_Notin_b_AtomVectors(AtomVector a, AtomVector b);
AtomVector FindOtherAtomsWithinMolecule(Atom *queryAtom);
 // A function that compares atom numbers to see which is higher:
bool compareAtomNumbers(Atom *a1, Atom *a2);
ResidueVector SortResidueNeighborsByAcendingConnectionAtomNumber(AtomVector neighboringAtoms);
class Branch // This is only temporary. This will be replaced by new Template Graph class.
{
public:
    Branch(Atom *rootAtom, int depth = 0) : depth_ (depth), maxDepth_ (depth), rootAtom_ (rootAtom) {branchFound = false;};
    inline Atom* GetEnd() {return endAtom_;}
    inline Atom* GetRoot() {return rootAtom_;}
    inline int GetDepth() {return depth_;}
    inline bool IsBranchFound() {return branchFound;}
    //inline void AddToPath(Atom* atom) {path_.push_back(atom);}
    //inline void SetEnd(Atom* atom, int depth) {endAtom_ = atom; depth_ = depth}
    inline void SetRoot(Atom* atom) {rootAtom_ = atom;}
    inline void SetEnd(Atom* atom) {endAtom_ = atom; branchFound = true;}
    inline void SetDepth(int depth) {depth_ = depth;}
    inline void ChangeDepth(int delta) 
    {
        depth_ += delta;
        if (depth_ > maxDepth_) 
            maxDepth_ = depth_;
    }
    inline bool AtMaxDepth() {return maxDepth_ == depth_;}
    
    // inline void RemoveBond(Atom *otherAtom) {atomNodePtr_->RemoveEdge(otherAtom->GetNode());}
private:
    int depth_;
    int maxDepth_;
    Atom *rootAtom_;
    Atom *endAtom_;
    bool branchFound;
    //AtomVector path_;
}; 
void FindEndsOfBranchesFromLinkageAtom(Atom *currentAtom, Atom* previousAtom, Branch *branch);
std::string GetNonCarbonHeavyAtomNumbered(AtomVector atoms, std::string queryNumber);
}
#endif // SELECTIONS_H
