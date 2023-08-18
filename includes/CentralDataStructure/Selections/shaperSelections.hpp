#ifndef INCLUDES_CENTRALDATASTRUCTURE_SELECTIONS_SHAPERSELECTIONS_HPP_
#define INCLUDES_CENTRALDATASTRUCTURE_SELECTIONS_SHAPERSELECTIONS_HPP_

#include "includes/CentralDataStructure/atom.hpp"
#include "includes/CentralDataStructure/residue.hpp"
#include "includes/CentralDataStructure/Shapers/residueLinkage.hpp"
#include <vector>
#include <string>
#include <algorithm> //find
#include "templatedSelections.hpp"

namespace cdsSelections
{
    // Helper class
    class Branch // This should all be replaced by better code in the Template Graph Library. Looking at you Preston.
    {
      public:
        Branch(cds::Atom* rootAtom, int depth = 0) : depth_(depth), maxDepth_(depth), rootAtom_(rootAtom)
        {
            branchFound = false;
        };

        inline cds::Atom* GetEnd()
        {
            return endAtom_;
        }

        inline cds::Atom* GetRoot()
        {
            return rootAtom_;
        }

        inline int GetDepth()
        {
            return depth_;
        }

        inline bool IsBranchFound()
        {
            return branchFound;
        }

        // inline void AddToPath(Atom* atom) {path_.push_back(atom);}
        inline void SetRoot(cds::Atom* atom)
        {
            rootAtom_ = atom;
        }

        inline void SetEnd(cds::Atom* atom)
        {
            endAtom_    = atom;
            branchFound = true;
        }

        inline void SetDepth(int depth)
        {
            depth_ = depth;
        }

        inline void ChangeDepth(int delta)
        {
            depth_ += delta;
            if (depth_ > maxDepth_)
            {
                maxDepth_ = depth_;
            }
        }

        inline bool AtMaxDepth()
        {
            return maxDepth_ == depth_;
        }

      private:
        int depth_;
        int maxDepth_;
        cds::Atom* rootAtom_ = nullptr;
        cds::Atom* endAtom_  = nullptr;
        bool branchFound;
        // AtomVector path_;
    };

    // Free functions
    std::vector<cds::Atom*> FindCyclePoints(cds::Atom* atom, cds::Residue* residue);
    void FindEndsOfBranchesFromLinkageAtom(cds::Atom* currentAtom, cds::Atom* previousAtom, cds::Residue* residue,
                                           Branch* branch);
    bool FindCyclePoint(cds::Atom* previous_atom, cds::Residue* residue, cds::Atom* current_atom,
                        std::vector<cds::Atom*>* atom_path, bool* found_cycle_point, cds::Atom*& cycle_point);
    cds::Atom* FindCyclePointNeighbor(const std::vector<cds::Atom*> atom_path, cds::Atom* cycle_point,
                                      cds::Residue* cyclePointResidue);
    bool FindPathBetweenTwoAtoms(cds::Atom* current_atom, cds::Residue* currentResidue, cds::Atom* target_atom,
                                 cds::Residue* targetResidue, std::vector<cds::Atom*>* atom_path, bool* found);
    void FindAtomsConnectingResidues(cds::Atom* current_atom, const cds::Residue* currentResidue,
                                     const cds::Residue* otherResidue, std::vector<cds::Atom*>* connecting_atoms,
                                     bool* found_neighbor);
    void ClearAtomLabels(cds::Residue* residue);
    // ResidueLinkages
    cds::ResidueLinkage* selectLinkageWithIndex(std::vector<cds::ResidueLinkage>& inputLinkages,
                                                const long long unsigned int indexQuery);
    std::vector<cds::ResidueLinkage> SplitLinkagesIntoPermutants(std::vector<cds::ResidueLinkage>& inputLinkages);
} // namespace cdsSelections
#endif
