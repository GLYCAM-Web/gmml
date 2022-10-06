#ifndef GMML_INCLUDES_CENTRAL_DATA_STRUCTURE_SHAPERS_RESIDUE_LINKAGE_HPP
#define GMML_INCLUDES_CENTRAL_DATA_STRUCTURE_SHAPERS_RESIDUE_LINKAGE_HPP
/*
 * This class figures out the rotatable bonds between two residues
 * Starts/ends at the CA atoms in proteins. Looks for cycles (as they aren't rotatable).
 * Stores each rotatable bond as a RotatableDihedral object.
 */
//#include "includes/CentralDataStructure/cdsAtom.hpp"
//#include "includes/CentralDataStructure/cdsResidue.hpp"
#include "includes/MolecularModeling/Selections/selections.hpp"
#include "includes/CentralDataStructure/Shapers/rotatableDihedral.hpp"
#include "includes/CodeUtils/logging.hpp"

#include <vector>

namespace cds
{
template<class residueT, class atomT>
class ResidueLinkage
{
public:
    //////////////////////////////////////////////////////////
    //                       CONSTRUCTOR                    //
    //////////////////////////////////////////////////////////
    ResidueLinkage(residueT *nonReducingResidue1, residueT *reducingResidue2, bool reverseAtomsThatMove = true);
    ResidueLinkage(residueT *nonReducingResidue1, residueT *reducingResidue2, std::vector<atomT*> alsoMovingAtoms, bool reverseAtomsThatMove = true);
    //////////////////////////////////////////////////////////
    //                       ACCESSOR                       //
    //////////////////////////////////////////////////////////
    std::vector<RotatableDihedral<atomT>> GetRotatableDihedrals() const;
    std::vector<RotatableDihedral<atomT>> GetRotatableDihedralsWithMultipleRotamers() const;
    int GetNumberOfRotatableDihedrals() const;
    int GetNumberOfShapes(const bool likelyShapesOnly = false) const;
    inline residueT* GetFromThisResidue1() const {return from_this_residue1_;}
    inline residueT* GetToThisResidue2() const {return to_this_residue2_;}
    bool CheckIfConformer() const;
    inline bool GetIfExtraAtoms() const {return isExtraAtoms_;}
    inline std::vector<atomT*> GetExtraAtoms() {return extraAtomsThatMove_;}
    inline unsigned long long GetIndex() const {return index_;}
    inline std::string GetName() const;
    //////////////////////////////////////////////////////////
    //                       MUTATOR                        //
    //////////////////////////////////////////////////////////
    inline void SetRotatableDihedrals(std::vector<RotatableDihedral<atomT>> rotatableDihedrals) {RotatableDihedrals_ = rotatableDihedrals;}
    //////////////////////////////////////////////////////////
    //                       FUNCTIONS                      //
    //////////////////////////////////////////////////////////
    void SetDefaultShapeUsingMetadata();
    void SetRandomShapeUsingMetadata(bool useRanges = true);
    void SetSpecificShapeUsingMetadata(int shapeNumber, bool useRanges = false);
    void SetSpecificShape(std::string dihedralName, std::string selectedRotamer);
    void SetCustomDihedralAngles(std::vector <double> dihedral_angles);
    void SetShapeToPrevious();
    void SetRandomDihedralAngles();
    void DetermineAtomsThatMove();
    // Simple meaning you only check each RotatableDihedral in series, not every combination.
    void SimpleWiggle(std::vector<atomT*>& overlapAtomSet1, std::vector<atomT*>& overlapAtomSet2, const int angleIncrement = 5);
    void SimpleWiggleCurrentRotamers(std::vector<atomT*>& overlapAtomSet1, std::vector<atomT*>& overlapAtomSet2, const int angleIncrement = 5);
    inline void SetIndex(unsigned long long index) {index_ = index;}
    //////////////////////////////////////////////////////////
    //                       DISPLAY FUNCTION               //
    //////////////////////////////////////////////////////////
    std::string Print() const;
private:
    //////////////////////////////////////////////////////////
    //                    PRIVATE FUNCTIONS                 //
    //////////////////////////////////////////////////////////
    std::vector<residueT*> GetResidues() const;
    inline bool GetIfReversedAtomsThatMove() const {return reverseAtomsThatMove_;}
    atomT* GetFromThisConnectionAtom1() const;
    atomT* GetToThisConnectionAtom2() const;
    void SetIfReversedAtomsThatMove(bool reversedAtomsThatMove);
    void AddExtraAtomsThatMove(std::vector<atomT*> extraAtoms);
    void InitializeClass(residueT* from_this_residue1, residueT* to_this_residue2, bool reverseAtomsThatMove);
    bool CheckIfViableLinkage() const;
    std::vector<RotatableDihedral<atomT>> FindRotatableDihedralsConnectingResidues(atomT* from_this_connection_atom1, atomT* to_this_connection_atom2);
    //std::vector<atomT*> DealWithBranchesFromLinkages(std::vector<atomT*> linearLinkageAtoms, Atom *cycle_point1, Atom *cycle_point2);
    // Previous function generates a list of linearly connected atoms that define the rotatable bonds
    // This function splits that list into groups of 4 and creates RotatableDihedral objects
    std::vector<RotatableDihedral<atomT>> SplitAtomVectorIntoRotatableDihedrals(std::vector<atomT*> atoms);
    gmml::MolecularMetadata::GLYCAM::DihedralAngleDataVector FindMetadata(const atomT* from_this_connection_atom1, const atomT* to_this_connection_atom2) const;
    void AddMetadataToRotatableDihedrals(gmml::MolecularMetadata::GLYCAM::DihedralAngleDataVector metadata);
    void SetResidues(residueT* residue1, residueT* residue2);
    void SetConnectionAtoms(residueT* residue1, residueT* residue2);
    void SetConformerUsingMetadata(bool useRanges = false, int conformerNumber = 0);
    unsigned long long GenerateIndex();
    std::string DetermineLinkageNameFromResidueNames() const;
    inline void SetName(std::string name) {name_ = name;}
    //////////////////////////////////////////////////////////
    //                       ATTRIBUTES                     //
    //////////////////////////////////////////////////////////
    residueT* from_this_residue1_;
    residueT* to_this_residue2_;
    atomT* from_this_connection_atom1_;
    atomT* to_this_connection_atom2_;
    std::vector<RotatableDihedral<atomT>> RotatableDihedrals_;
    bool reverseAtomsThatMove_;
    std::vector<atomT*> extraAtomsThatMove_;
    bool isExtraAtoms_ = true;
    unsigned long long index_;
    std::string name_; //e.g. "DGalpb1-6DGlcpNAc"

};

template<class residueT, class atomT>
inline ResidueLinkage<residueT, atomT>::ResidueLinkage(residueT *nonReducingResidue1, residueT *reducingResidue2, bool reverseAtomsThatMove)
{
    isExtraAtoms_ = false;
    this->InitializeClass(nonReducingResidue1, reducingResidue2, reverseAtomsThatMove);
}

template<class residueT, class atomT>
inline ResidueLinkage<residueT, atomT>::ResidueLinkage(residueT *nonReducingResidue1, residueT *reducingResidue2, std::vector<atomT*> alsoMovingAtoms, bool reverseAtomsThatMove)
{ // Order of calling functions is important!
    this->AddExtraAtomsThatMove(alsoMovingAtoms);
    this->InitializeClass(nonReducingResidue1, reducingResidue2, reverseAtomsThatMove);
}

template<class residueT, class atomT>
std::vector<residueT*> ResidueLinkage<residueT, atomT>::GetResidues() const
{
    std::vector<residueT*> residues {from_this_residue1_, to_this_residue2_};
    return residues;
}

template<class residueT, class atomT>
std::vector<RotatableDihedral<atomT>> ResidueLinkage<residueT, atomT>::GetRotatableDihedralsWithMultipleRotamers() const
{
    std::vector<RotatableDihedral<atomT>> returningDihedrals;
    for (auto &entry : this->GetRotatableDihedrals())
    {
        if (entry.GetNumberOfRotamers() > 1)
        {
            returningDihedrals.push_back(entry);
        }
    }
    return returningDihedrals;
}

template<class residueT, class atomT>
std::vector<RotatableDihedral<atomT>> ResidueLinkage<residueT, atomT>::GetRotatableDihedrals() const
{
    if (RotatableDihedrals_.empty())
    {
        std::stringstream ss;
        ss << "Error: RotatableDihedrals in this linkage is empty: " << from_this_residue1_->getId() << "-" << to_this_residue2_->getId() << std::endl;
        gmml::log(__LINE__,__FILE__,gmml::ERR, ss.str());
        throw std::runtime_error(ss.str());
    }
    return RotatableDihedrals_;
}

template<class residueT, class atomT>
int ResidueLinkage<residueT, atomT>::GetNumberOfShapes( const bool likelyShapesOnly) const// Can have conformers (sets of rotamers) or permutations of rotamers
{
    int numberOfShapes = 1;
    if ( (RotatableDihedrals_.empty()) || (RotatableDihedrals_.at(0).GetMetadata().empty()) )
    {
        return numberOfShapes;
    }
    if (RotatableDihedrals_.at(0).GetMetadata().at(0).rotamer_type_.compare("permutation")==0)
    {
        for (auto &entry : RotatableDihedrals_)
        {
            numberOfShapes = (numberOfShapes * entry.GetNumberOfRotamers(likelyShapesOnly));
        }
    }
    else if (RotatableDihedrals_.at(0).GetMetadata().at(0).rotamer_type_.compare("conformer")==0)
    { // Conformer should mean that each dihedral will have the same number of metadata entries.
        //numberOfShapes = RotatableDihedrals_.size(); // This was correct for ASN for the wrong reason. 4 conformers and 4 dihedrals...
        numberOfShapes = RotatableDihedrals_.at(0).GetNumberOfRotamers(likelyShapesOnly);
    }
    return numberOfShapes;
}

template<class residueT, class atomT>
bool ResidueLinkage<residueT, atomT>::CheckIfConformer() const
{
    if (RotatableDihedrals_.empty())
    {
        std::string errorMessage = "Error in ResidueLinkage::checkIfConformer as RotatableDihedrals_.empty()\n";
        gmml::log(__LINE__,__FILE__,gmml::ERR, errorMessage);
        throw std::runtime_error(errorMessage);
    }
    else if (RotatableDihedrals_.at(0).GetMetadata().empty())
    {
        std::string errorMessage = "Error in ResidueLinkage::checkIfConformer as RotatableDihedrals_.at(0).GetMetadata().empty()\n";
        gmml::log(__LINE__,__FILE__,gmml::ERR, errorMessage);
        throw std::runtime_error(errorMessage);
    }
    else
    {
        if (RotatableDihedrals_.at(0).GetMetadata().at(0).rotamer_type_.compare("permutation")==0)
            return false;
        else
            return true;
    }
    return false; //Default to shut up the compiler. Shut up compiler gawd.
}

template<class residueT, class atomT>
std::string ResidueLinkage<residueT, atomT>::GetName() const
{
    if (!name_.empty())
    {
        return name_;
    }
    return this->DetermineLinkageNameFromResidueNames();
}

template<class residueT, class atomT>
void ResidueLinkage<residueT, atomT>::SetDefaultShapeUsingMetadata()
{
    for (auto &entry : RotatableDihedrals_)
    {
        entry.SetSpecificAngleEntryUsingMetadata(); // Default is first entry
    }
}

template<class residueT, class atomT>
void ResidueLinkage<residueT, atomT>::SetRandomShapeUsingMetadata(bool useRanges)
{
    if (RotatableDihedrals_.at(0).GetMetadata().at(0).rotamer_type_.compare("permutation")==0)
    {
        for (auto &entry : RotatableDihedrals_)
        {
            entry.SetRandomAngleEntryUsingMetadata(useRanges);
        }
    }
    else if (RotatableDihedrals_.at(0).GetMetadata().at(0).rotamer_type_.compare("conformer")==0)
    {
        int numberOfConformers = RotatableDihedrals_.at(0).GetMetadata().size();
        std::uniform_int_distribution<> distr(0, (numberOfConformers - 1)); // define the range
        int randomlySelectedConformerNumber = distr(rng);
        for (auto &entry : RotatableDihedrals_)
        {
            entry.SetSpecificAngleEntryUsingMetadata(useRanges, randomlySelectedConformerNumber);
        }
    }
}

template<class residueT, class atomT>
void ResidueLinkage<residueT, atomT>::SetSpecificShapeUsingMetadata(int shapeNumber, bool useRanges)
{
    for (auto &entry : RotatableDihedrals_)
    {
        entry.SetSpecificAngleEntryUsingMetadata(useRanges, shapeNumber);
    }
}

template<class residueT, class atomT>
void ResidueLinkage<residueT, atomT>::SetSpecificShape(std::string dihedralName, std::string selectedRotamer)
{
    for (auto &RotatableDihedral : RotatableDihedrals_)
    {
        // This will call RotatableDihedrals that don't have dihedralName (phi,psi), and nothing will happen. Hmmm.
        if (RotatableDihedral.SetSpecificShape(dihedralName, selectedRotamer))
            return; // Return once you manage to set a shape.
    }
    std::string errorMessage = "Did not set " + dihedralName + " to " + selectedRotamer + " as requested in ResidueLinkage::SetSpecificShape()";
    gmml::log(__LINE__,__FILE__,gmml::ERR, errorMessage);
    throw std::runtime_error(errorMessage);
}

template<class residueT, class atomT>
void ResidueLinkage<residueT, atomT>::SetCustomDihedralAngles(std::vector <double> dihedral_angles)
{
    if(dihedral_angles.size() == RotatableDihedrals_.size())
    {
        std::vector <double>::iterator dihedral_angle_iterator = dihedral_angles.begin();
        for (auto &RotatableDihedral : RotatableDihedrals_)
        {
            RotatableDihedral.SetDihedralAngle(*dihedral_angle_iterator);
            ++dihedral_angle_iterator;
        }
    }
    else
    {
        std::stringstream ss;
        ss << "ERROR attempted to set " << dihedral_angles.size() << " dihedrals for a residue linkage with " << RotatableDihedrals_.size() << " rotatable dihedrals\n These numbers need to match for this function to work\nSomething has gone wrong.";
        gmml::log(__LINE__,__FILE__,gmml::ERR, ss.str());
        throw std::runtime_error(ss.str());
    }
    return;
}

template<class residueT, class atomT>
void ResidueLinkage<residueT, atomT>::SetShapeToPrevious()
{
    for(typename std::vector<RotatableDihedral<atomT>>::iterator RotatableDihedral = RotatableDihedrals_.begin(); RotatableDihedral != RotatableDihedrals_.end(); ++RotatableDihedral)
    {
        RotatableDihedral->SetDihedralAngleToPrevious();
    }
    return;
}

template<class residueT, class atomT>
void ResidueLinkage<residueT, atomT>::SetRandomDihedralAngles()
{
    for(typename std::vector<RotatableDihedral<atomT>>::iterator RotatableDihedral = RotatableDihedrals_.begin(); RotatableDihedral != RotatableDihedrals_.end(); ++RotatableDihedral)
    {
        RotatableDihedral->RandomizeDihedralAngle();
    }
    return;
}

template<class residueT, class atomT>
void ResidueLinkage<residueT, atomT>::DetermineAtomsThatMove()
{
    for(typename std::vector<RotatableDihedral<atomT>>::iterator RotatableDihedral = RotatableDihedrals_.begin(); RotatableDihedral != RotatableDihedrals_.end(); ++RotatableDihedral)
    {
        RotatableDihedral->DetermineAtomsThatMove();
    }
    return;
}


template<class residueT, class atomT>
void ResidueLinkage<residueT, atomT>::SimpleWiggle(std::vector<atomT*>& overlapAtomSet1, std::vector<atomT*>& overlapAtomSet2, const int angleIncrement)
{
    for(auto &RotatableDihedral : this->GetRotatableDihedrals())
    {
        RotatableDihedral.WiggleUsingAllRotamers(overlapAtomSet1, overlapAtomSet2, angleIncrement);
    }
}

template<class residueT, class atomT>
void ResidueLinkage<residueT, atomT>::SimpleWiggleCurrentRotamers(std::vector<atomT*>& overlapAtomSet1, std::vector<atomT*>& overlapAtomSet2, const int angleIncrement)
{
    for(auto &RotatableDihedral : this->GetRotatableDihedrals())
     {
//        std::cout << "SimpleWiggling current rotamer: " << RotatableDihedral.Print() << "\n";
        RotatableDihedral.WiggleWithinCurrentRotamer(overlapAtomSet1, overlapAtomSet2, angleIncrement);
     }
}

template<class residueT, class atomT>
std::string ResidueLinkage<residueT, atomT>::Print() const
{
    std::stringstream ss;
    ss << "ResidueLinkage Index: " << this->GetIndex() << ", Name: " << this->GetName() << ", NumberOfShapes: " << this->GetNumberOfShapes()
              << ", ids: " << this->GetFromThisResidue1()->getId() << "@" << this->GetFromThisConnectionAtom1()->GetName()
              << " -- " << this->GetToThisResidue2()->getId() << "@" << this->GetToThisConnectionAtom2()->GetName() << "\n";
    gmml::log(__LINE__, __FILE__, gmml::INF, ss.str());
    for(auto & rotatableDihedral : this->GetRotatableDihedrals() )
    {
         ss << rotatableDihedral.Print();
    }
    return ss.str();
}

// PRIVATE

template<class residueT, class atomT>
void ResidueLinkage<residueT, atomT>::InitializeClass(residueT* from_this_residue1, residueT* to_this_residue2, bool reverseAtomsThatMove)
{
    //set local debug flag
    int local_debug = -1;


    this->SetResidues(from_this_residue1, to_this_residue2);
    this->SetIfReversedAtomsThatMove(reverseAtomsThatMove);
    this->SetConnectionAtoms(from_this_residue1_, to_this_residue2_);
    if(local_debug > 0)
    {
        // std::cout << "Maybe Finding connection between " << from_this_residue1->getId() << " :: " << to_this_residue2->getId() << std::endl;
        gmml::log(__LINE__, __FILE__, gmml::INF, "Maybe Finding connection between " + from_this_residue1->getId() + " :: " + to_this_residue2->getId());
    }
    if(this->CheckIfViableLinkage())
    {
//        std::cout << "Finding connection between " << from_this_residue1->getId() << " :: " << to_this_residue2->getId() << std::endl;
//        std::cout << "Connection atoms are from: " << from_this_connection_atom1_->getId() << " to " << to_this_connection_atom2_->getId() << std::endl;
        if(local_debug > 0)
        {
            gmml::log(__LINE__, __FILE__, gmml::INF, "Finding connection between " + from_this_residue1->getId() + " :: " + to_this_residue2->getId());
            gmml::log(__LINE__, __FILE__, gmml::INF, "Connection atoms are from: " + from_this_connection_atom1_->getId() + " to " + to_this_connection_atom2_->getId());
        }
        RotatableDihedrals_ = this->FindRotatableDihedralsConnectingResidues(from_this_connection_atom1_, to_this_connection_atom2_);
//        std::cout << "Finding metadata for " << from_this_residue1->getId() << " :: " << to_this_residue2->getId() << std::endl;
        if(local_debug > 0)
        {
            gmml::log(__LINE__, __FILE__, gmml::INF, "Finding metadata for " + from_this_residue1->getId() + " :: " + to_this_residue2->getId());
        }
        gmml::MolecularMetadata::GLYCAM::DihedralAngleDataVector metadata = this->FindMetadata(from_this_connection_atom1_, to_this_connection_atom2_);
      //  std::cout << "Metadata found:\n";
      //  for (auto &dihedralAngleData : metadata)
      //  {
      //      std::cout << dihedralAngleData.print() << std::endl;
      //  }
        if(local_debug > 0)
        {
            gmml::log(__LINE__, __FILE__, gmml::INF, "Metadata found:");
            for (auto &dihedralAngleData : metadata)
            {
                gmml::log(__LINE__, __FILE__, gmml::INF, dihedralAngleData.print());
            }
        }
        this->AddMetadataToRotatableDihedrals(metadata);
    }
    this->SetIndex(this->GenerateIndex());
    return;
}

template<class residueT, class atomT>
bool ResidueLinkage<residueT, atomT>::CheckIfViableLinkage() const
{
    for(auto &residue : this->GetResidues())
    {
        if (residue->GetAtoms().size() <= 1)
        { // If either linkage has only 1 atom, return false. Should not set dihedral.
            return false;
        }
    }
    return true;
}

template<class residueT, class atomT>
std::vector<RotatableDihedral<atomT>> ResidueLinkage<residueT, atomT>::FindRotatableDihedralsConnectingResidues(atomT* from_this_connection_atom1, atomT* to_this_connection_atom2)
{
    // Going to ignore tags etc.
    // Given two residues that are connected. Find connecting atoms.
    // Search neighbors other than connected atom. Ie search out in both directions, but remain within same residue.
    // Warning, residue may have fused cycles!
    // Will fail for non-protein residues without cycles. As don't have a non-rotatable bond to anchor from. Can code that later (and deal with branches from these residues).
//    std::cout << "Finding rot bonds for " << from_this_connection_atom1->GetResidue()->getId() << " and " << to_this_connection_atom2->GetResidue()->getId() << "\n";

    std::vector<atomT*> from_this_residue1_cycle_points = selection::FindCyclePoints(from_this_connection_atom1);
//    std::cout << "Moving onto second residue.\n";
    std::vector<atomT*> to_this_residue2_cycle_points = selection::FindCyclePoints(to_this_connection_atom2);
    // Need to reverse one of these, so when concatenated, they are ordered ok. This might not be ok.
    //std::reverse(to_this_residue2_cycle_points.begin(), to_this_residue2_cycle_points.end());
    std::reverse(from_this_residue1_cycle_points.begin(), from_this_residue1_cycle_points.end());
    // Now concatenate:
    from_this_residue1_cycle_points.insert( from_this_residue1_cycle_points.end(), to_this_residue2_cycle_points.begin(), to_this_residue2_cycle_points.end() );
    // Now that have a list of rotation points. Split into pairs and find rotatable bonds between them
    bool found = false;
    std::vector<atomT*> connecting_atoms = {from_this_connection_atom1, to_this_connection_atom2};
//     std::cout << "cycle point atoms are:\n";
//     for(auto & atom : from_this_residue1_cycle_points)
//     {
//         std::cout << atom->getId() << "\n";
//     }
//     std::cout << "\n";
    std::vector<RotatableDihedral<atomT>> rotatableDihedralsInBranches;
    for(int i = 0; i < from_this_residue1_cycle_points.size(); i = i+2)
    {
        atomT* cycle_point1 = from_this_residue1_cycle_points.at(i);
        atomT* cycle_point2 = from_this_residue1_cycle_points.at(i+1);

        found = false;
        connecting_atoms.clear();
//        std::cout << "Finding Path between:" << cycle_point1->getId() << " and " << cycle_point2->getId() << "\n";
        selection::FindPathBetweenTwoAtoms(cycle_point1, cycle_point2, &connecting_atoms, &found);
        selection::ClearAtomDescriptions(cycle_point1->GetResidue());
        selection::ClearAtomDescriptions(cycle_point2->GetResidue());
        // Find neighboring atoms needed to define dihedral. Pass in connecting atoms so don't find any of those.
        atomT* neighbor1 =  selection::FindCyclePointNeighbor(connecting_atoms, cycle_point1);
        atomT* neighbor2 =  selection::FindCyclePointNeighbor(connecting_atoms, cycle_point2);
        // Insert these neighbors into list of connecting atoms, at beginning and end of vector.
        // connecting_atoms gets populated as it falls out, so list is reversed from what you'd expect
        std::reverse(connecting_atoms.begin(), connecting_atoms.end());
        connecting_atoms.insert(connecting_atoms.begin(), neighbor1);
        connecting_atoms.push_back(neighbor2);

//         std::cout << "Updated Path between:\n " << cycle_point1->getId() << " and " << cycle_point2->getId() << "\n";
//         for (const auto& atom : connecting_atoms)
//             std::cout << atom->getId() << "\n";
//         std::cout << "\n";
        selection::ClearAtomDescriptions(cycle_point1->GetResidue());
        selection::ClearAtomDescriptions(cycle_point2->GetResidue());
        // This mess was made to address the branching in 2-7 and 2-8 linkages.
        // These branches are long enough that they need default torsions set.
        if (connecting_atoms.size() > 4) // Otherwise there are no torsions
        { // Only viable linkages. Throw if not >4?
            for(typename std::vector<atomT*>::iterator it = connecting_atoms.begin()+1; it != connecting_atoms.end()-1; ++it)
            { //
                atomT* connectionAtom = (*it);
                if ( (connectionAtom != cycle_point1) && (connectionAtom != cycle_point2) )
                {
                    for(auto &neighbor : connectionAtom->GetNode()->GetNodeNeighbors()) // Need an interator
                    {
                        if (std::find(connecting_atoms.begin(), connecting_atoms.end(), neighbor) == connecting_atoms.end()) // if not in the vector
                        {
                            selection::Branch branch(connectionAtom);
                            selection::FindEndsOfBranchesFromLinkageAtom(neighbor, connectionAtom, &branch);
                            if (branch.IsBranchFound())
                            {
                                found = false;
                                std::vector<atomT*> foundPath;
                                // This fills in foundPath:
                                selection::FindPathBetweenTwoAtoms(branch.GetRoot(), branch.GetEnd(), &foundPath, &found);
                                atomT* neighbor = selection::FindCyclePointNeighbor(foundPath, branch.GetRoot());
                                //foundPath.insert(foundPath.begin(), neighbor);
                                foundPath.push_back(neighbor);
//                                 std::cout << "Found atoms:\n";
//                                 for (auto &atom: foundPath)
//                                     std::cout << atom->getId() << "\n";
                                std::vector<RotatableDihedral<atomT>> temp = this->SplitAtomVectorIntoRotatableDihedrals(foundPath);
                                rotatableDihedralsInBranches.insert( rotatableDihedralsInBranches.end(), temp.begin(), temp.end() );
                            }
                        }
                    }
                }
            }
        } // End dealing with branching linkages
//         std::cout << "These are the assigned branched RotatableDihedrals:\n";
//         for (auto &dihedral : rotatableDihedralsInBranches)
//             dihedral.Print();
    }
    std::vector<RotatableDihedral<atomT>> RotatableDihedrals = this->SplitAtomVectorIntoRotatableDihedrals(connecting_atoms);
    // Add any linkage branches (in 2-7 and 2-8) to the rest.
    RotatableDihedrals.insert( RotatableDihedrals.end(), rotatableDihedralsInBranches.begin(), rotatableDihedralsInBranches.end() );
//     std::cout << "These are the assigned RotatableDihedrals:\n";
//     for (auto &dihedral : RotatableDihedrals)
//         dihedral.Print();
    return RotatableDihedrals;
}

template<class residueT, class atomT>
std::vector<RotatableDihedral<atomT>> ResidueLinkage<residueT, atomT>::SplitAtomVectorIntoRotatableDihedrals(std::vector<atomT*> atoms)
{
    //Ok looking for sets of four atoms, but shifting along vector by one atom for each dihedral.
    // So four atoms will make one rotatable bond, five will make two bonds, six will make three etc.
    std::vector<RotatableDihedral<atomT>> RotatableDihedrals_generated;
    if(atoms.size() < 4)
    {
        std::stringstream ss;
        ss << "ERROR in ResidueLinkage::SplitAtomVectorIntoRotatableDihedrals, not enough atoms in atom vector: " << atoms.size() << "\n";
        ss << "This should be 4 or something is very wrong\n";
        gmml::log(__LINE__,__FILE__,gmml::ERR,ss.str());
        throw std::runtime_error(ss.str());
    }
    else
    {
        for(typename std::vector<atomT*>::iterator it1 = atoms.begin(); it1 != (atoms.end()-3); ++it1)
        {
            atomT* atom1 = *it1;
            atomT* atom2 = *(it1+1);
            atomT* atom3 = *(it1+2);
            atomT* atom4 = *(it1+3);
            RotatableDihedrals_generated.emplace_back(atom1, atom2, atom3, atom4, this->GetIfReversedAtomsThatMove());
        }
    }
    return RotatableDihedrals_generated;
}

template<class residueT, class atomT>
gmml::MolecularMetadata::GLYCAM::DihedralAngleDataVector ResidueLinkage<residueT, atomT>::FindMetadata(const atomT* from_this_connection_atom1, const atomT* to_this_connection_atom2) const
{
    gmml::MolecularMetadata::GLYCAM::DihedralAngleDataContainer DihedralAngleMetadata;
    gmml::MolecularMetadata::GLYCAM::DihedralAngleDataVector matching_entries = DihedralAngleMetadata.GetEntriesForLinkage(from_this_connection_atom1, to_this_connection_atom2);
    //std::cout << "Found these " << matching_entries.size() << " entries:\n";
    // for (const auto& entry : matching_entries)
    // {
    //     std::cout << entry.index_ << " : " << entry.atom1_ << ", " << entry.atom2_ << ", " << entry.atom3_ << ", " << entry.atom4_ << ", " << entry.default_angle_value_ << "\n";
    // }
    if (matching_entries.empty())
    {
        std::stringstream ss;
        ss << "No Metadata entries found for connection between " << from_this_connection_atom1->getId() << " and " << to_this_connection_atom2->getId() << "\n";
        ss << "Note that order should be reducing atom - anomeric atom\n";
        gmml::log(__LINE__,__FILE__,gmml::ERR,ss.str());
        throw std::runtime_error(ss.str());
    }
    return matching_entries;
}

template<class residueT, class atomT>
void ResidueLinkage<residueT, atomT>::AddMetadataToRotatableDihedrals(gmml::MolecularMetadata::GLYCAM::DihedralAngleDataVector metadata)
{
    // Adding another sanity check to this insanity
    if (RotatableDihedrals_.size() > metadata.size())
    {
        std::stringstream ss;
        ss << "Problem with the metadata found in gmml for this linkage. The number of identified rotatable dihedrals: " << RotatableDihedrals_.size() << " is greater than the number of metadata items: " << metadata.size() << " found for this linkage:\n" << this->Print();
        gmml::log(__LINE__,__FILE__,gmml::WAR, ss.str());
        throw std::runtime_error(ss.str());
    }
    for (auto & rotatableDihedral : RotatableDihedrals_)
    {
        rotatableDihedral.ClearMetadata(); // First clear any metadata already in place.
    }
    for (const auto& entry : metadata)
    {
        int vector_position = (entry.number_of_bonds_from_anomeric_carbon_ - 1); // vectors start at 0.
//        std::cout << "Adding to position: "<< vector_position << " in vector of size: " << RotatableDihedrals_.size() << std::endl;
        if (vector_position < RotatableDihedrals_.size())
        {
            RotatableDihedrals_.at(vector_position).AddMetadata(entry);
//            std::cout << "Set " << entry.dihedral_angle_name_ << " meta data for "
//                      << RotatableDihedrals_.at(vector_position).GetAtoms().at(0)->getId() << ", "
//                      << RotatableDihedrals_.at(vector_position).GetAtoms().at(1)->getId() << ", "
//                      << RotatableDihedrals_.at(vector_position).GetAtoms().at(2)->getId() << ", "
//                      << RotatableDihedrals_.at(vector_position).GetAtoms().at(3)->getId() << ". ";
//            std::cout << "Index: " << entry.index_ << ", angle: " << entry.default_angle_value_ << "\n";
            //RotatableDihedrals_.at(vector_position).Print();
        }
        else
        {
            std::string message = "Tried to add metadata to a rotatable bond that does not exist.\nCheck both dihedralangledata metadata and ResidueLinkage::FindRotatableDihedralsConnectingResidues.\nNote this is normal for a sialic acid with multiple 2-7, 2-8 and or 2-9 linkages\n";
            gmml::log(__LINE__,__FILE__,gmml::WAR, message);
            //throw std::runtime_error(message);
        }
    }
//    std::cout << "\n";
    return;
}


template<class residueT, class atomT>
void ResidueLinkage<residueT, atomT>::SetResidues(residueT* residue1, residueT* residue2)
{
    from_this_residue1_ = residue1;
    to_this_residue2_ = residue2;
}

template<class residueT, class atomT>
void ResidueLinkage<residueT, atomT>::SetConnectionAtoms(residueT* residue1, residueT* residue2)
{
    std::vector<atomT*> connecting_atoms;
    bool found = false;
    selection::FindAtomsConnectingResidues(residue1->GetAtoms().at(0), residue2, &connecting_atoms, &found);
    from_this_connection_atom1_ = connecting_atoms.at(0);
    to_this_connection_atom2_ = connecting_atoms.at(1);
}

template<class residueT, class atomT>
void ResidueLinkage<residueT, atomT>::SetConformerUsingMetadata(bool useRanges, int conformerNumber)
{
    for (auto &entry : RotatableDihedrals_)
    {
        entry.SetSpecificAngleEntryUsingMetadata(useRanges, conformerNumber);
    }
}

template<class residueT, class atomT>
unsigned long long  ResidueLinkage<residueT, atomT>::GenerateIndex()
{
    static unsigned long long s_ResidueLinkageIndex = 0; // static keyword means it is created only once and persists beyond scope of code block.
    return s_ResidueLinkageIndex++; // makes copy of s_AtomIndex, increments the real s_AtomIndex, then returns the value in the copy
}


} //namespace
#endif // GMML_INCLUDES_GEOMETRYTOPOLOGY_RESIDUELINKAGES_RESIDUE_LINKAGE_HPP
