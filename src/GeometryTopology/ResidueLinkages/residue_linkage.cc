#include <random>
#include <algorithm> // reverse
#include "includes/GeometryTopology/ResidueLinkages/residue_linkage.hpp"
#include "includes/MolecularModeling/overlaps.hpp"
#include "includes/External_Libraries/PCG/pcg_extras.h"
#include "includes/External_Libraries/PCG/pcg_random.h"
#include "includes/MolecularMetadata/GLYCAM/glycam06residuecodes.hpp" // For lookup in GetName function
#include "includes/CodeUtils/logging.hpp"

// Seed with a real random value, if available
static pcg_extras::seed_seq_from<std::random_device> seed_source;
// Make a random number engine
static pcg32 rng(seed_source);

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
Residue_linkage::Residue_linkage(Residue* nonReducingResidue1, Residue* reducingResidue2, bool reverseAtomsThatMove)
{
    isExtraAtoms_ = false;
    this->InitializeClass(nonReducingResidue1, reducingResidue2, reverseAtomsThatMove);
}

Residue_linkage::Residue_linkage(Residue* nonReducingResidue1, Residue* reducingResidue2, AtomVector alsoMovingAtoms,
                                 bool reverseAtomsThatMove)
{ // Order of calling functions is important!
    this->AddExtraAtomsThatMove(alsoMovingAtoms);
    this->InitializeClass(nonReducingResidue1, reducingResidue2, reverseAtomsThatMove);
}

//////////////////////////////////////////////////////////
//                       ACCESSOR                       //
//////////////////////////////////////////////////////////

ResidueVector Residue_linkage::GetResidues() const
{
    ResidueVector residues {from_this_residue1_, to_this_residue2_};
    return residues;
}

bool Residue_linkage::GetIfReversedAtomsThatMove() const
{
    return reverseAtomsThatMove_;
}

std::vector<Rotatable_dihedral> Residue_linkage::GetRotatableDihedralsWithMultipleRotamers() const
{
    std::vector<Rotatable_dihedral> returningDihedrals;
    for (auto& entry : this->GetRotatableDihedrals())
    {
        if (entry.GetNumberOfRotamers() > 1)
        {
            returningDihedrals.push_back(entry);
        }
    }
    return returningDihedrals;
}

std::vector<Rotatable_dihedral> Residue_linkage::GetRotatableDihedrals() const
{
    if (rotatable_dihedrals_.empty())
    {
        std::stringstream ss;
        ss << "Error: rotatable_dihedrals in this linkage is empty: " << from_this_residue1_->GetId() << "-"
           << to_this_residue2_->GetId() << std::endl;
        gmml::log(__LINE__, __FILE__, gmml::ERR, ss.str());
        throw std::runtime_error(ss.str());
    }
    return rotatable_dihedrals_;
}

int Residue_linkage::GetNumberOfShapes(
    const bool likelyShapesOnly) const // Can have conformers (sets of rotamers) or permutations of rotamers
{
    int numberOfShapes = 1;
    if ((rotatable_dihedrals_.empty()) || (rotatable_dihedrals_.at(0).GetMetadata().empty()))
    {
        return numberOfShapes;
    }
    if (rotatable_dihedrals_.at(0).GetMetadata().at(0).rotamer_type_.compare("permutation") == 0)
    {
        for (auto& entry : rotatable_dihedrals_)
        {
            numberOfShapes = (numberOfShapes * entry.GetNumberOfRotamers(likelyShapesOnly));
        }
    }
    else if (rotatable_dihedrals_.at(0).GetMetadata().at(0).rotamer_type_.compare("conformer") == 0)
    { // Conformer should mean that each dihedral will have the same number of metadata entries.
        // numberOfShapes = rotatable_dihedrals_.size(); // This was correct for ASN for the wrong reason. 4 conformers
        // and 4 dihedrals...
        numberOfShapes = rotatable_dihedrals_.at(0).GetNumberOfRotamers(likelyShapesOnly);
    }
    return numberOfShapes;
}

Residue* Residue_linkage::GetFromThisResidue1() const
{
    return from_this_residue1_;
}

Residue* Residue_linkage::GetToThisResidue2() const
{
    return to_this_residue2_;
}

Atom* Residue_linkage::GetFromThisConnectionAtom1() const
{
    return from_this_connection_atom1_;
}

Atom* Residue_linkage::GetToThisConnectionAtom2() const
{
    return to_this_connection_atom2_;
}

bool Residue_linkage::CheckIfConformer() const
{
    if (rotatable_dihedrals_.empty())
    {
        std::string errorMessage = "Error in Residue_linkage::checkIfConformer as rotatable_dihedrals_.empty()\n";
        gmml::log(__LINE__, __FILE__, gmml::ERR, errorMessage);
        throw std::runtime_error(errorMessage);
    }
    else if (rotatable_dihedrals_.at(0).GetMetadata().empty())
    {
        std::string errorMessage =
            "Error in Residue_linkage::checkIfConformer as rotatable_dihedrals_.at(0).GetMetadata().empty()\n";
        gmml::log(__LINE__, __FILE__, gmml::ERR, errorMessage);
        throw std::runtime_error(errorMessage);
    }
    else
    {
        if (rotatable_dihedrals_.at(0).GetMetadata().at(0).rotamer_type_.compare("permutation") == 0)
        {
            return false;
        }
        else
        {
            return true;
        }
    }
    return false; // Default to shut up the compiler. Shut up compiler gawd.
}

bool Residue_linkage::GetIfExtraAtoms() const
{
    return isExtraAtoms_;
}

AtomVector Residue_linkage::GetExtraAtoms()
{
    return extraAtomsThatMove_;
}

unsigned long long Residue_linkage::GetIndex() const
{
    return index_;
}

std::string Residue_linkage::GetName() const
{
    if (!name_.empty())
    {
        return name_;
    }
    return this->DetermineLinkageNameFromResidueNames();
}

std::string Residue_linkage::DetermineLinkageNameFromResidueNames() const
{
    gmml::MolecularMetadata::GLYCAM::Glycam06ResidueNamesToCodesLookupContainer nameLookup;
    std::string residue1Name = nameLookup.GetResidueForCode(this->GetFromThisResidue1()->GetName());
    std::string residue2Name = nameLookup.GetResidueForCode(this->GetToThisResidue2()->GetName());
    // std::cout << this->GetFromThisConnectionAtom1()->GetName() << std::endl;
    // std::cout << this->GetToThisConnectionAtom2()->GetName() << std::endl;
    std::string atom1Name    = this->GetFromThisConnectionAtom1()->GetName();
    std::string atom2Name    = this->GetToThisConnectionAtom2()->GetName();
    char link1               = *atom1Name.rbegin(); //
    char link2               = *atom2Name.rbegin(); // Messy for Acetyl.
    std::stringstream linkageName;
    linkageName << residue1Name << link1 << "-" << link2 << residue2Name;
    return linkageName.str();
}

int Residue_linkage::GetNumberOfRotatableDihedrals() const
{
    return rotatable_dihedrals_.size();
}

//////////////////////////////////////////////////////////
//                       MUTATOR                        //
//////////////////////////////////////////////////////////

void Residue_linkage::SetRotatableDihedrals(std::vector<Rotatable_dihedral> rotatableDihedrals)
{
    rotatable_dihedrals_ = rotatableDihedrals;
}

void Residue_linkage::SetIfReversedAtomsThatMove(bool reversedAtomsThatMove)
{
    reverseAtomsThatMove_ = reversedAtomsThatMove;
}

void Residue_linkage::AddExtraAtomsThatMove(AtomVector extraAtoms)
{
    extraAtomsThatMove_ = extraAtoms;
    isExtraAtoms_       = true;
}

//////////////////////////////////////////////////////////
//                       FUNCTIONS                      //
//////////////////////////////////////////////////////////
void Residue_linkage::SetDefaultShapeUsingMetadata()
{
    for (auto& entry : rotatable_dihedrals_)
    {
        entry.SetSpecificAngleEntryUsingMetadata(); // Default is first entry
    }
}

void Residue_linkage::SetRandomShapeUsingMetadata(bool useRanges)
{
    if (rotatable_dihedrals_.at(0).GetMetadata().at(0).rotamer_type_.compare("permutation") == 0)
    {
        for (auto& entry : rotatable_dihedrals_)
        {
            entry.SetRandomAngleEntryUsingMetadata(useRanges);
        }
    }
    else if (rotatable_dihedrals_.at(0).GetMetadata().at(0).rotamer_type_.compare("conformer") == 0)
    {
        int numberOfConformers = rotatable_dihedrals_.at(0).GetMetadata().size();
        std::uniform_int_distribution<> distr(0, (numberOfConformers - 1)); // define the range
        int randomlySelectedConformerNumber = distr(rng);
        for (auto& entry : rotatable_dihedrals_)
        {
            entry.SetSpecificAngleEntryUsingMetadata(useRanges, randomlySelectedConformerNumber);
        }
    }
}

// This is dangerous if used by non-conformer linkage. Change it to check.
void Residue_linkage::SetSpecificShapeUsingMetadata(int shapeNumber, bool useRanges)
{
    for (auto& entry : rotatable_dihedrals_)
    {
        entry.SetSpecificAngleEntryUsingMetadata(useRanges, shapeNumber);
    }
}

void Residue_linkage::SetSpecificShape(std::string dihedralName, std::string selectedRotamer)
{
    for (auto& rotatable_dihedral : rotatable_dihedrals_)
    {
        // This will call rotatable_dihedrals that don't have dihedralName (phi,psi), and nothing will happen. Hmmm.
        if (rotatable_dihedral.SetSpecificShape(dihedralName, selectedRotamer))
        {
            return; // Return once you manage to set a shape.
        }
    }
    std::string errorMessage = "Did not set " + dihedralName + " to " + selectedRotamer +
                               " as requested in Residue_linkage::SetSpecificShape()";
    gmml::log(__LINE__, __FILE__, gmml::ERR, errorMessage);
    throw std::runtime_error(errorMessage);
}

void Residue_linkage::SetCustomDihedralAngles(std::vector<double> dihedral_angles)
{
    if (dihedral_angles.size() == rotatable_dihedrals_.size())
    {
        std::vector<double>::iterator dihedral_angle_iterator = dihedral_angles.begin();
        for (auto& rotatable_dihedral : rotatable_dihedrals_)
        {
            rotatable_dihedral.SetDihedralAngle(*dihedral_angle_iterator);
            ++dihedral_angle_iterator;
        }
    }
    else
    {
        std::stringstream ss;
        ss << "ERROR attempted to set " << dihedral_angles.size() << " dihedrals for a residue linkage with "
           << rotatable_dihedrals_.size()
           << " rotatable dihedrals\n These numbers need to match for this function to work\nSomething has gone wrong.";
        gmml::log(__LINE__, __FILE__, gmml::ERR, ss.str());
        throw std::runtime_error(ss.str());
    }
    return;
}

void Residue_linkage::SetShapeToPrevious()
{
    for (std::vector<Rotatable_dihedral>::iterator rotatable_dihedral = rotatable_dihedrals_.begin();
         rotatable_dihedral != rotatable_dihedrals_.end(); ++rotatable_dihedral)
    {
        rotatable_dihedral->SetDihedralAngleToPrevious();
    }
    return;
}

void Residue_linkage::SetRandomDihedralAngles()
{
    for (std::vector<Rotatable_dihedral>::iterator rotatable_dihedral = rotatable_dihedrals_.begin();
         rotatable_dihedral != rotatable_dihedrals_.end(); ++rotatable_dihedral)
    {
        rotatable_dihedral->RandomizeDihedralAngle();
    }
    return;
}

void Residue_linkage::DetermineAtomsThatMove()
{
    for (std::vector<Rotatable_dihedral>::iterator rotatable_dihedral = rotatable_dihedrals_.begin();
         rotatable_dihedral != rotatable_dihedrals_.end(); ++rotatable_dihedral)
    {
        rotatable_dihedral->DetermineAtomsThatMove();
    }
    return;
}

void Residue_linkage::SimpleWiggle(AtomVector& overlapAtomSet1, AtomVector& overlapAtomSet2, const int angleIncrement)
{
    for (auto& rotatable_dihedral : this->GetRotatableDihedrals())
    {
        rotatable_dihedral.WiggleUsingAllRotamers(overlapAtomSet1, overlapAtomSet2, angleIncrement);
    }
}

void Residue_linkage::SimpleWiggleCurrentRotamers(AtomVector& overlapAtomSet1, AtomVector& overlapAtomSet2,
                                                  const int angleIncrement)
{
    for (auto& rotatable_dihedral : this->GetRotatableDihedrals())
    {
        //        std::cout << "SimpleWiggling current rotamer: " << rotatable_dihedral.Print() << "\n";
        rotatable_dihedral.WiggleWithinCurrentRotamer(overlapAtomSet1, overlapAtomSet2, angleIncrement);
    }
}

void Residue_linkage::SetIndex(unsigned long long index)
{
    index_ = index;
}

//////////////////////////////////////////////////////////
//                       DISPLAY FUNCTION               //
//////////////////////////////////////////////////////////
std::string Residue_linkage::Print() const
{
    std::stringstream ss;
    ss << "Residue_linkage Index: " << this->GetIndex() << ", Name: " << this->GetName()
       << ", NumberOfShapes: " << this->GetNumberOfShapes() << ", ids: " << this->GetFromThisResidue1()->GetId() << "@"
       << this->GetFromThisConnectionAtom1()->GetName() << " -- " << this->GetToThisResidue2()->GetId() << "@"
       << this->GetToThisConnectionAtom2()->GetName() << "\n";
    gmml::log(__LINE__, __FILE__, gmml::INF, ss.str());
    for (auto& rotatableDihedral : this->GetRotatableDihedrals())
    {
        ss << rotatableDihedral.Print();
    }
    return ss.str();
}

//////////////////////////////////////////////////////////
//                    PRIVATE FUNCTIONS                 //
//////////////////////////////////////////////////////////
void Residue_linkage::InitializeClass(Residue* from_this_residue1, Residue* to_this_residue2, bool reverseAtomsThatMove)
{
    // set local debug flag
    int local_debug = -1;

    this->SetResidues(from_this_residue1, to_this_residue2);
    this->SetIfReversedAtomsThatMove(reverseAtomsThatMove);
    this->SetConnectionAtoms(from_this_residue1_, to_this_residue2_);
    if (local_debug > 0)
    {
        // std::cout << "Maybe Finding connection between " << from_this_residue1->GetId() << " :: " <<
        // to_this_residue2->GetId() << std::endl;
        gmml::log(__LINE__, __FILE__, gmml::INF,
                  "Maybe Finding connection between " + from_this_residue1->GetId() +
                      " :: " + to_this_residue2->GetId());
    }
    if (this->CheckIfViableLinkage())
    {
        //        std::cout << "Finding connection between " << from_this_residue1->GetId() << " :: " <<
        //        to_this_residue2->GetId() << std::endl; std::cout << "Connection atoms are from: " <<
        //        from_this_connection_atom1_->GetId() << " to " << to_this_connection_atom2_->GetId() << std::endl;
        if (local_debug > 0)
        {
            gmml::log(__LINE__, __FILE__, gmml::INF,
                      "Finding connection between " + from_this_residue1->GetId() + " :: " + to_this_residue2->GetId());
            gmml::log(__LINE__, __FILE__, gmml::INF,
                      "Connection atoms are from: " + from_this_connection_atom1_->GetId() + " to " +
                          to_this_connection_atom2_->GetId());
        }
        rotatable_dihedrals_ =
            this->FindRotatableDihedralsConnectingResidues(from_this_connection_atom1_, to_this_connection_atom2_);
        //        std::cout << "Finding metadata for " << from_this_residue1->GetId() << " :: " <<
        //        to_this_residue2->GetId() << std::endl;
        if (local_debug > 0)
        {
            gmml::log(__LINE__, __FILE__, gmml::INF,
                      "Finding metadata for " + from_this_residue1->GetId() + " :: " + to_this_residue2->GetId());
        }
        gmml::MolecularMetadata::GLYCAM::DihedralAngleDataVector metadata = this->FindMetadata();
        //  std::cout << "Metadata found:\n";
        //  for (auto &dihedralAngleData : metadata)
        //  {
        //      std::cout << dihedralAngleData.print() << std::endl;
        //  }
        if (local_debug > 0)
        {
            gmml::log(__LINE__, __FILE__, gmml::INF, "Metadata found:");
            for (auto& dihedralAngleData : metadata)
            {
                gmml::log(__LINE__, __FILE__, gmml::INF, dihedralAngleData.print());
            }
        }
        this->AddMetadataToRotatableDihedrals(metadata);
    }
    this->SetIndex(this->GenerateIndex());
}

bool Residue_linkage::CheckIfViableLinkage() const
{
    for (auto& residue : this->GetResidues())
    {
        if (residue->GetAtoms().size() <= 1)
        { // If either linkage has only 1 atom, return false. Should not set dihedral.
            return false;
        }
    }
    return true;
}

std::vector<Rotatable_dihedral>
Residue_linkage::FindRotatableDihedralsConnectingResidues(Atom* from_this_connection_atom1,
                                                          Atom* to_this_connection_atom2)
{
    // Going to ignore tags etc.
    // Given two residues that are connected. Find connecting atoms.
    // Search neighbors other than connected atom. Ie search out in both directions, but remain within same residue.
    // Warning, residue may have fused cycles!
    // Will fail for non-protein residues without cycles. As don't have a non-rotatable bond to anchor from. Can code
    // that later (and deal with branches from these residues).
    //    std::cout << "Finding rot bonds for " << from_this_connection_atom1->GetResidue()->GetId() << " and " <<
    //    to_this_connection_atom2->GetResidue()->GetId() << "\n";

    AtomVector from_this_residue1_cycle_points = selection::FindCyclePoints(from_this_connection_atom1);
    //    std::cout << "Moving onto second residue.\n";
    AtomVector to_this_residue2_cycle_points   = selection::FindCyclePoints(to_this_connection_atom2);
    // Need to reverse one of these, so when concatenated, they are ordered ok. This might not be ok.
    // std::reverse(to_this_residue2_cycle_points.begin(), to_this_residue2_cycle_points.end());
    std::reverse(from_this_residue1_cycle_points.begin(), from_this_residue1_cycle_points.end());
    // Now concatenate:
    from_this_residue1_cycle_points.insert(from_this_residue1_cycle_points.end(), to_this_residue2_cycle_points.begin(),
                                           to_this_residue2_cycle_points.end());
    // Now that have a list of rotation points. Split into pairs and find rotatable bonds between them
    bool found                  = false;
    AtomVector connecting_atoms = {from_this_connection_atom1, to_this_connection_atom2};
    //     std::cout << "cycle point atoms are:\n";
    //     for(auto & atom : from_this_residue1_cycle_points)
    //     {
    //         std::cout << atom->GetId() << "\n";
    //     }
    //     std::cout << "\n";
    std::vector<Rotatable_dihedral> rotatableDihedralsInBranches;
    for (int i = 0; i < from_this_residue1_cycle_points.size(); i = i + 2)
    {
        Atom* cycle_point1 = from_this_residue1_cycle_points.at(i);
        Atom* cycle_point2 = from_this_residue1_cycle_points.at(i + 1);

        found = false;
        connecting_atoms.clear();
        //        std::cout << "Finding Path between:" << cycle_point1->GetId() << " and " << cycle_point2->GetId() <<
        //        "\n";
        selection::FindPathBetweenTwoAtoms(cycle_point1, cycle_point2, &connecting_atoms, &found);
        selection::ClearAtomDescriptions(cycle_point1->GetResidue());
        selection::ClearAtomDescriptions(cycle_point2->GetResidue());
        // Find neighboring atoms needed to define dihedral. Pass in connecting atoms so don't find any of those.
        Atom* neighbor1 = selection::FindCyclePointNeighbor(connecting_atoms, cycle_point1);
        Atom* neighbor2 = selection::FindCyclePointNeighbor(connecting_atoms, cycle_point2);
        // Insert these neighbors into list of connecting atoms, at beginning and end of vector.
        // connecting_atoms gets populated as it falls out, so list is reversed from what you'd expect
        std::reverse(connecting_atoms.begin(), connecting_atoms.end());
        connecting_atoms.insert(connecting_atoms.begin(), neighbor1);
        connecting_atoms.push_back(neighbor2);

        //         std::cout << "Updated Path between:\n " << cycle_point1->GetId() << " and " << cycle_point2->GetId()
        //         << "\n"; for (const auto& atom : connecting_atoms)
        //             std::cout << atom->GetId() << "\n";
        //         std::cout << "\n";
        selection::ClearAtomDescriptions(cycle_point1->GetResidue());
        selection::ClearAtomDescriptions(cycle_point2->GetResidue());
        // This mess was made to address the branching in 2-7 and 2-8 linkages.
        // These branches are long enough that they need default torsions set.
        if (connecting_atoms.size() > 4) // Otherwise there are no torsions
        {                                // Only viable linkages. Throw if not >4?
            for (AtomVector::iterator it = connecting_atoms.begin() + 1; it != connecting_atoms.end() - 1; ++it)
            { //
                Atom* connectionAtom = (*it);
                if ((connectionAtom != cycle_point1) && (connectionAtom != cycle_point2))
                {
                    for (auto& neighbor : connectionAtom->GetNode()->GetNodeNeighbors()) // Need an interator
                    {
                        if (std::find(connecting_atoms.begin(), connecting_atoms.end(), neighbor) ==
                            connecting_atoms.end()) // if not in the vector
                        {
                            selection::Branch branch(connectionAtom);
                            selection::FindEndsOfBranchesFromLinkageAtom(neighbor, connectionAtom, &branch);
                            if (branch.IsBranchFound())
                            {
                                found = false;
                                AtomVector foundPath;
                                // This fills in foundPath:
                                selection::FindPathBetweenTwoAtoms(branch.GetRoot(), branch.GetEnd(), &foundPath,
                                                                   &found);
                                Atom* neighbor = selection::FindCyclePointNeighbor(foundPath, branch.GetRoot());
                                // foundPath.insert(foundPath.begin(), neighbor);
                                foundPath.push_back(neighbor);
                                //                                 std::cout << "Found atoms:\n";
                                //                                 for (auto &atom: foundPath)
                                //                                     std::cout << atom->GetId() << "\n";
                                std::vector<Rotatable_dihedral> temp =
                                    this->SplitAtomVectorIntoRotatableDihedrals(foundPath);
                                rotatableDihedralsInBranches.insert(rotatableDihedralsInBranches.end(), temp.begin(),
                                                                    temp.end());
                            }
                        }
                    }
                }
            }
        } // End dealing with branching linkages
          //         std::cout << "These are the assigned branched rotatable_dihedrals:\n";
          //         for (auto &dihedral : rotatableDihedralsInBranches)
          //             dihedral.Print();
    }
    std::vector<Rotatable_dihedral> rotatable_dihedrals = this->SplitAtomVectorIntoRotatableDihedrals(connecting_atoms);
    // Add any linkage branches (in 2-7 and 2-8) to the rest.
    rotatable_dihedrals.insert(rotatable_dihedrals.end(), rotatableDihedralsInBranches.begin(),
                               rotatableDihedralsInBranches.end());
    //     std::cout << "These are the assigned rotatable_dihedrals:\n";
    //     for (auto &dihedral : rotatable_dihedrals)
    //         dihedral.Print();
    return rotatable_dihedrals;
}

std::vector<Rotatable_dihedral> Residue_linkage::SplitAtomVectorIntoRotatableDihedrals(AtomVector atoms)
{
    // Ok looking for sets of four atoms, but shifting along vector by one atom for each dihedral.
    //  So four atoms will make one rotatable bond, five will make two bonds, six will make three etc.
    std::vector<Rotatable_dihedral> rotatable_dihedrals_generated;
    if (atoms.size() < 4)
    {
        std::stringstream ss;
        ss << "ERROR in Residue_linkage::SplitAtomVectorIntoRotatableDihedrals, not enough atoms in atom vector: "
           << atoms.size() << "\n";
        ss << "This should be 4 or something is very wrong\n";
        gmml::log(__LINE__, __FILE__, gmml::ERR, ss.str());
        throw std::runtime_error(ss.str());
    }
    else
    {
        for (AtomVector::iterator it1 = atoms.begin(); it1 != (atoms.end() - 3); ++it1)
        {
            Atom* atom1 = *it1;
            Atom* atom2 = *(it1 + 1);
            Atom* atom3 = *(it1 + 2);
            Atom* atom4 = *(it1 + 3);
            rotatable_dihedrals_generated.emplace_back(atom1, atom2, atom3, atom4, this->GetIfReversedAtomsThatMove());
        }
    }
    return rotatable_dihedrals_generated;
}

gmml::MolecularMetadata::GLYCAM::DihedralAngleDataVector Residue_linkage::FindMetadata() const
{
    gmml::MolecularMetadata::GLYCAM::DihedralAngleDataContainer DihedralAngleMetadata;
    gmml::MolecularMetadata::GLYCAM::DihedralAngleDataVector matching_entries =
        DihedralAngleMetadata.GetEntriesForLinkage(
            this->GetFromThisConnectionAtom1()->GetName(), this->GetFromThisResidue1()->GetName(),
            this->GetToThisConnectionAtom2()->GetName(), this->GetToThisResidue2()->GetName());
    // std::cout << "Found these " << matching_entries.size() << " entries:\n";
    //  for (const auto& entry : matching_entries)
    //  {
    //      std::cout << entry.index_ << " : " << entry.atom1_ << ", " << entry.atom2_ << ", " << entry.atom3_ << ", "
    //      << entry.atom4_ << ", " << entry.default_angle_value_ << "\n";
    //  }
    if (matching_entries.empty())
    {
        std::stringstream ss;
        ss << "No Metadata entries found for connection between " << this->GetFromThisConnectionAtom1()->GetId()
           << " and " << this->GetToThisConnectionAtom2()->GetId() << "\n";
        ss << "Note that order should be reducing atom - anomeric atom\n";
        gmml::log(__LINE__, __FILE__, gmml::ERR, ss.str());
        throw std::runtime_error(ss.str());
    }
    return matching_entries;
}

// This is funky. Either better checking or a redesign?
// Update March 2022: No need to throw in case of branched 2-7 or 2-8 when there is something on the 2-9 so the linkage
// finder algo doesn't detect those as rotatable by this linkage, but the meta data for those dihedrals is still found.
// I think separating out the adding metadata information into a separate step from the linkage finding step was a
// mistake. Also trying to make it so generic that it doesn't need to look at atom names is also hurting when trying to
// figure out errors.. I would like to throw when it's not a 2-7/2-8 and it's found too much metadata.
void Residue_linkage::AddMetadataToRotatableDihedrals(gmml::MolecularMetadata::GLYCAM::DihedralAngleDataVector metadata)
{
    // Adding another sanity check to this insanity
    if (rotatable_dihedrals_.size() > metadata.size())
    {
        std::stringstream ss;
        ss << "Problem with the metadata found in gmml for this linkage. The number of identified rotatable dihedrals: "
           << rotatable_dihedrals_.size() << " is greater than the number of metadata items: " << metadata.size()
           << " found for this linkage:\n"
           << this->Print();
        gmml::log(__LINE__, __FILE__, gmml::WAR, ss.str());
        throw std::runtime_error(ss.str());
    }
    for (auto& rotatableDihedral : rotatable_dihedrals_)
    {
        rotatableDihedral.ClearMetadata(); // First clear any metadata already in place.
    }
    for (const auto& entry : metadata)
    {
        int vector_position = (entry.number_of_bonds_from_anomeric_carbon_ - 1); // vectors start at 0.
        //        std::cout << "Adding to position: "<< vector_position << " in vector of size: " <<
        //        rotatable_dihedrals_.size() << std::endl;
        if (vector_position < rotatable_dihedrals_.size())
        {
            rotatable_dihedrals_.at(vector_position).AddMetadata(entry);
            //            std::cout << "Set " << entry.dihedral_angle_name_ << " meta data for "
            //                      << rotatable_dihedrals_.at(vector_position).GetAtoms().at(0)->GetId() << ", "
            //                      << rotatable_dihedrals_.at(vector_position).GetAtoms().at(1)->GetId() << ", "
            //                      << rotatable_dihedrals_.at(vector_position).GetAtoms().at(2)->GetId() << ", "
            //                      << rotatable_dihedrals_.at(vector_position).GetAtoms().at(3)->GetId() << ". ";
            //            std::cout << "Index: " << entry.index_ << ", angle: " << entry.default_angle_value_ << "\n";
            // rotatable_dihedrals_.at(vector_position).Print();
        }
        else
        {
            std::string message =
                "Tried to add metadata to a rotatable bond that does not exist.\nCheck both dihedralangledata metadata "
                "and Residue_linkage::FindRotatableDihedralsConnectingResidues.\nNote this is normal for a sialic acid "
                "with multiple 2-7, 2-8 and or 2-9 linkages\n";
            gmml::log(__LINE__, __FILE__, gmml::WAR, message);
            // throw std::runtime_error(message);
        }
    }
    //    std::cout << "\n";
    return;
}

void Residue_linkage::SetResidues(Residue* residue1, Residue* residue2)
{
    from_this_residue1_ = residue1;
    to_this_residue2_   = residue2;
}

void Residue_linkage::SetConnectionAtoms(Residue* residue1, Residue* residue2)
{
    AtomVector connecting_atoms;
    bool found = false;
    selection::FindAtomsConnectingResidues(residue1->GetAtoms().at(0), residue2, &connecting_atoms, &found);
    from_this_connection_atom1_ = connecting_atoms.at(0);
    to_this_connection_atom2_   = connecting_atoms.at(1);
}

void Residue_linkage::SetConformerUsingMetadata(bool useRanges, int conformerNumber)
{
    for (auto& entry : rotatable_dihedrals_)
    {
        entry.SetSpecificAngleEntryUsingMetadata(useRanges, conformerNumber);
    }
}

unsigned long long Residue_linkage::GenerateIndex()
{
    static unsigned long long s_ResidueLinkageIndex =
        0; // static keyword means it is created only once and persists beyond scope of code block.
    return s_ResidueLinkageIndex++; // makes copy of s_AtomIndex, increments the real s_AtomIndex, then returns the
                                    // value in the copy
}

void Residue_linkage::SetName(std::string name)
{
    name_ = name;
}
