#include "includes/CentralDataStructure/Shapers/residueLinkage.hpp"
#include "includes/CentralDataStructure/Selections/shaperSelections.hpp"
#include "includes/MolecularMetadata/GLYCAM/glycam06residuecodes.hpp" // For lookup in GetName function

using cds::ResidueLinkage;
using cds::RotatableDihedral;


ResidueLinkage::ResidueLinkage(cds::Residue *nonReducingResidue1, cds::Residue *reducingResidue2, bool reverseAtomsThatMove)
{
    isExtraAtoms_ = false;
    this->InitializeClass(nonReducingResidue1, reducingResidue2, reverseAtomsThatMove);
}

ResidueLinkage::ResidueLinkage(cds::Residue *nonReducingResidue1, cds::Residue *reducingResidue2, std::vector<cds::Atom*> alsoMovingAtoms, bool reverseAtomsThatMove)
{ // Order of calling functions is important!
    this->AddExtraAtomsThatMove(alsoMovingAtoms);
    this->InitializeClass(nonReducingResidue1, reducingResidue2, reverseAtomsThatMove);
}

std::vector<cds::Residue*> ResidueLinkage::GetResidues() const
{
    std::vector<cds::Residue*> residues {from_this_residue1_, to_this_residue2_};
    return residues;
}

std::vector<RotatableDihedral> ResidueLinkage::GetRotatableDihedralsWithMultipleRotamers() const
{
    std::vector<RotatableDihedral> returningDihedrals;
    for (auto &entry : this->GetRotatableDihedrals())
    {
        if (entry.GetNumberOfRotamers() > 1)
        {
            returningDihedrals.push_back(entry);
        }
    }
    return returningDihedrals;
}

std::vector<RotatableDihedral> ResidueLinkage::GetRotatableDihedrals() const
{
    if (rotatableDihedrals_.empty())
    {
        std::stringstream ss;
        ss << "Error: RotatableDihedrals in this linkage is empty: " << from_this_residue1_->getStringId() << "-" << to_this_residue2_->getStringId() << std::endl;
        gmml::log(__LINE__,__FILE__,gmml::ERR, ss.str());
        throw std::runtime_error(ss.str());
    }
    return rotatableDihedrals_;
}

int ResidueLinkage::GetNumberOfShapes( const bool likelyShapesOnly) const// Can have conformers (sets of rotamers) or permutations of rotamers
{
    int numberOfShapes = 1;
    if ( (rotatableDihedrals_.empty()) || (rotatableDihedrals_.at(0).GetMetadata().empty()) )
    {
        return numberOfShapes;
    }
    if (rotatableDihedrals_.at(0).GetMetadata().at(0).rotamer_type_.compare("permutation")==0)
    {
        for (auto &entry : rotatableDihedrals_)
        {
            numberOfShapes = (numberOfShapes * entry.GetNumberOfRotamers(likelyShapesOnly));
        }
    }
    else if (rotatableDihedrals_.at(0).GetMetadata().at(0).rotamer_type_.compare("conformer")==0)
    { // Conformer should mean that each dihedral will have the same number of metadata entries.
        //numberOfShapes = RotatableDihedrals_.size(); // This was correct for ASN for the wrong reason. 4 conformers and 4 dihedrals...
        numberOfShapes = rotatableDihedrals_.at(0).GetNumberOfRotamers(likelyShapesOnly);
    }
    return numberOfShapes;
}

bool ResidueLinkage::CheckIfConformer() const
{
    if (rotatableDihedrals_.empty())
    {
        std::string errorMessage = "Error in ResidueLinkage::checkIfConformer as RotatableDihedrals_.empty()\n";
        gmml::log(__LINE__,__FILE__,gmml::ERR, errorMessage);
        throw std::runtime_error(errorMessage);
    }
    else if (rotatableDihedrals_.at(0).GetMetadata().empty())
    {
        std::string errorMessage = "Error in ResidueLinkage::checkIfConformer as RotatableDihedrals_.at(0).GetMetadata().empty()\n";
        gmml::log(__LINE__,__FILE__,gmml::ERR, errorMessage);
        throw std::runtime_error(errorMessage);
    }
    else
    {
        if (rotatableDihedrals_.at(0).GetMetadata().at(0).rotamer_type_.compare("permutation")==0)
            return false;
        else
            return true;
    }
    return false; //Default to shut up the compiler. Shut up compiler gawd.
}

std::string ResidueLinkage::GetName() const
{
    if (!name_.empty())
    {
        return name_;
    }
    return this->DetermineLinkageNameFromResidueNames();
}

std::string ResidueLinkage::DetermineLinkageNameFromResidueNames() const
{
    gmml::MolecularMetadata::GLYCAM::Glycam06ResidueNamesToCodesLookupContainer nameLookup;
    std::string residue1Name = nameLookup.GetResidueForCode(this->GetFromThisResidue1()->getName());
    std::string residue2Name = nameLookup.GetResidueForCode(this->GetToThisResidue2()->getName());
   // std::cout << this->GetFromThisConnectionAtom1()->GetName() << std::endl;
   // std::cout << this->GetToThisConnectionAtom2()->GetName() << std::endl;
    std::string atom1Name = this->GetFromThisConnectionAtom1()->getName();
    std::string atom2Name = this->GetToThisConnectionAtom2()->getName();
    char link1 = *atom1Name.rbegin(); //
    char link2 = *atom2Name.rbegin(); // Messy for Acetyl.
    std::stringstream linkageName;
    linkageName << residue1Name << link1 << "-" << link2 << residue2Name;
    return linkageName.str();
}

unsigned long int ResidueLinkage::GetNumberOfRotatableDihedrals() const
{
   return rotatableDihedrals_.size();
}
//////////////////////////////////////////////////////////
//                       MUTATOR                        //
//////////////////////////////////////////////////////////

void ResidueLinkage::AddExtraAtomsThatMove(std::vector<cds::Atom*> extraAtoms)
{
    extraAtomsThatMove_ = extraAtoms;
    isExtraAtoms_ = true;
}

void ResidueLinkage::SetDefaultShapeUsingMetadata()
{
    for (auto &entry : rotatableDihedrals_)
    {
        entry.SetSpecificAngleEntryUsingMetadata(); // Default is first entry
    }
}

void ResidueLinkage::SetRandomShapeUsingMetadata(bool useRanges)
{
    if (rotatableDihedrals_.at(0).GetMetadata().at(0).rotamer_type_.compare("permutation")==0)
    {
        for (auto &entry : rotatableDihedrals_)
        {
            entry.SetRandomAngleEntryUsingMetadata(useRanges);
        }
    }
    else if (rotatableDihedrals_.at(0).GetMetadata().at(0).rotamer_type_.compare("conformer")==0)
    {
        int numberOfConformers = rotatableDihedrals_.at(0).GetMetadata().size();
        std::uniform_int_distribution<> distr(0, (numberOfConformers - 1)); // define the range
        int randomlySelectedConformerNumber = distr(rng);
        for (auto &entry : rotatableDihedrals_)
        {
            entry.SetSpecificAngleEntryUsingMetadata(useRanges, randomlySelectedConformerNumber);
        }
    }
}

void ResidueLinkage::SetSpecificShapeUsingMetadata(int shapeNumber, bool useRanges)
{
    for (auto &entry : rotatableDihedrals_)
    {
        entry.SetSpecificAngleEntryUsingMetadata(useRanges, shapeNumber);
    }
}

void ResidueLinkage::SetSpecificShape(std::string dihedralName, std::string selectedRotamer)
{
    for (auto &RotatableDihedral : rotatableDihedrals_)
    {
        // This will call RotatableDihedrals that don't have dihedralName (phi,psi), and nothing will happen. Hmmm.
        if (RotatableDihedral.SetSpecificShape(dihedralName, selectedRotamer))
            return; // Return once you manage to set a shape.
    }
    std::string errorMessage = "Did not set " + dihedralName + " to " + selectedRotamer + " as requested in ResidueLinkage::SetSpecificShape()";
    gmml::log(__LINE__,__FILE__,gmml::ERR, errorMessage);
    throw std::runtime_error(errorMessage);
}

void ResidueLinkage::SetCustomDihedralAngles(std::vector <double> dihedral_angles)
{
    if(dihedral_angles.size() == rotatableDihedrals_.size())
    {
        std::vector <double>::iterator dihedral_angle_iterator = dihedral_angles.begin();
        for (auto &RotatableDihedral : rotatableDihedrals_)
        {
            RotatableDihedral.SetDihedralAngle(*dihedral_angle_iterator);
            ++dihedral_angle_iterator;
        }
    }
    else
    {
        std::stringstream ss;
        ss << "ERROR attempted to set " << dihedral_angles.size() << " dihedrals for a residue linkage with " << rotatableDihedrals_.size() << " rotatable dihedrals\n These numbers need to match for this function to work\nSomething has gone wrong.";
        gmml::log(__LINE__,__FILE__,gmml::ERR, ss.str());
        throw std::runtime_error(ss.str());
    }
    return;
}

void ResidueLinkage::SetShapeToPrevious()
{
    for(typename std::vector<RotatableDihedral>::iterator RotatableDihedral = rotatableDihedrals_.begin(); RotatableDihedral != rotatableDihedrals_.end(); ++RotatableDihedral)
    {
        RotatableDihedral->SetDihedralAngleToPrevious();
    }
    return;
}

void ResidueLinkage::SetRandomDihedralAngles()
{
    for(typename std::vector<RotatableDihedral>::iterator RotatableDihedral = rotatableDihedrals_.begin(); RotatableDihedral != rotatableDihedrals_.end(); ++RotatableDihedral)
    {
        RotatableDihedral->RandomizeDihedralAngle();
    }
    return;
}

void ResidueLinkage::DetermineAtomsThatMove()
{
    for(typename std::vector<RotatableDihedral>::iterator RotatableDihedral = rotatableDihedrals_.begin(); RotatableDihedral != rotatableDihedrals_.end(); ++RotatableDihedral)
    {
        RotatableDihedral->DetermineAtomsThatMove();
    }
    return;
}

void ResidueLinkage::SimpleWiggle(std::vector<cds::Atom*>& overlapAtomSet1, std::vector<cds::Atom*>& overlapAtomSet2, const int angleIncrement)
{
    for(auto &RotatableDihedral : this->GetRotatableDihedrals())
    {
        RotatableDihedral.WiggleUsingAllRotamers(overlapAtomSet1, overlapAtomSet2, angleIncrement);
    }
}

void ResidueLinkage::SimpleWiggleCurrentRotamers(std::vector<cds::Atom*>& overlapAtomSet1, std::vector<cds::Atom*>& overlapAtomSet2, const int angleIncrement)
{
    for(auto &RotatableDihedral : this->GetRotatableDihedrals())
     {
//        std::cout << "SimpleWiggling current rotamer: " << RotatableDihedral.Print() << "\n";
        RotatableDihedral.WiggleWithinCurrentRotamer(overlapAtomSet1, overlapAtomSet2, angleIncrement);
     }
}

void ResidueLinkage::SimpleWiggleCurrentRotamers(std::vector<cds::Residue*>& overlapSet1, std::vector<cds::Residue*>& overlapSet2, const int angleIncrement)
{
    for(auto &RotatableDihedral : this->GetRotatableDihedrals())
     {
//        std::cout << "SimpleWiggling current rotamer: " << RotatableDihedral.Print() << "\n";
        RotatableDihedral.WiggleWithinCurrentRotamer(overlapSet1, overlapSet2, angleIncrement);
     }
}

std::string ResidueLinkage::Print() const
{
    std::stringstream ss;
    ss << "ResidueLinkage Index: " << this->GetIndex() << ", Name: " << this->GetName() << ", NumberOfShapes: " << this->GetNumberOfShapes()
              << ", ids: " << this->GetFromThisResidue1()->getStringId() << "@" << this->GetFromThisConnectionAtom1()->getName()
              << " -- " << this->GetToThisResidue2()->getStringId() << "@" << this->GetToThisConnectionAtom2()->getName() << "\n";
    gmml::log(__LINE__, __FILE__, gmml::INF, ss.str());
    for(auto & rotatableDihedral : this->GetRotatableDihedrals() )
    {
         ss << rotatableDihedral.Print();
    }
    return ss.str();
}

// PRIVATE
void ResidueLinkage::InitializeClass(cds::Residue* from_this_residue1, cds::Residue* to_this_residue2, bool reverseAtomsThatMove)
{
    //set local debug flag
    int local_debug = -1;
    gmml::log(__LINE__, __FILE__, gmml::INF, "Maybe Finding connection between " + from_this_residue1->getStringId() + " :: " + to_this_residue2->getStringId());
    this->SetResidues(from_this_residue1, to_this_residue2);
    this->SetIfReversedAtomsThatMove(reverseAtomsThatMove);
    this->SetConnectionAtoms(from_this_residue1_, to_this_residue2_);
    if(local_debug > 0)
    {
        // std::cout << "Maybe Finding connection between " << from_this_residue1->getId() << " :: " << to_this_residue2->getId() << std::endl;
        gmml::log(__LINE__, __FILE__, gmml::INF, "Maybe Finding connection between " + from_this_residue1->getStringId() + " :: " + to_this_residue2->getStringId());
    }
    if(this->CheckIfViableLinkage())
    {
//        std::cout << "Finding connection between " << from_this_residue1->getId() << " :: " << to_this_residue2->getId() << std::endl;
//        std::cout << "Connection atoms are from: " << from_this_connection_atom1_->getId() << " to " << to_this_connection_atom2_->getId() << std::endl;
        if(local_debug > 0)
        {
            gmml::log(__LINE__, __FILE__, gmml::INF, "Finding connection between " + from_this_residue1->getStringId() + " :: " + to_this_residue2->getStringId());
            gmml::log(__LINE__, __FILE__, gmml::INF, "Connection atoms are from: " + from_this_connection_atom1_->getId() + " to " + to_this_connection_atom2_->getId());
        }
        rotatableDihedrals_ = this->FindRotatableDihedralsConnectingResidues(from_this_connection_atom1_, to_this_connection_atom2_);
//        std::cout << "Finding metadata for " << from_this_residue1->getId() << " :: " << to_this_residue2->getId() << std::endl;
        if(local_debug > 0)
        {
            gmml::log(__LINE__, __FILE__, gmml::INF, "Finding metadata for " + from_this_residue1->getStringId() + " :: " + to_this_residue2->getStringId());
        }
        gmml::MolecularMetadata::GLYCAM::DihedralAngleDataVector metadata = this->FindMetadata();
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

bool ResidueLinkage::CheckIfViableLinkage() const
{
    for(auto &residue : this->GetResidues())
    {
        if (residue->getAtoms().size() <= 1)
        { // If either linkage has only 1 atom, return false. Should not set dihedral.
            return false;
        }
    }
    return true;
}

std::vector<RotatableDihedral> ResidueLinkage::FindRotatableDihedralsConnectingResidues(cds::Atom* from_this_connection_atom1, cds::Atom* to_this_connection_atom2)
{
    // Going to ignore tags etc.
    // Given two residues that are connected. Find connecting atoms.
    // Search neighbors other than connected atom. Ie search out in both directions, but remain within same residue.
    // Warning, residue may have fused cycles!
    // Will fail for non-protein residues without cycles. As don't have a non-rotatable bond to anchor from. Can code that later (and deal with branches from these residues).
//    std::cout << "Finding rot bonds for " << from_this_connection_atom1->GetResidue()->getId() << " and " << to_this_connection_atom2->GetResidue()->getId() << "\n";

    std::vector<cds::Atom*> from_this_residue1_cycle_points = cdsSelections::FindCyclePoints(from_this_connection_atom1, this->GetFromThisResidue1());
//    std::cout << "Moving onto second residue.\n";
    std::vector<cds::Atom*> to_this_residue2_cycle_points = cdsSelections::FindCyclePoints(to_this_connection_atom2, this->GetToThisResidue2());
    // Need to reverse one of these, so when concatenated, they are ordered ok. This might not be ok.
    //std::reverse(to_this_residue2_cycle_points.begin(), to_this_residue2_cycle_points.end());
    std::reverse(from_this_residue1_cycle_points.begin(), from_this_residue1_cycle_points.end());
    // Now concatenate:
    from_this_residue1_cycle_points.insert( from_this_residue1_cycle_points.end(), to_this_residue2_cycle_points.begin(), to_this_residue2_cycle_points.end() );
    // Now that have a list of rotation points. Split into pairs and find rotatable bonds between them
    bool found = false;
    std::vector<cds::Atom*> connecting_atoms = {from_this_connection_atom1, to_this_connection_atom2};
//     std::cout << "cycle point atoms are:\n";
//     for(auto & atom : from_this_residue1_cycle_points)
//     {
//         std::cout << atom->getId() << "\n";
//     }
//     std::cout << "\n";
    std::vector<RotatableDihedral> rotatableDihedralsInBranches;
    for(long unsigned int i = 0; i < from_this_residue1_cycle_points.size(); i = i+2)
    {
        cds::Atom* cycle_point1 = from_this_residue1_cycle_points.at(i);
        cds::Atom* cycle_point2 = from_this_residue1_cycle_points.at(i+1);

        found = false;
        connecting_atoms.clear();
     //   std::cout << "Finding Path between:" << cycle_point1->getId() << " and " << cycle_point2->getId() << "\n";
        cdsSelections::FindPathBetweenTwoAtoms(cycle_point1, this->GetFromThisResidue1(), cycle_point2, this->GetToThisResidue2(), &connecting_atoms, &found);
        cdsSelections::ClearAtomLabels(this->GetFromThisResidue1()); //ToDo change to free function or memeber function that clears labels.
        cdsSelections::ClearAtomLabels(this->GetToThisResidue2());
        // Find neighboring atoms needed to define dihedral. Pass in connecting atoms so don't find any of those.
        cds::Atom* neighbor1 =  cdsSelections::FindCyclePointNeighbor(connecting_atoms, cycle_point1, this->GetFromThisResidue1());
        cds::Atom* neighbor2 =  cdsSelections::FindCyclePointNeighbor(connecting_atoms, cycle_point2, this->GetToThisResidue2());
        // Insert these neighbors into list of connecting atoms, at beginning and end of vector.
        // connecting_atoms gets populated as it falls out, so list is reversed from what you'd expect
        std::reverse(connecting_atoms.begin(), connecting_atoms.end());
        connecting_atoms.insert(connecting_atoms.begin(), neighbor1);
        connecting_atoms.push_back(neighbor2);

        // std::cout << "Updated Path between:\n " << cycle_point1->getId() << " and " << cycle_point2->getId() << "\n";
//         for (const auto& atom : connecting_atoms)
//             std::cout << atom->getId() << "\n";
//         std::cout << "\n";
        cdsSelections::ClearAtomLabels(this->GetFromThisResidue1()); //ToDo change to free function or member function that clears labels.
        cdsSelections::ClearAtomLabels(this->GetToThisResidue2());
        // This mess was made to address the branching in 2-7 and 2-8 linkages.
        // These branches are long enough that they need default torsions set.
        if (connecting_atoms.size() > 4) // Otherwise there are no torsions
        { // Only viable linkages. Throw if not >4?
            for(typename std::vector<cds::Atom*>::iterator it = connecting_atoms.begin()+1; it != connecting_atoms.end()-1; ++it)
            { //
                cds::Atom* connectionAtom = (*it);
                if ( (connectionAtom != cycle_point1) && (connectionAtom != cycle_point2) )
                {
                    for(auto &neighbor : connectionAtom->getNeighbors()) // Need an interator
                    {
                        if (std::find(connecting_atoms.begin(), connecting_atoms.end(), neighbor) == connecting_atoms.end()) // if not in the vector
                        {
                            cdsSelections::Branch branch(connectionAtom);
                            cdsSelections::FindEndsOfBranchesFromLinkageAtom(neighbor, connectionAtom, this->GetToThisResidue2(), &branch);
                            if (branch.IsBranchFound())
                            {
                                found = false;
                                std::vector<cds::Atom*> foundPath;
                                // This fills in foundPath:
                                cdsSelections::FindPathBetweenTwoAtoms(branch.GetRoot(), this->GetFromThisResidue1(), branch.GetEnd(), this->GetToThisResidue2(), &foundPath, &found);
                                cds::Atom* neighbor = cdsSelections::FindCyclePointNeighbor(foundPath, branch.GetRoot(), this->GetToThisResidue2());
                                //foundPath.insert(foundPath.begin(), neighbor);
                                foundPath.push_back(neighbor);
                               //  std::cout << "Found atoms:\n";
                           //      for (auto &atom: foundPath)
                             //        std::cout << atom->getId() << "\n";
                                std::vector<RotatableDihedral> temp = this->SplitAtomVectorIntoRotatableDihedrals(foundPath);
                                rotatableDihedralsInBranches.insert( rotatableDihedralsInBranches.end(), temp.begin(), temp.end() );
                            }
                        }
                    }
                }
            }
        } // End dealing with branching linkages
//         std::cout << "These are the assigned branched RotatableDihedrals:\n";
//         for (auto &dihedral : rotatableDihedralsInBranches)
//             std::cout << dihedral.Print();
    }
    std::vector<RotatableDihedral> RotatableDihedrals = this->SplitAtomVectorIntoRotatableDihedrals(connecting_atoms);
    // Add any linkage branches (in 2-7 and 2-8) to the rest.
    RotatableDihedrals.insert( RotatableDihedrals.end(), rotatableDihedralsInBranches.begin(), rotatableDihedralsInBranches.end() );
//     std::cout << "These are the assigned RotatableDihedrals:\n";
//     for (auto &dihedral : RotatableDihedrals)
//         std::cout << dihedral.Print();
    return RotatableDihedrals;
}

std::vector<RotatableDihedral> ResidueLinkage::SplitAtomVectorIntoRotatableDihedrals(std::vector<cds::Atom*> atoms)
{
    //Ok looking for sets of four atoms, but shifting along vector by one atom for each dihedral.
    // So four atoms will make one rotatable bond, five will make two bonds, six will make three etc.
    std::vector<RotatableDihedral> RotatableDihedrals_generated;
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
//        std::cout << "Creating rotatable dihedrals between: " << this->GetFromThisResidue1()->getName() << " and " << this->GetToThisResidue2()->getName() << "\n";
        for(typename std::vector<cds::Atom*>::iterator it1 = atoms.begin(); it1 != (atoms.end()-3); ++it1)
        {
            cds::Atom* atom1 = *it1;
            cds::Atom* atom2 = *(it1+1);
            cds::Atom* atom3 = *(it1+2);
            cds::Atom* atom4 = *(it1+3);
            RotatableDihedrals_generated.emplace_back(atom1, atom2, atom3, atom4, this->GetIfReversedAtomsThatMove());
        }
    }
    return RotatableDihedrals_generated;
}

gmml::MolecularMetadata::GLYCAM::DihedralAngleDataVector ResidueLinkage::FindMetadata() const
{
    gmml::MolecularMetadata::GLYCAM::DihedralAngleDataContainer DihedralAngleMetadata;
    gmml::MolecularMetadata::GLYCAM::DihedralAngleDataVector matching_entries = DihedralAngleMetadata.GetEntriesForLinkage(this->GetFromThisConnectionAtom1()->getName(),
            this->GetFromThisResidue1()->getName(),
            this->GetToThisConnectionAtom2()->getName(),
            this->GetToThisResidue2()->getName() );
    //std::cout << "Found these " << matching_entries.size() << " entries:\n";
    // for (const auto& entry : matching_entries)
    // {
    //     std::cout << entry.index_ << " : " << entry.atom1_ << ", " << entry.atom2_ << ", " << entry.atom3_ << ", " << entry.atom4_ << ", " << entry.default_angle_value_ << "\n";
    // }
    if (matching_entries.empty())
    {
        matching_entries = DihedralAngleMetadata.GetEntriesForLinkage(this->GetToThisConnectionAtom2()->getName(),
                this->GetToThisResidue2()->getName(),
                this->GetFromThisConnectionAtom1()->getName(),
                this->GetFromThisResidue1()->getName() );
        std::reverse(matching_entries.begin(), matching_entries.end()); // I think this will work...
    }
    if (matching_entries.empty())
    {
        std::stringstream ss;
        ss << "No Metadata entries found for connection between " << this->GetFromThisResidue1()->getName() << "@" << this->GetFromThisConnectionAtom1()->getName() << " and " << this->GetToThisResidue2()->getName() << "@" << GetToThisConnectionAtom2()->getName() << "\n";
        ss << "Note that order should be reducing atom - anomeric atom, but I've tried reversing the order and it didn't fix the issue.\n";
        gmml::log(__LINE__,__FILE__,gmml::ERR,ss.str());
        throw std::runtime_error(ss.str());
    }
    return matching_entries;
}

void ResidueLinkage::AddMetadataToRotatableDihedrals(gmml::MolecularMetadata::GLYCAM::DihedralAngleDataVector metadata)
{
    // Adding another sanity check to this insanity
    if (rotatableDihedrals_.size() > metadata.size())
    {
        std::stringstream ss;
        ss << "Problem with the metadata found in gmml for this linkage. The number of identified rotatable dihedrals: " << rotatableDihedrals_.size() << " is greater than the number of metadata items: " << metadata.size() << " found for this linkage:\n" << this->Print();
        gmml::log(__LINE__,__FILE__,gmml::WAR, ss.str());
        throw std::runtime_error(ss.str());
    }
    for (auto & rotatableDihedral : rotatableDihedrals_)
    {
        rotatableDihedral.ClearMetadata(); // First clear any metadata already in place.
    }
    for (const auto& entry : metadata)
    { // ToDo: This is insane, make sane.
        unsigned long int vector_position = (entry.number_of_bonds_from_anomeric_carbon_ - 1); // vectors start at 0.
//        std::cout << "Adding to position: "<< vector_position << " in vector of size: " << RotatableDihedrals_.size() << std::endl;
        if (vector_position < rotatableDihedrals_.size())
        {
            RotatableDihedral& currentRotatableDihedral = rotatableDihedrals_.at(vector_position);
            currentRotatableDihedral.AddMetadata(entry);
            if(entry.dihedral_angle_name_ == "Psi" && entry.atom4_.at(0) == 'H' )
            { // If it's a psi angle and is supposed to be defined by a H...
                if(!currentRotatableDihedral.IsThereHydrogenForPsiAngle())
                {
                    this->GetToThisResidue2()->addAtom(currentRotatableDihedral.CreateHydrogenAtomForPsiAngle());
                }
            }
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
            std::string message = "Tried to add metadata to a rotatable bond that does not exist.\nCheck both dihedralangledata metadata and ResidueLinkage::FindRotatableDihedralsConnectingResidues.\nNote this is normal for a sialic acid with multiple 2-7, 2-8 and or 2-9 linkages and this warning can be ignored\n";
            gmml::log(__LINE__,__FILE__,gmml::WAR, message);
            //throw std::runtime_error(message);
        }
    }
//    std::cout << "\n";
    return;
}

void ResidueLinkage::SetResidues(cds::Residue* residue1, cds::Residue* residue2)
{
    from_this_residue1_ = residue1;
    to_this_residue2_ = residue2;
}

void ResidueLinkage::SetConnectionAtoms(cds::Residue* residue1, cds::Residue* residue2)
{
    std::vector<cds::Atom*> connecting_atoms;
    bool found = false;
    cdsSelections::FindAtomsConnectingResidues(residue1->getAtoms().at(0), residue1, residue2, &connecting_atoms, &found);
    if (connecting_atoms.size() < 2)
    {
        throw std::runtime_error("Two residues passed into ResidueLinkage that have no connection atoms.");
    }
    from_this_connection_atom1_ = connecting_atoms.at(0);
    to_this_connection_atom2_ = connecting_atoms.at(1);
}

void ResidueLinkage::SetConformerUsingMetadata(bool useRanges, int conformerNumber)
{
    for (auto &entry : rotatableDihedrals_)
    {
        entry.SetSpecificAngleEntryUsingMetadata(useRanges, conformerNumber);
    }
}


unsigned long long  ResidueLinkage::GenerateIndex()
{
    static unsigned long long s_ResidueLinkageIndex = 0; // static keyword means it is created only once and persists beyond scope of code block.
    return s_ResidueLinkageIndex++; // makes copy of s_AtomIndex, increments the real s_AtomIndex, then returns the value in the copy
}
