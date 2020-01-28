#include <random>
#include "../../../includes/GeometryTopology/ResidueLinkages/residue_linkage.hpp"
#include "../../../includes/MolecularModeling/overlaps.hpp"
#include "../../../includes/External_Libraries/PCG/pcg_extras.hpp"
#include "../../../includes/External_Libraries/PCG/pcg_random.hpp"
#include "../../../includes/MolecularMetadata/GLYCAM/glycam06residuecodes.hpp" // For lookup in GetName function

// Seed with a real random value, if available
static pcg_extras::seed_seq_from<std::random_device> seed_source;
// Make a random number engine
static pcg32 rng(seed_source);

//////////////////////////////////////////////////////////
//                    TYPE DEFINITION                   //
//////////////////////////////////////////////////////////

typedef std::vector<Rotatable_dihedral> RotatableDihedralVector;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////

Residue_linkage::Residue_linkage() {} // Do nothin

Residue_linkage::Residue_linkage(Residue *nonReducingResidue1, Residue *reducingResidue2, bool reverseAtomsThatMove)
{
    isExtraAtoms_ = false;
    this->InitializeClass(nonReducingResidue1, reducingResidue2, reverseAtomsThatMove);
}

Residue_linkage::Residue_linkage(Residue *nonReducingResidue1, Residue *reducingResidue2, AtomVector alsoMovingAtoms, bool reverseAtomsThatMove)
{ // Order of calling functions is important!
    this->AddExtraAtomsThatMove(alsoMovingAtoms);
    this->InitializeClass(nonReducingResidue1, reducingResidue2, reverseAtomsThatMove);
}
//////////////////////////////////////////////////////////
//                       ACCESSOR                       //
//////////////////////////////////////////////////////////

ResidueVector Residue_linkage::GetResidues()
{
    ResidueVector residues {from_this_residue1_, to_this_residue2_};
    return residues;
}

bool Residue_linkage::GetIfReversedAtomsThatMove()
{
    return reverseAtomsThatMove_;
}


RotatableDihedralVector Residue_linkage::GetRotatableDihedralsWithMultipleRotamers()
{
    RotatableDihedralVector returningDihedrals;
    for (auto &entry : rotatable_dihedrals_)
    {
        if (entry.GetNumberOfRotamers() > 1)
        {
            returningDihedrals.push_back(entry);
        }
    }
    return returningDihedrals;
}

RotatableDihedralVector Residue_linkage::GetRotatableDihedrals() const
{
    if (rotatable_dihedrals_.empty())
    {
        // LOG IT
        //std::cerr << "rotatable_dihedrals in this linkage is empty: " << from_this_residue1_->GetId() << "-" << to_this_residue2_->GetId() << std::endl;
    }
    return rotatable_dihedrals_;
}

int Residue_linkage::GetNumberOfShapes() // Can have conformers (sets of rotamers) or permutations of rotamers
{
    int numberOfShapes = 1;
    if (rotatable_dihedrals_.at(0).GetMetadata().at(0).rotamer_type_.compare("permutation")==0)
    {
        for (auto &entry : rotatable_dihedrals_)
        {
            numberOfShapes = (numberOfShapes * entry.GetNumberOfRotamers());
        }
    }
    else if (rotatable_dihedrals_.at(0).GetMetadata().at(0).rotamer_type_.compare("conformer")==0)
    {
        numberOfShapes = rotatable_dihedrals_.size();
    }
    return numberOfShapes;
}

Residue* Residue_linkage::GetFromThisResidue1()
{
    return from_this_residue1_;
}

Residue* Residue_linkage::GetToThisResidue2()
{
    return to_this_residue2_;
}

Atom* Residue_linkage::GetFromThisConnectionAtom1()
{
    return from_this_connection_atom1_;
}

Atom* Residue_linkage::GetToThisConnectionAtom2()
{
    return to_this_connection_atom2_;
}

bool Residue_linkage::CheckIfConformer()
{
    if (rotatable_dihedrals_.empty())
        std::cerr << "Error in Residue_linkage::checkIfConformer as rotatable_dihedrals_.empty()\n";
    else if (rotatable_dihedrals_.at(0).GetMetadata().empty())
        std::cerr << "Error in Residue_linkage::checkIfConformer as rotatable_dihedrals_.at(0).GetMetadata().empty()\n";
    else
    {
        if (rotatable_dihedrals_.at(0).GetMetadata().at(0).rotamer_type_.compare("permutation")==0)
            return false;
        else
            return true;
    }
    return false; //Default to shut up the compiler
}

bool Residue_linkage::GetIfExtraAtoms()
{
    return isExtraAtoms_;
}

AtomVector Residue_linkage::GetExtraAtoms()
{
    return extraAtomsThatMove_;
}

unsigned long long Residue_linkage::GetIndex()
{
    return index_;
}

std::string Residue_linkage::GetName()
{
    if (!name_.empty())
    {
        return name_;
    }
    return this->DetermineLinkageNameFromResidueNames();
}

std::string Residue_linkage::DetermineLinkageNameFromResidueNames()
{
    gmml::MolecularMetadata::GLYCAM::Glycam06ResidueNamesToCodesLookupContainer nameLookup;
    std::string residue1Name = nameLookup.GetResidueForCode(this->GetFromThisResidue1()->GetName());
    std::string residue2Name = nameLookup.GetResidueForCode(this->GetToThisResidue2()->GetName());
   // std::cout << this->GetFromThisConnectionAtom1()->GetName() << std::endl;
   // std::cout << this->GetToThisConnectionAtom2()->GetName() << std::endl;
    std::string atom1Name = this->GetFromThisConnectionAtom1()->GetName();
    std::string atom2Name = this->GetToThisConnectionAtom2()->GetName();
    char link1 = *atom1Name.rbegin(); //
    char link2 = *atom2Name.rbegin(); // Messy for Acetyl.
    std::stringstream linkageName;
    linkageName << residue1Name << link1 << "-" << link2 << residue2Name;
    this->SetName(linkageName.str());
    return linkageName.str();
}

int Residue_linkage::GetNumberOfRotatableDihedrals()
{
   return rotatable_dihedrals_.size();
}

//////////////////////////////////////////////////////////
//                       MUTATOR                        //
//////////////////////////////////////////////////////////

void Residue_linkage::SetRotatableDihedrals(RotatableDihedralVector rotatableDihedrals)
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
    isExtraAtoms_ = true;
}


//////////////////////////////////////////////////////////
//                       FUNCTIONS                      //
//////////////////////////////////////////////////////////

void Residue_linkage::GenerateAllShapesUsingMetadata()
{

}

void Residue_linkage::SetDefaultShapeUsingMetadata()
{
    for (auto &entry : rotatable_dihedrals_)
    {
        entry.SetSpecificAngleEntryUsingMetadata(); // Default is first entry
    }
}

void Residue_linkage::SetRandomShapeUsingMetadata(bool useRanges)
{
    if (rotatable_dihedrals_.at(0).GetMetadata().at(0).rotamer_type_.compare("permutation")==0)
    {
        for (auto &entry : rotatable_dihedrals_)
        {
            entry.SetRandomAngleEntryUsingMetadata(useRanges);
        }
    }
    else if (rotatable_dihedrals_.at(0).GetMetadata().at(0).rotamer_type_.compare("conformer")==0)
    {
        int numberOfConformers = rotatable_dihedrals_.at(0).GetMetadata().size();
        std::uniform_int_distribution<> distr(0, (numberOfConformers - 1)); // define the range
        int randomlySelectedConformerNumber = distr(rng);
        for (auto &entry : rotatable_dihedrals_)
        {
            entry.SetSpecificAngleEntryUsingMetadata(useRanges, randomlySelectedConformerNumber);
        }
    }
}

// This is dangerous if used by non-conformer linkage. Change it to check.
void Residue_linkage::SetSpecificShapeUsingMetadata(int shapeNumber, bool useRanges)
{
    for (auto &entry : rotatable_dihedrals_)
    {
        entry.SetSpecificAngleEntryUsingMetadata(useRanges, shapeNumber);
    }
}

void Residue_linkage::SetCustomDihedralAngles(std::vector <double> dihedral_angles)
{
    if(dihedral_angles.size() == rotatable_dihedrals_.size())
    {
        std::vector <double>::iterator dihedral_angle_iterator = dihedral_angles.begin();
        for (auto &rotatable_dihedral : rotatable_dihedrals_)
        {
            rotatable_dihedral.SetDihedralAngle(*dihedral_angle_iterator);
            ++dihedral_angle_iterator;
        }
    }
    else
    {   // Really need to figure out this throwing exceptions lark.
        std::cerr << "ERROR; attempted to set dihedral angles for set of dihedrals but with mismatching number of bonds to angles\n" << std::endl;
    }
}

void Residue_linkage::SetShapeToPrevious()
{
    for(RotatableDihedralVector::iterator rotatable_dihedral = rotatable_dihedrals_.begin(); rotatable_dihedral != rotatable_dihedrals_.end(); ++rotatable_dihedral)
    {
        rotatable_dihedral->SetDihedralAngleToPrevious();
    }
}

void Residue_linkage::SetRandomDihedralAngles()
{
    for(RotatableDihedralVector::iterator rotatable_dihedral = rotatable_dihedrals_.begin(); rotatable_dihedral != rotatable_dihedrals_.end(); ++rotatable_dihedral)
    {
        rotatable_dihedral->RandomizeDihedralAngle();
    }
}

void Residue_linkage::DetermineAtomsThatMove()
{
    for(RotatableDihedralVector::iterator rotatable_dihedral = rotatable_dihedrals_.begin(); rotatable_dihedral != rotatable_dihedrals_.end(); ++rotatable_dihedral)
    {
        rotatable_dihedral->DetermineAtomsThatMove();
    }
}

void Residue_linkage::SimpleWiggle(AtomVector overlapAtomSet1, AtomVector overlapAtomSet2, double overlapTolerance, int angleIncrement)
{
    double current_overlap = gmml::CalculateAtomicOverlapsBetweenNonBondedAtoms(overlapAtomSet1, overlapAtomSet2);
    double lowest_overlap = current_overlap;
    RotatableDihedralVector rotatable_bond_vector = this->GetRotatableDihedrals();
    for(auto &rotatable_dihedral : rotatable_bond_vector)
    {
        double best_dihedral_angle = rotatable_dihedral.CalculateDihedralAngle();
  //      std::cout << "Starting new linkage with best angle as " << best_dihedral_angle << "\n";
        gmml::MolecularMetadata::GLYCAM::DihedralAngleDataVector metadata_entries = rotatable_dihedral.GetMetadata();
        for(auto &metadata : metadata_entries)
        {
            double lower_bound = (metadata.default_angle_value_ - metadata.lower_deviation_);
            double upper_bound = (metadata.default_angle_value_ + metadata.upper_deviation_);
            double current_dihedral = lower_bound;
            while(current_dihedral <= upper_bound )
            {
                rotatable_dihedral.SetDihedralAngle(current_dihedral);
                current_overlap = gmml::CalculateAtomicOverlapsBetweenNonBondedAtoms(overlapAtomSet1, overlapAtomSet2);
           //     std::cout << "Dihedral(best): " << current_dihedral << "(" << best_dihedral_angle << ")" <<  ". Overlap(best): " << current_overlap << "(" << lowest_overlap << ")" << "\n";
                if (lowest_overlap >= (current_overlap + 0.01)) // 0.01 otherwise rounding errors
                {
                   // rotatable_dihedral.Print();
                    lowest_overlap = current_overlap;
                    best_dihedral_angle = current_dihedral;
           //         std::cout << "Best angle is now " << best_dihedral_angle << "\n";
                }
                // Perfer angles closer to default.
                else if ( (lowest_overlap == current_overlap) && (std::abs(metadata.default_angle_value_ - best_dihedral_angle )
                                                                  > std::abs(metadata.default_angle_value_ - current_dihedral)) )
                {
                    best_dihedral_angle = current_dihedral;
                }
                current_dihedral += angleIncrement; // increment
            }
        }
   //     std::cout << "Setting best angle as " << best_dihedral_angle << "\n";
        rotatable_dihedral.SetDihedralAngle(best_dihedral_angle);
        if(lowest_overlap <= overlapTolerance)
            return;
    }
    return; // Note possibility of earlier return above
}

void Residue_linkage::SetIndex(unsigned long long index)
{
    index_ = index;
}

void Residue_linkage::SetName(std::string name)
{
    name_ = name;
}

//////////////////////////////////////////////////////////
//                       DISPLAY FUNCTION               //
//////////////////////////////////////////////////////////

void Residue_linkage::Print()
{
    for(RotatableDihedralVector::iterator rotatable_dihedral = rotatable_dihedrals_.begin(); rotatable_dihedral != rotatable_dihedrals_.end(); ++rotatable_dihedral)
    {
        rotatable_dihedral->Print();
    }
}

//////////////////////////////////////////////////////////
//                       OPERATORS                      //
//////////////////////////////////////////////////////////

std::ostream& operator<<(std::ostream& os, const Residue_linkage& residue_linkage)
{
    RotatableDihedralVector rotatable_dihedrals = residue_linkage.GetRotatableDihedrals();
    for(RotatableDihedralVector::iterator rotatable_dihedral = rotatable_dihedrals.begin(); rotatable_dihedral != rotatable_dihedrals.end(); ++rotatable_dihedral)
    {
        os << (*rotatable_dihedral);
    }
    return os;
} // operator<<


//////////////////////////////////////////////////////////
//                    PRIVATE FUNCTIONS                 //
//////////////////////////////////////////////////////////

void Residue_linkage::InitializeClass(Residue *from_this_residue1, Residue *to_this_residue2, bool reverseAtomsThatMove)
{
    this->SetResidues(from_this_residue1, to_this_residue2);
    this->SetIfReversedAtomsThatMove(reverseAtomsThatMove);
    this->SetConnectionAtoms(from_this_residue1_, to_this_residue2_);
    if(this->CheckIfViableLinkage())
    {
       // std::cout << "Finding connection between " << from_this_residue1->GetId() << " :: " << to_this_residue2->GetId() << std::endl;
        rotatable_dihedrals_ = this->FindRotatableDihedralsConnectingResidues(from_this_connection_atom1_, to_this_connection_atom2_);
        gmml::MolecularMetadata::GLYCAM::DihedralAngleDataVector metadata = this->FindMetadata(from_this_connection_atom1_, to_this_connection_atom2_);
        this->AddMetadataToRotatableDihedrals(metadata);
    }
    this->SetIndex(this->GenerateIndex());
}

bool Residue_linkage::CheckIfViableLinkage()
{
    bool viable = true; // Assume ok.
    for(auto &residue : this->GetResidues())
    {
        int heavyAtomCount = 0;
        for(auto &atom : residue->GetAtoms())
            if (atom->GetName().at(0) != 'H')
                ++heavyAtomCount;
        if (heavyAtomCount <= 2)
            viable = false; // if either residue has two few atoms, it's not a viable linkage
    }
    return viable;
}

RotatableDihedralVector Residue_linkage::FindRotatableDihedralsConnectingResidues(Atom *from_this_connection_atom1, Atom *to_this_connection_atom2)
{
    // Going to ignore tags etc.
    // Given two residues that are connected. Find connecting atoms.
    // Search neighbors other than connected atom. Ie search out in both directions, but remain within same residue.
    // Warning, residue may have fused cycles!
    // Will fail for non-protein residues without cycles. As don't have a non-rotatable bond to anchor from. Can code that later (and deal with branches from these residues).
   // std::cout << "Finding rot bonds for " << from_this_connection_atom1->GetResidue()->GetId() << " and " << to_this_connection_atom2->GetResidue()->GetId() << "\n";

    AtomVector from_this_residue1_cycle_points = selection::FindCyclePoints(from_this_connection_atom1);
  //  std::cout << "Moving onto second residue.\n";
    AtomVector to_this_residue2_cycle_points = selection::FindCyclePoints(to_this_connection_atom2);
    // Need to reverse one of these, so when concatenated, they are ordered ok. This might not be ok.
//    std::reverse(to_this_residue2_cycle_points.begin(), to_this_residue2_cycle_points.end());
    std::reverse(from_this_residue1_cycle_points.begin(), from_this_residue1_cycle_points.end());
    // Now concatenate:
    from_this_residue1_cycle_points.insert( from_this_residue1_cycle_points.end(), to_this_residue2_cycle_points.begin(), to_this_residue2_cycle_points.end() );
    // Now that have a list of rotation points. Split into pairs and find rotatable bonds between them
    bool found = false;
    AtomVector connecting_atoms = {from_this_connection_atom1, to_this_connection_atom2};
  //  std::cout << "cycle point atoms are:\n";
  //  for(auto & atom : from_this_residue1_cycle_points)
   //     std::cout << atom->GetId();
 //   std::cout << "\n";

    for(int i = 0; i < from_this_residue1_cycle_points.size(); i = i+2)
    {
   //     std::cout << "Oh ya, this seems like a great place to crash right now\n";
        Atom *cycle_point1 = from_this_residue1_cycle_points.at(i);
        Atom *cycle_point2 = from_this_residue1_cycle_points.at(i+1);

        found = false;
        connecting_atoms.clear();
  //      std::cout << "Finding Path between:" << cycle_point1->GetId() << cycle_point2->GetId();
        selection::FindPathBetweenTwoAtoms(cycle_point1, cycle_point2, &connecting_atoms, &found);
        selection::ClearAtomDescriptions(cycle_point1->GetResidue());
        selection::ClearAtomDescriptions(cycle_point2->GetResidue());
        // Find neighboring atoms needed to define dihedral. Pass in connecting atoms so don't find any of those.
        Atom *neighbor1 =  selection::FindCyclePointNeighbor(connecting_atoms, cycle_point1);
        Atom *neighbor2 =  selection::FindCyclePointNeighbor(connecting_atoms, cycle_point2);
        // Insert these neighbors into list of connecting atoms, at beginning and end of vector.
        // connecting_atoms gets populated as it falls out, so list is reveresed from what you'd expect
        std::reverse(connecting_atoms.begin(), connecting_atoms.end());
        connecting_atoms.insert(connecting_atoms.begin(), neighbor1);
        connecting_atoms.push_back(neighbor2);

       // std::cout << "Updated Path between:\n " << cycle_point1->GetId() << cycle_point2->GetId();
     //   for (const auto& atom : connecting_atoms)
         //   std::cout << atom->GetId();
    }
    RotatableDihedralVector rotatable_dihedrals = this->SplitAtomVectorIntoRotatableDihedrals(connecting_atoms);
    return rotatable_dihedrals;
}

RotatableDihedralVector Residue_linkage::SplitAtomVectorIntoRotatableDihedrals(AtomVector atoms)
{
    //Ok looking for sets of four atoms, but shifting along vector by one atom for each dihedral.
    // So four atoms will make one rotatable bond, five will make two bonds, six will make three etc.
    RotatableDihedralVector rotatable_dihedrals_generated;
    if(atoms.size() < 4)
    {
//        std::cout << "ERROR; in Residue_linkage::SplitAtomVectorIntoRotatableDihedrals, not enough atoms in atom vector: " << atoms.size() << std::endl;
    }
    else
    {
        for(AtomVector::iterator it1 = atoms.begin(); it1 != (atoms.end()-3); ++it1)
        {
            Atom *atom1 = *it1;
            Atom *atom2 = *(it1+1);
            Atom *atom3 = *(it1+2);
            Atom *atom4 = *(it1+3);
            rotatable_dihedrals_generated.emplace_back(atom1, atom2, atom3, atom4, this->GetIfReversedAtomsThatMove());
        }
    }
    return rotatable_dihedrals_generated;
}

gmml::MolecularMetadata::GLYCAM::DihedralAngleDataVector Residue_linkage::FindMetadata(Atom *from_this_connection_atom1, Atom *to_this_connection_atom2)
{
    gmml::MolecularMetadata::GLYCAM::DihedralAngleDataContainer DihedralAngleMetadata;
    gmml::MolecularMetadata::GLYCAM::DihedralAngleDataVector matching_entries = DihedralAngleMetadata.GetEntriesForLinkage(from_this_connection_atom1, to_this_connection_atom2);
//    std::cout << "Found these " << matching_entries.size() << " entries:\n";
//    for (const auto& entry : matching_entries)
//    {
//        std::cout << entry.index_ << " : " << entry.atom1_ << ", " << entry.atom2_ << ", " << entry.atom3_ << ", " << entry.atom4_ << ", " << entry.default_angle_value_ << "\n";
//    }
    if (matching_entries.empty())
    {
        std::cerr << "No Metadata entries found for connection between " << from_this_connection_atom1->GetId() << " and " << to_this_connection_atom2 << "\n";
        std::cerr << "Note that order should be reducing atom - anomeric atom\n";
    }
    return matching_entries;
}

void Residue_linkage::AddMetadataToRotatableDihedrals(gmml::MolecularMetadata::GLYCAM::DihedralAngleDataVector metadata)
{
    for (const auto& entry : metadata)
    {
//        int bond_number = int (entry.number_of_bonds_from_anomeric_carbon_); // typecast to an int
        int vector_position = (entry.number_of_bonds_from_anomeric_carbon_ - 1); // vectors start at 0.
//        std::cout << "Adding to position: "<< vector_position << " in vector of size: " << rotatable_dihedrals_.size() << std::endl;
        if (vector_position <= rotatable_dihedrals_.size())
        {
            rotatable_dihedrals_.at(vector_position).AddMetadata(entry);
//            std::cout << "Set " << entry.dihedral_angle_name_ << " meta data for "
//                      << rotatable_dihedrals_.at(vector_position).GetAtoms().at(0)->GetId() << ", "
//                      << rotatable_dihedrals_.at(vector_position).GetAtoms().at(1)->GetId() << ", "
//                      << rotatable_dihedrals_.at(vector_position).GetAtoms().at(2)->GetId() << ", "
//                      << rotatable_dihedrals_.at(vector_position).GetAtoms().at(3)->GetId() << "\n";
//            std::cout << "Added " << entry.index_ << " = " << entry.default_angle_value_ << " to: \n";
//            rotatable_dihedrals_.at(vector_position).Print();
        }
        else
        {
           std::cerr << "Huge problem in residue_linkage.cpp AddMetadataToRotatableDihedrals. Tried to add metadata to a rotatable bond that does not exist.\n"
                         "Check both dihedralangledata metadata and Residue_linkage::FindRotatableDihedralsConnectingResidues." << std::endl;
        }
    }
}

void Residue_linkage::SetResidues(Residue *residue1, Residue *residue2)
{
    from_this_residue1_ = residue1;
    to_this_residue2_ = residue2;
}

void Residue_linkage::SetConnectionAtoms(Residue *residue1, Residue *residue2)
{
    AtomVector connecting_atoms;
    bool found = false;
    selection::FindAtomsConnectingResidues(residue1->GetAtoms().at(0), residue2, &connecting_atoms, &found);
    from_this_connection_atom1_ = connecting_atoms.at(0);
    to_this_connection_atom2_ = connecting_atoms.at(1);
}

void Residue_linkage::SetConformerUsingMetadata(bool useRanges, int conformerNumber)
{
    for (auto &entry : rotatable_dihedrals_)
    {
        entry.SetSpecificAngleEntryUsingMetadata(useRanges, conformerNumber);
    }
}

unsigned long long Residue_linkage::GenerateIndex()
{
    static unsigned long long s_ResidueLinkageIndex = 0; // static keyword means it is created only once and persists beyond scope of code block.
    return s_ResidueLinkageIndex++; // makes copy of s_AtomIndex, increments the real s_AtomIndex, then returns the value in the copy
}
