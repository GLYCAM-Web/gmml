#include <random>
#include "../../../includes/GeometryTopology/ResidueLinkages/rotatable_dihedral.hpp"
#include "../../../includes/MolecularModeling/atomnode.hpp" // For UpdateAtomsIfPsi
#include "../../../includes/utils.hpp"
#include "../../../includes/External_Libraries/PCG/pcg_random.hpp"
#include "../../../includes/GeometryTopology/geometrytopology.hpp"

// Seed with a real random value, if available
static pcg_extras::seed_seq_from<std::random_device> seed_source;
// Make a random number engine
static pcg32 rng(seed_source);

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////

// Of the next three forms, I'll probably use only one and delete the others later
Rotatable_dihedral::Rotatable_dihedral(Atom *atom1, Atom *atom2, Atom *atom3, Atom *atom4, bool reverseAtomsThatMove)
{
    AtomVector atoms {atom1, atom2, atom3, atom4};
    this->Initialize(atoms, reverseAtomsThatMove);
}
Rotatable_dihedral::Rotatable_dihedral(Atom *atom1, Atom *atom2, Atom *atom3, Atom *atom4, AtomVector extraAtomsThatMove, bool reverseAtomsThatMove)
{
    AtomVector atoms {atom1, atom2, atom3, atom4};
    this->Initialize(atoms, reverseAtomsThatMove);
    this->AddExtraAtomsThatMove(extraAtomsThatMove);
}

//Rotatable_dihedral::Rotatable_dihedral(AtomVector atoms, bool reverseAtomsThatMove)
//{
//    this->Initialize(atoms, reverseAtomsThatMove);
//}

//Rotatable_dihedral::Rotatable_dihedral(AtomVector atoms, AtomVector atoms_that_move)
//{
//    this->SetAtoms(atoms);
//    this->SetAtomsThatMove(atoms_that_move);
//}

//////////////////////////////////////////////////////////
//                       ACCESSOR                       //
//////////////////////////////////////////////////////////

double Rotatable_dihedral::CalculateDihedralAngle() const
{
    GeometryTopology::Coordinate* a1 = atom1_->GetCoordinate();
    GeometryTopology::Coordinate* a2 = atom2_->GetCoordinate();
    GeometryTopology::Coordinate* a3 = atom3_->GetCoordinate();
    GeometryTopology::Coordinate* a4 = atom4_->GetCoordinate();

    GeometryTopology::Coordinate b1 = a2;
    b1.operator -(*a1);
    GeometryTopology::Coordinate b2 = a3;
    b2.operator -(*a2);
    GeometryTopology::Coordinate b3 = a4;
    b3.operator -(*a3);
    GeometryTopology::Coordinate b4 = b2;
    b4.operator *(-1);

    GeometryTopology::Coordinate b2xb3 = b2;
    b2xb3.CrossProduct(b3);

    GeometryTopology::Coordinate b1_m_b2n = b1;
    b1_m_b2n.operator *(b2.length());

    GeometryTopology::Coordinate b1xb2 = b1;
    b1xb2.CrossProduct(b2);

    double current_dihedral_angle = atan2(b1_m_b2n.DotProduct(b2xb3), b1xb2.DotProduct(b2xb3));

    return (current_dihedral_angle * (180 / gmml::PI_RADIAN) ); // Convert to degrees
}

AtomVector Rotatable_dihedral::GetAtoms() const
{
    AtomVector atoms = {atom1_, atom2_, atom3_, atom4_};
    return atoms;
}

AtomVector Rotatable_dihedral::GetAtomsThatMove()
{
    return atoms_that_move_;
}

bool Rotatable_dihedral::GetIsAtomsThatMoveReversed()
{
    return isAtomsThatMoveReversed_;
}

double Rotatable_dihedral::GetPreviousDihedralAngle()
{
    return previous_dihedral_angle_;
}

DihedralAngleDataVector Rotatable_dihedral::GetMetadata()
{
    return assigned_metadata_;
}

int Rotatable_dihedral::GetNumberOfRotamers()
{
    if (assigned_metadata_.empty())
    {
        std::cerr << "Error in Rotatable_dihedral::GetNumberOfRotamers; no metadata has been set.\n";
        return 0;
    }
    else
    {
        return assigned_metadata_.size();
    }
}

std::vector<double> Rotatable_dihedral::GetAllPossibleAngleValues(int interval)
{
    std::vector<double> allPossibleAngleValues;
    if (assigned_metadata_.empty())
    {
        std::cerr << "Error in Rotatable_dihedral::GetAllPossibleAngleValues; no metadata has been set.\n";
    }
    else
    {
        gmml::MolecularMetadata::GLYCAM::DihedralAngleDataVector metadata_entries = this->GetMetadata();
        //std::cout << metadata_entries.at(0).dihedral_angle_name_;
        for(auto &metadata : metadata_entries)
        {
            double lower_bound = (metadata.default_angle_value_ - metadata.lower_deviation_);
            double upper_bound = (metadata.default_angle_value_ + metadata.upper_deviation_);
            double current_dihedral = lower_bound;
            while(current_dihedral <= upper_bound)
            {
                allPossibleAngleValues.push_back(current_dihedral);
                current_dihedral += interval; // increment
            }
        }
    }
    return allPossibleAngleValues;
}

//////////////////////////////////////////////////////////
//                       MUTATOR                        //
//////////////////////////////////////////////////////////

void Rotatable_dihedral::DetermineAtomsThatMove()
{
    // In keeping with giving residues as GlcNAc1-4Gal, and wanting the moving atoms to be in the opposite direction (defaults):
    AtomVector atoms_that_move;
    if (this->GetIsAtomsThatMoveReversed())
    {
        atoms_that_move.push_back(atom3_);
        atom2_->FindConnectedAtoms(atoms_that_move);
    }
    else
    {
        atoms_that_move.push_back(atom2_);
        atom3_->FindConnectedAtoms(atoms_that_move);
    }
    this->SetAtomsThatMove(atoms_that_move);
}

void Rotatable_dihedral::AddExtraAtomsThatMove(AtomVector extraAtoms)
{
    extra_atoms_that_move_.insert(extra_atoms_that_move_.end(), extraAtoms.begin(), extraAtoms.end());
    return;
}

// Only this function uses radians. Everything else, in and out, should be degrees.
void Rotatable_dihedral::SetDihedralAngle(double dihedral_angle)
{
    if(!this->wasEverRotated_)
    {
        this->SetWasEverRotated(true);
        this->DetermineAtomsThatMove();
    }
    GeometryTopology::Coordinate* a1 = atom1_->GetCoordinate();
    GeometryTopology::Coordinate* a2 = atom2_->GetCoordinate();
    GeometryTopology::Coordinate* a3 = atom3_->GetCoordinate();
    GeometryTopology::Coordinate* a4 = atom4_->GetCoordinate();
    if(this->GetIsAtomsThatMoveReversed())
    { // A function that was agnostic of this would be great.
        a1 = atom4_->GetCoordinate();
        a2 = atom3_->GetCoordinate();
        a3 = atom2_->GetCoordinate();
        a4 = atom1_->GetCoordinate();
    }
    GeometryTopology::Coordinate b1 = a2;
    b1.operator -(*a1);
    GeometryTopology::Coordinate b2 = a3;
    b2.operator -(*a2);
    GeometryTopology::Coordinate b3 = a4;
    b3.operator -(*a3);
    GeometryTopology::Coordinate b4 = b2;
    b4.operator *(-1);

    GeometryTopology::Coordinate b2xb3 = b2;
    b2xb3.CrossProduct(b3);

    GeometryTopology::Coordinate b1_m_b2n = b1;
    b1_m_b2n.operator *(b2.length());

    GeometryTopology::Coordinate b1xb2 = b1;
    b1xb2.CrossProduct(b2);

    double current_dihedral = atan2(b1_m_b2n.DotProduct(b2xb3), b1xb2.DotProduct(b2xb3));
    double** dihedral_angle_matrix = gmml::GenerateRotationMatrix(&b4, a2, current_dihedral - gmml::ConvertDegree2Radian(dihedral_angle));
    this->RecordPreviousDihedralAngle(gmml::ConvertRadian2Degree(current_dihedral));

    // Yo you should add something here that checks if atoms_that_move_ is set. Yeah you.

//    std::cout << "Setting dihedral for " << atom1_->GetId() << ":"  << atom2_->GetId() << ":"  << atom3_->GetId() << ":"  << atom4_->GetId() <<  ": " << dihedral_angle << "\n";
    for(AtomVector::iterator it = atoms_that_move_.begin(); it != atoms_that_move_.end(); it++)
    {
        Atom *atom = *it;
    //    std::cout << ", " << atom->GetId();
        GeometryTopology::Coordinate* atom_coordinate = atom->GetCoordinate();
        GeometryTopology::Coordinate result;
        result.SetX(dihedral_angle_matrix[0][0] * atom_coordinate->GetX() + dihedral_angle_matrix[0][1] * atom_coordinate->GetY() +
                dihedral_angle_matrix[0][2] * atom_coordinate->GetZ() + dihedral_angle_matrix[0][3]);
        result.SetY(dihedral_angle_matrix[1][0] * atom_coordinate->GetX() + dihedral_angle_matrix[1][1] * atom_coordinate->GetY() +
                dihedral_angle_matrix[1][2] * atom_coordinate->GetZ() + dihedral_angle_matrix[1][3]);
        result.SetZ(dihedral_angle_matrix[2][0] * atom_coordinate->GetX() + dihedral_angle_matrix[2][1] * atom_coordinate->GetY() +
                dihedral_angle_matrix[2][2] * atom_coordinate->GetZ() + dihedral_angle_matrix[2][3]);

        atom->GetCoordinate()->SetX(result.GetX());
        atom->GetCoordinate()->SetY(result.GetY());
        atom->GetCoordinate()->SetZ(result.GetZ());
    }
    if(!extra_atoms_that_move_.empty())
    {
        for(AtomVector::iterator it = extra_atoms_that_move_.begin(); it != extra_atoms_that_move_.end(); it++)
        {
            Atom *atom = *it;
        //    std::cout << ", " << atom->GetId();
            GeometryTopology::Coordinate* atom_coordinate = atom->GetCoordinate();
            GeometryTopology::Coordinate result;
            result.SetX(dihedral_angle_matrix[0][0] * atom_coordinate->GetX() + dihedral_angle_matrix[0][1] * atom_coordinate->GetY() +
                    dihedral_angle_matrix[0][2] * atom_coordinate->GetZ() + dihedral_angle_matrix[0][3]);
            result.SetY(dihedral_angle_matrix[1][0] * atom_coordinate->GetX() + dihedral_angle_matrix[1][1] * atom_coordinate->GetY() +
                    dihedral_angle_matrix[1][2] * atom_coordinate->GetZ() + dihedral_angle_matrix[1][3]);
            result.SetZ(dihedral_angle_matrix[2][0] * atom_coordinate->GetX() + dihedral_angle_matrix[2][1] * atom_coordinate->GetY() +
                    dihedral_angle_matrix[2][2] * atom_coordinate->GetZ() + dihedral_angle_matrix[2][3]);

            atom->GetCoordinate()->SetX(result.GetX());
            atom->GetCoordinate()->SetY(result.GetY());
            atom->GetCoordinate()->SetZ(result.GetZ());
        }
    }
  //  std::cout << std::endl;
    return;
}


void Rotatable_dihedral::SetDihedralAngleToPrevious()
{
    this->SetDihedralAngle(this->GetPreviousDihedralAngle());
}


double Rotatable_dihedral::RandomizeDihedralAngle()
{
    return Rotatable_dihedral::RandomizeDihedralAngleWithinRange(0.0, 360.0);
    //return (rand() % 360) + 1 - 180; // Can get same one everytime for testing
}

double Rotatable_dihedral::RandomizeDihedralAngleWithinRange(double min, double max)
{
//    std::random_device rd1; // obtain a random number from hardware
//    std::mt19937 eng1(rd1()); // seed the generator
    std::uniform_real_distribution<> angle_distribution(min, max); // define the range

//    double random_angle = angle_distribution(eng1);
    double random_angle = angle_distribution(rng);
   // std::cout << "Random angle is: " << random_angle << "\n";

    /*******************************************/
    /*               IMPORTANT                 */
    /*******************************************/

    this->SetDihedralAngle(random_angle); // THIS IS IMPORTANT!!! THIS SHOULD BE SEPARATED?!?! The two other functions call this one. Seems fragile.

    /*******************************************/
    /*               IMPORTANT                 */
    /*******************************************/

    return random_angle;
    //return rand() % (max + 1 - min) + min; // Can get same one everytime for testing
}

// OG for the CB, I need this class to know which metadata index it's currently set to. So this won't work for that.
double Rotatable_dihedral::RandomizeDihedralAngleWithinRanges(std::vector<std::pair<double,double>> ranges)
{
    // For usage, can do ranges.emplace_back(min, max);
    // Pass in a vector of pairs of ranges.
    // First select one of those ranges.
    // Then create an angle within the selected range.

    // Rando stuff from slack overflow:
//    std::random_device rd; // obtain a random number from hardware
//    std::mt19937 eng(rd()); // seed the generator
    std::uniform_int_distribution<> distr(0, (ranges.size() - 1)); // define the range

    // Select one of the ranges
    int range_selection = distr(rng);
    //std::cout << "Randomly selected range number " << range_selection << "\n";

    // create an angle within the selected range
    return this->RandomizeDihedralAngleWithinRange(ranges.at(range_selection).first, ranges.at(range_selection).second);
}

void Rotatable_dihedral::SetMetadata(DihedralAngleDataVector metadataVector)
{
    assigned_metadata_ = metadataVector;
    this->UpdateAtomsIfPsi();
}

void Rotatable_dihedral::AddMetadata(DihedralAngleData metadata)
{
    assigned_metadata_.push_back(metadata);
    this->UpdateAtomsIfPsi();
}

void Rotatable_dihedral::SetRandomAngleEntryUsingMetadata(bool useRanges)
{
    if (assigned_metadata_.empty())
    {
        std::cerr << "Error in Rotatable_dihedral::SetRandomAngleEntryUsingMetadata; no metadata has been set for:.\n"
                  << atom1_->GetId() << atom2_->GetId() << atom3_->GetId() << atom4_->GetId();
    }
    else if(assigned_metadata_.size() == 1)
    {
        const auto& entry = assigned_metadata_.at(0);
        if (useRanges)
        {
            double lower = (entry.default_angle_value_ - entry.lower_deviation_) ;
            double upper = (entry.default_angle_value_ + entry.upper_deviation_) ;
            this->RandomizeDihedralAngleWithinRange(lower, upper);
        }
        else
        {
            this->SetDihedralAngle(entry.default_angle_value_);
        }
    } // I need different behaviour here than previously. If use_ranges is false, just take the first of multiple rotamers, don't randomize.
    else if(assigned_metadata_.size() >= 2) // Some dihedral angles have multiple rotamers, thus mulitple ranges to select from. e.g -60, 60, 180
    {
        // OG update for carb builder (CB). This class needs to know what metadata entry it's currently set to.
        if (useRanges)
        {
            // first randomly pick one of the meta data entries
            std::uniform_int_distribution<> distr(0, (assigned_metadata_.size() - 1)); // define the range
            const auto& randomEntry = assigned_metadata_.at(distr(rng));
            double lower = (randomEntry.default_angle_value_ - randomEntry.lower_deviation_) ;
            double upper = (randomEntry.default_angle_value_ + randomEntry.upper_deviation_) ;
            this->RandomizeDihedralAngleWithinRange(lower, upper);
        }
        else
        { // Just set it to the first entries default value
            const auto& entry = assigned_metadata_.at(0);
            this->RandomizeDihedralAngleWithinRange(entry.default_angle_value_, entry.default_angle_value_);
        }
    }
    return;
}

void Rotatable_dihedral::SetSpecificAngleEntryUsingMetadata(bool useRanges, int angleEntryNumber)
{
    if (assigned_metadata_.empty())
    {
        std::cerr << "Error in Rotatable_dihedral::SetSpecificAngleUsingMetadata; no metadata has been set.\n";
    }
    else if (assigned_metadata_.size() <= angleEntryNumber)
    {
         std::cerr << "Error in Rotatable_dihedral::SetSpecificAngleUsingMetadata; angleEntryNumber of " << angleEntryNumber << " is too large as metadata.size() is " << assigned_metadata_.size() << ".\n";
    }
    else
    {
        const auto& entry = assigned_metadata_.at(angleEntryNumber);
        if (useRanges)
        {
            double lower = (entry.default_angle_value_ - entry.lower_deviation_) ;
            double upper = (entry.default_angle_value_ + entry.upper_deviation_) ;
            this->RandomizeDihedralAngleWithinRange(lower, upper);
        }
        else
        {
            this->SetDihedralAngle(entry.default_angle_value_);
        }
    }
}

//////////////////////////////////////////////////////////
//                  PRIVATE FUNCTIONS                   //
//////////////////////////////////////////////////////////

void Rotatable_dihedral::Initialize(AtomVector atoms, bool reverseAtomsThatMove)
{
    this->SetWasEverRotated(false);
    this->SetAtoms(atoms);
    this->SetIsAtomsThatMoveReversed(reverseAtomsThatMove);
    // this->DetermineAtomsThatMove(); // Will be done the first time SetDihedralAngle is called.
}

void Rotatable_dihedral::SetAtoms(AtomVector atoms)
{
    atom1_ = atoms.at(0);
    atom2_ = atoms.at(1);
    atom3_ = atoms.at(2);
    atom4_ = atoms.at(3);
}

void Rotatable_dihedral::SetAtomsThatMove(AtomVector atoms)
{
    atoms_that_move_ = atoms;
//    std::cout << "Set the following to move:\n";
//    for (auto &moving_atom : atoms_that_move_)
//    {
//        std::cout << moving_atom->GetId() << ", ";
//    }
//    std::cout << "\n";
}

void Rotatable_dihedral::SetIsAtomsThatMoveReversed(bool isAtomsThatMoveReversed)
{
    isAtomsThatMoveReversed_ = isAtomsThatMoveReversed;
}

void Rotatable_dihedral::RecordPreviousDihedralAngle(double dihedral_angle)
{
    previous_dihedral_angle_ = dihedral_angle;
}

void Rotatable_dihedral::UpdateAtomsIfPsi()
{
    for(auto &entry : assigned_metadata_)
    {
        // If it's a psi angle and is supposed to be defined by a H...
        if ((entry.dihedral_angle_name_.compare("Psi")==0) && (entry.atom4_.at(0)=='H'))
        {// Find the neighbor of current atom3 which is a hydrogen, and set it to be the new atom4;
            bool foundHydrogen = false;
            for(auto &neighbor : atom3_->GetNode()->GetNodeNeighbors())
            {
                if(neighbor->GetName().at(0)=='H')
                {
//                    std::cout << "In ";
//                    this->Print();
//                    std::cout << "Replaced atom4_ with " << neighbor->GetId() << "\n";
                    atom4_ = neighbor;
                    foundHydrogen = true;
                }
            }
            if (!foundHydrogen)
            {
                atom4_ = Rotatable_dihedral::CreateHydrogenAtomForPsi(atom3_);
            }
        }
    }
    return;
}

Atom* Rotatable_dihedral::CreateHydrogenAtomForPsi(Atom *centralAtom)
{
    if(centralAtom->GetNode()->GetNodeNeighbors().size() != 3)
    {
        std::cerr << "Error in Rotatable_dihedral::CreateHydrogenAtomForPsi. centralAtom neighbors = " <<
                  centralAtom->GetNode()->GetNodeNeighbors().size() << " for " << centralAtom->GetId() << std::endl;
    }
    GeometryTopology::CoordinateVector threeNeighborCoords;
    for (auto &neighbor : centralAtom->GetNode()->GetNodeNeighbors())
    {
        threeNeighborCoords.push_back(neighbor->GetCoordinate());
    }
    GeometryTopology::Coordinate newCoord = GeometryTopology::CreateMissingCoordinateForTetrahedralAtom(centralAtom->GetCoordinate(), threeNeighborCoords);
    Atom *newAtom = new Atom(centralAtom->GetResidue(), "HHH", newCoord);
    centralAtom->GetResidue()->AddAtom(newAtom);
    centralAtom->GetNode()->AddNodeNeighbor(newAtom);
    // Have to create a new node for this new atom.
    MolecularModeling::AtomNode *leakyNode = new MolecularModeling::AtomNode();
    newAtom->SetNode(leakyNode); // A jest!
    newAtom->GetNode()->AddNodeNeighbor(centralAtom);
    return newAtom;
}

void Rotatable_dihedral::SetWasEverRotated(bool wasEverRotated)
{
    wasEverRotated_ = wasEverRotated;
}

bool Rotatable_dihedral::CheckIfEverRotated()
{
    return wasEverRotated_;
}


//////////////////////////////////////////////////////////
//                       DISPLAY FUNCTION               //
//////////////////////////////////////////////////////////

void Rotatable_dihedral::Print()
{
//    std::cout << atom1_->GetName() << ", " << atom2_->GetName() << ", " << atom3_->GetName() << ", " << atom4_->GetName() << ": " << this->CalculateDihedralAngle()  << ".\n";
////    for(AtomVector::iterator it1 = atoms_that_move_.begin(); it1 != atoms_that_move_.end(); ++it1)
////    {
////        Atom *atom = *it1;
////        std::cout << atom->GetName() << ", ";
////    }
//    std::cout << std::endl;
}

//////////////////////////////////////////////////////////
//                       OPERATORS                      //
//////////////////////////////////////////////////////////

std::ostream& operator<<(std::ostream& os, const Rotatable_dihedral& rotatable_dihedral)
{
    AtomVector atoms = rotatable_dihedral.GetAtoms();
    os << atoms.at(0)->GetName() << ", " << atoms.at(1)->GetName() << ", " << atoms.at(2)->GetName() << ", " << atoms.at(3)->GetName() << ": " << rotatable_dihedral.CalculateDihedralAngle() << ".\n";
    return os;
} // operator<<
