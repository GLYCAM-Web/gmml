#include <random>
#include "../../../includes/GeometryTopology/ResidueLinkages/rotatable_dihedral.h"
#include "../../../includes/MolecularModeling/atomnode.hpp" // For UpdateAtomsIfPsi
#include "../../../includes/utils.hpp"
#include "../../../includes/External_Libraries/PCG/pcg_random.hpp"

// Seed with a real random value, if available
static pcg_extras::seed_seq_from<std::random_device> seed_source;
// Make a random number engine
static pcg32 rng(seed_source);

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////

// Of the next three forms, I'll probably use only one and delete the others later
Rotatable_dihedral::Rotatable_dihedral(Atom *atom1, Atom *atom2, Atom *atom3, Atom *atom4)
{
    AtomVector atoms {atom1, atom2, atom3, atom4};
    this->Initialize(atoms);
}

Rotatable_dihedral::Rotatable_dihedral(AtomVector atoms)
{
    this->Initialize(atoms);
}

Rotatable_dihedral::Rotatable_dihedral(AtomVector atoms, AtomVector atoms_that_move)
{
    this->SetAtoms(atoms);
    this->SetAtomsThatMove(atoms_that_move);
}

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

double Rotatable_dihedral::GetPreviousDihedralAngle()
{
    return previous_dihedral_angle_;
}

DihedralAngleDataVector Rotatable_dihedral::GetMetadata()
{
    return assigned_metadata_;
}

//////////////////////////////////////////////////////////
//                       MUTATOR                        //
//////////////////////////////////////////////////////////

void Rotatable_dihedral::DetermineAtomsThatMove()
{
    // In keeping with giving residues as GlcNAc1-4Gal, and wanting the moving atoms to be in the opposite direction:
    AtomVector atoms_that_move;
    atoms_that_move.push_back(atom3_);
    atom2_->FindConnectedAtoms(atoms_that_move);
    this->SetAtomsThatMove(atoms_that_move);
}


// Only this function uses radians. Everything else, in and out, should be degrees.
void Rotatable_dihedral::SetDihedralAngle(double dihedral_angle)
{
    GeometryTopology::Coordinate* a1 = atom4_->GetCoordinate();
    GeometryTopology::Coordinate* a2 = atom3_->GetCoordinate();
    GeometryTopology::Coordinate* a3 = atom2_->GetCoordinate();
    GeometryTopology::Coordinate* a4 = atom1_->GetCoordinate();

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

    //std::cout << "For " << atom1_->GetId() << ":"  << atom2_->GetId() << ":"  << atom3_->GetId() << ":"  << atom4_->GetId() <<  ".\nMoving: ";
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
    //std::cout << std::endl;
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

void Rotatable_dihedral::SetDihedralAngleUsingMetadata(bool use_ranges)
{
    if (assigned_metadata_.empty())
    {
        std::cout << "Error in Rotatable_dihedral::SetDihedralAngleUsingMetadata; no metadata has been set for:.\n"
                  << atom1_->GetId() << atom2_->GetId() << atom3_->GetId() << atom4_->GetId();
    }
    else if(assigned_metadata_.size() == 1)
    {
        for (const auto& entry : assigned_metadata_) // Some dihedral angles have only one rotamer i.e. one angle with a set of ranges. e.g 60 +/- 20.
        {
            double lower = entry.default_angle_value_;
            double upper = entry.default_angle_value_;
            if (use_ranges)
            {
                lower = (entry.default_angle_value_ - entry.lower_deviation_) ;
                upper = (entry.default_angle_value_ + entry.upper_deviation_) ;
            }
            this->RandomizeDihedralAngleWithinRange(lower, upper);
        }
    }
    else if(assigned_metadata_.size() >= 2) // Some dihedral angles have multiple rotamers, thus mulitple ranges to select from. e.g -60, 60, 180
    {
        std::vector<std::pair<double,double>> ranges;
        for (const auto& entry : assigned_metadata_)
        {
            double lower = entry.default_angle_value_;
            double upper = entry.default_angle_value_;
            if (use_ranges)
            {
                lower = (entry.default_angle_value_ - entry.lower_deviation_) ;
                upper = (entry.default_angle_value_ + entry.upper_deviation_) ;
            }
            ranges.emplace_back(lower, upper);
        }
        this->RandomizeDihedralAngleWithinRanges(ranges); // Pass all of the ranges into the function. It will randomly select an angle from within the ranges.
    }
    return;
}

//////////////////////////////////////////////////////////
//                  PRIVATE FUNCTIONS                   //
//////////////////////////////////////////////////////////

void Rotatable_dihedral::Initialize(AtomVector atoms)
{
    this->SetAtoms(atoms);
    this->DetermineAtomsThatMove();
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

void Rotatable_dihedral::RecordPreviousDihedralAngle(double dihedral_angle)
{
    previous_dihedral_angle_ = dihedral_angle;
}

void Rotatable_dihedral::UpdateAtomsIfPsi()
{
    for(auto &entry : assigned_metadata_)
    {
        // If it's a psi angle and is supposed to be defined by a H...
        if ((entry.dihedral_angle_name_.compare("psi")==0) && (entry.atom4_.at(0)=='H'))
        {// Find the neighbor of current atom3 which is a hydrogen, and set it to be the new atom4;
            for(auto &neighbor : atom3_->GetNode()->GetNodeNeighbors())
            {
                if(neighbor->GetName().at(0)=='H')
                {
                    std::cout << "In ";
                    this->Print();
                    std::cout << "Replaced atom4_ with " << neighbor->GetId() << "\n";
                    atom4_ = neighbor;
                }
            }
        }
    }
    return;
}


//////////////////////////////////////////////////////////
//                       DISPLAY FUNCTION               //
//////////////////////////////////////////////////////////

void Rotatable_dihedral::Print()
{
    std::cout << atom1_->GetName() << ", " << atom2_->GetName() << ", " << atom3_->GetName() << ", " << atom4_->GetName() << ": " << this->CalculateDihedralAngle()  << ".\n";
//    for(AtomVector::iterator it1 = atoms_that_move_.begin(); it1 != atoms_that_move_.end(); ++it1)
//    {
//        Atom *atom = *it1;
//        std::cout << atom->GetName() << ", ";
//    }
    std::cout << std::endl;
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
