#include "includes/CentralDataStructure/Shapers/rotatableDihedral.hpp"
#include "includes/CentralDataStructure/Shapers/atomToCoordinateInterface.hpp" // getCoordinatesFromAtoms
#include "includes/CentralDataStructure/Shapers/shapers.hpp"
#include "includes/CentralDataStructure/Selections/atomSelections.hpp" //FindConnectedAtoms
#include "includes/CentralDataStructure/Measurements/measurements.hpp"

using cds::RotatableDihedral;

RotatableDihedral::RotatableDihedral(cds::Atom* atom1, cds::Atom* atom2, cds::Atom* atom3, cds::Atom* atom4,
                                     bool reverseAtomsThatMove)
{
    std::vector<cds::Atom*> atoms {atom1, atom2, atom3, atom4};
    this->Initialize(atoms, reverseAtomsThatMove);
}

RotatableDihedral::RotatableDihedral(cds::Atom* atom1, cds::Atom* atom2, cds::Atom* atom3, cds::Atom* atom4,
                                     std::vector<cds::Atom*> extraAtomsThatMove, bool reverseAtomsThatMove)
{
    std::vector<cds::Atom*> atoms {atom1, atom2, atom3, atom4};
    this->Initialize(atoms, reverseAtomsThatMove);
    this->AddExtraAtomsThatMove(extraAtomsThatMove);
}

double RotatableDihedral::CalculateDihedralAngle(const std::string type) const
{
    if (type == "glycamReport")
    {
        if ((this->GetName() == "Phi") && (atom1_->getName() == "C1"))
        { // ToDO there will be an issue with C1-O5 linkages unless atom knows which residue it is in.
            // Coordinate* o5Coord = atom1_->GetResidue()->GetAtom("O5")->getCoordinate();
            Coordinate* o5Coord = cdsSelections::getNeighborNamed(atom1_, "O5")->getCoordinate();
            return cds::CalculateDihedralAngle(o5Coord, atom2_->getCoordinate(), atom3_->getCoordinate(),
                                               atom4_->getCoordinate());
        }
    }
    return cds::CalculateDihedralAngle(atom1_->getCoordinate(), atom2_->getCoordinate(), atom3_->getCoordinate(),
                                       atom4_->getCoordinate());
}

std::vector<cds::Atom*> RotatableDihedral::GetAtoms() const
{
    std::vector<cds::Atom*> atoms = {atom1_, atom2_, atom3_, atom4_};
    return atoms;
}

int RotatableDihedral::GetNumberOfRotamers(bool likelyShapesOnly) const
{
    if (this->GetMetadata().empty())
    {
        gmml::log(__LINE__, __FILE__, gmml::ERR,
                  "Error in RotatableDihedral::GetNumberOfRotamers; no metadata has been set.\n");
        return 0;
    }
    else
    {
        int count = 0;
        for (auto& entry : this->GetMetadata())
        {
            if ((entry.weight_ < 0.01) && (likelyShapesOnly)) // I'm hardcoding it, I don't care.
            {}                                                // Do nought.
            else
            {
                count++;
            }
        }
        return count;
    }
}

DihedralAngleDataVector RotatableDihedral::GetLikelyMetadata() const
{
    DihedralAngleDataVector returningMetadata;
    for (auto& entry : this->GetMetadata())
    {
        if (entry.weight_ >= 0.01) // HARDCODE EVERYTHING.
        {
            // std::cout << "Likely entry for: " << entry.dihedral_angle_name_ << " weigh: " << entry.weight_ << "\n";
            returningMetadata.push_back(entry);
        }
        // else {std::cout << "UnLikely entry for: " << entry.dihedral_angle_name_ << " weigh: " << entry.weight_ <<
        // "\n";}
    }
    return returningMetadata;
}

std::vector<double> RotatableDihedral::GetAllPossibleAngleValues(const int interval) const
{
    std::vector<double> allPossibleAngleValues;
    if (assigned_metadata_.empty())
    {
        std::stringstream ss;
        ss << "Error in RotatableDihedral::GetAllPossibleAngleValues; no metadata has been set for:\n"
           << atom1_->getId() << " " << atom2_->getId() << " " << atom3_->getId() << " " << atom4_->getId() << "\n";
        gmml::log(__LINE__, __FILE__, gmml::ERR, ss.str());
        throw std::runtime_error(ss.str());
    }
    else
    {
        gmml::MolecularMetadata::GLYCAM::DihedralAngleDataVector metadata_entries = this->GetMetadata();
        // std::cout << metadata_entries.at(0).dihedral_angle_name_;
        for (auto& metadata : metadata_entries)
        {
            double lower_bound      = (metadata.default_angle_value_ - metadata.lower_deviation_);
            double upper_bound      = (metadata.default_angle_value_ + metadata.upper_deviation_);
            double current_dihedral = lower_bound;
            while (current_dihedral <= upper_bound)
            {
                allPossibleAngleValues.push_back(current_dihedral);
                current_dihedral += interval; // increment
            }
        }
    }
    return allPossibleAngleValues;
}

std::string RotatableDihedral::GetName() const
{
    if (this->GetLikelyMetadata().empty())
    {
        return "Boo";
    }
    else
    {
        return this->GetLikelyMetadata().at(0).dihedral_angle_name_;
    }
}

void RotatableDihedral::DetermineAtomsThatMove()
{
    // In keeping with giving residues as GlcNAc1-4Gal, and wanting the moving atoms to be in the opposite direction
    // (defaults):
    std::vector<cds::Atom*> atoms_that_move;
    if (this->GetIsAtomsThatMoveReversed())
    {
        //        std::cout << "Blocking access via " << atom3_->getName() << ", we will search outward from " <<
        //        atom2_->getName() << "\n";
        atoms_that_move.push_back(atom3_);
        cdsSelections::FindConnectedAtoms(atoms_that_move, atom2_);
    }
    else
    {
        //        std::cout << "Blocking access via " << atom2_->getName() << ", we will search outward from " <<
        //        atom3_->getName() << "\n";
        atoms_that_move.push_back(atom2_);
        cdsSelections::FindConnectedAtoms(atoms_that_move, atom3_);
    }
    //    std::cout << "Have determined that the following " << this->GetName() << " angle atoms will move:\n";
    // for (auto& atom : atoms_that_move)
    //{
    //        std::cout << atom->getName() << ", ";
    //}
    //    std::cout << "\n";
    this->SetAtomsThatMove(atoms_that_move);
}

void RotatableDihedral::SetAtomsThatMove(std::vector<cds::Atom*>& atoms)
{
    coordinatesThatMove_.clear();
    // coordinatesThatMove_.reserve(atoms.size());
    this->AddExtraAtomsThatMove(atoms);
    return;
}

void RotatableDihedral::AddExtraAtomsThatMove(std::vector<cds::Atom*>& extraAtoms)
{
    for (auto& atom : extraAtoms)
    {
        coordinatesThatMove_.push_back(atom->getCoordinate());
    }
    return;
}

void RotatableDihedral::SetDihedralAngle(const double dihedral_angle)
{
    if (!this->wasEverRotated_)
    {
        this->SetWasEverRotated(true);
        this->DetermineAtomsThatMove();
        // std::cout << "Rotatabledihedral: " << this->getIndex() << " was never rotated\n";
    }
    Coordinate* a1 = atom1_->getCoordinate();
    Coordinate* a2 = atom2_->getCoordinate();
    Coordinate* a3 = atom3_->getCoordinate();
    Coordinate* a4 = atom4_->getCoordinate();
    if (this->GetIsAtomsThatMoveReversed())
    { // A function that was agnostic of this would be great.
        a1 = atom4_->getCoordinate();
        a2 = atom3_->getCoordinate();
        a3 = atom2_->getCoordinate();
        a4 = atom1_->getCoordinate();
    }
    this->RecordPreviousDihedralAngle(this->CalculateDihedralAngle());
    // gmml::log(__LINE__,__FILE__,gmml::INF, "Here my dude");
    // std::stringstream ss;
    // ss << "Setting dihedral for: " << atom1_->getId() << " "  << atom2_->getId() << " "  << atom3_->getId() << " "
    //            << atom4_->getId() <<  "  " << dihedral_angle << "\n";
    // gmml::log(__LINE__,__FILE__,gmml::INF, ss.str());
    cds::SetDihedralAngle(a1, a2, a3, a4, dihedral_angle, this->GetCoordinatesThatMove());
    return;
}

void RotatableDihedral::SetDihedralAngleToPrevious()
{
    this->SetDihedralAngle(this->GetPreviousDihedralAngle());
}

double RotatableDihedral::RandomizeDihedralAngle()
{
    return RotatableDihedral::RandomizeDihedralAngleWithinRange(0.0, 360.0);
}

double RotatableDihedral::RandomizeDihedralAngleWithinRange(double min, double max)
{
    std::uniform_real_distribution<> angle_distribution(min, max); // define the range
    double random_angle = angle_distribution(rng);
    /*******************************************/
    /*               IMPORTANT                 */
    /*******************************************/

    this->SetDihedralAngle(random_angle); // THIS IS IMPORTANT!!! THIS SHOULD BE SEPARATED?!?! The two other functions
                                          // call this one. Seems fragile.

    /*******************************************/
    /*               IMPORTANT                 */
    /*******************************************/
    return random_angle;
}

double RotatableDihedral::RandomizeDihedralAngleWithinRanges(std::vector<std::pair<double, double>> ranges)
{
    std::uniform_int_distribution<> distr(0, (ranges.size() - 1)); // define the range
    // Select one of the ranges
    int range_selection = distr(rng);
    // create an angle within the selected range
    return this->RandomizeDihedralAngleWithinRange(ranges.at(range_selection).first, ranges.at(range_selection).second);
}

void RotatableDihedral::AddMetadata(DihedralAngleData metadata)
{
    assigned_metadata_.push_back(metadata);
    return;
}

void RotatableDihedral::SetRandomAngleEntryUsingMetadata(bool useRanges)
{
    if (assigned_metadata_.empty())
    {
        std::stringstream ss;
        ss << "Error in RotatableDihedral::SetRandomAngleEntryUsingMetadata; no metadata has been set for:.\n"
           << atom1_->getId() << atom2_->getId() << atom3_->getId() << atom4_->getId();
        gmml::log(__LINE__, __FILE__, gmml::ERR, ss.str());
        throw std::runtime_error(ss.str());
    }
    else
    { // first randomly pick one of the meta data entries. If there is only one entry, it will randomly always be it.
        std::uniform_int_distribution<> distr(0, (assigned_metadata_.size() - 1)); // define the range
        DihedralAngleData& entry = assigned_metadata_.at(distr(rng));
        this->SetCurrentMetaData(entry);
        if (useRanges)
        {
            double lower = (entry.default_angle_value_ - entry.lower_deviation_);
            double upper = (entry.default_angle_value_ + entry.upper_deviation_);
            this->RandomizeDihedralAngleWithinRange(lower, upper);
        }
        else
        { // Just set it to the entries default value
            this->RandomizeDihedralAngleWithinRange(entry.default_angle_value_, entry.default_angle_value_);
        }
    }
    return;
}

void RotatableDihedral::SetSpecificAngleEntryUsingMetadata(bool useRanges, long unsigned int angleEntryNumber)
{
    if (assigned_metadata_.empty())
    {
        std::stringstream ss;
        ss << "Error in RotatableDihedral::SetSpecificAngleUsingMetadata; no metadata has been set for:\n"
           << atom1_->getId() << " " << atom2_->getId() << " " << atom3_->getId() << " " << atom4_->getId() << "\n";
        gmml::log(__LINE__, __FILE__, gmml::ERR, ss.str());
        throw std::runtime_error(ss.str());
    }
    else if (assigned_metadata_.size() <= angleEntryNumber)
    {
        std::stringstream ss;
        ss << "Error in RotatableDihedral::SetSpecificAngleUsingMetadata; angleEntryNumber of " << angleEntryNumber
           << " is too large as metadata.size() is " << assigned_metadata_.size() << ".\n";
        gmml::log(__LINE__, __FILE__, gmml::ERR, ss.str());
        throw std::runtime_error(ss.str());
    }
    else
    {
        DihedralAngleData& entry = assigned_metadata_.at(angleEntryNumber);
        this->SetCurrentMetaData(entry);
        if (useRanges)
        {
            double lower = (entry.default_angle_value_ - entry.lower_deviation_);
            double upper = (entry.default_angle_value_ + entry.upper_deviation_);
            this->RandomizeDihedralAngleWithinRange(lower, upper);
            //            std::stringstream ss;
            //            ss << entry.dihedral_angle_name_ << " was set to " <<  entry.default_angle_value_ << "+/-" <<
            //            entry.lower_deviation_ << ". Atoms: "
            //                    << atom1_->getId() << " " << atom2_->getId() << " " << atom3_->getId() << " " <<
            //                    atom4_->getId() << "\n";
            //            gmml::log(__LINE__,__FILE__,gmml::INF, ss.str());
        }
        else
        {
            this->SetDihedralAngle(entry.default_angle_value_);
            //            std::stringstream ss;
            //            ss << entry.dihedral_angle_name_ << " was set to " <<  entry.default_angle_value_ << ". Atoms:
            //            "
            //                    << atom1_->getId() << " " << atom2_->getId() << " " << atom3_->getId() << " " <<
            //                    atom4_->getId() << "\n";
            //            gmml::log(__LINE__,__FILE__,gmml::INF, ss.str());
        }
    }
    return;
}

bool RotatableDihedral::SetSpecificShape(std::string dihedralName, std::string selectedRotamer)
{
    if (assigned_metadata_.empty())
    {
        std::string message =
            "Error in RotatableDihedral::SetSpecificAngleUsingMetadata; no metadata has been set for:\n" +
            atom1_->getId() + " " + atom2_->getId() + " " + atom3_->getId() + " " + atom4_->getId();
        message += ".\nA rotatable dihedral named: " + dihedralName +
                   " was requested to be set to the rotamer named: " + selectedRotamer;
        gmml::log(__LINE__, __FILE__, gmml::ERR, message);
        throw std::runtime_error(message);
    }
    // gmml::log(__LINE__,__FILE__,gmml::INF, "Made it here with " + dihedralName + " and " + selectedRotamer);
    if (dihedralName == this->GetMetadata().at(0).dihedral_angle_name_)
    {
        for (auto& metadata : this->GetMetadata())
        {
            //                std::stringstream ss;
            //                ss << dihedralName << ": " <<  selectedRotamer <<  " vs metadata " <<
            //                metadata.dihedral_angle_name_ << ": " << metadata.rotamer_name_;
            //                gmml::log(__LINE__,__FILE__,gmml::INF, ss.str());
            if (metadata.rotamer_name_ == selectedRotamer)
            {
                this->SetDihedralAngle(metadata.default_angle_value_);
                this->SetCurrentMetaData(metadata);
                //                    ss << "Bingo! Setting " << dihedralName << " to " <<
                //                    metadata.default_angle_value_ << " for " << atom1_->getId() << " " <<
                //                    atom2_->getId() << " " << atom3_->getId() << " " << atom4_->getId() << "\n";
                //                    gmml::log(__LINE__,__FILE__,gmml::INF, ss.str());
                return true;
            }
        }
    }

    return false;
}

void RotatableDihedral::WiggleUsingAllRotamers(std::vector<cds::Residue*>& overlapSet1,
                                               std::vector<cds::Residue*>& overlapSet2, const int& angleIncrement)
{

    double bestDihedral        = this->CalculateDihedralAngle();
    unsigned int lowestOverlap = cds::CountOverlappingAtoms(overlapSet2, overlapSet1);
    for (auto& metadata : this->GetMetadata())
    {
        double lowerBound = (metadata.default_angle_value_ - metadata.lower_deviation_);
        double upperBound = (metadata.default_angle_value_ + metadata.upper_deviation_);
        unsigned int newOverlap =
            this->WiggleWithinRangesDistanceCheck(overlapSet1, overlapSet2, angleIncrement, lowerBound, upperBound);
        if (lowestOverlap >= (newOverlap + 1))
        {
            lowestOverlap = newOverlap;
            bestDihedral  = this->CalculateDihedralAngle();
            this->SetCurrentMetaData(metadata);
        }
        else
        {
            this->SetDihedralAngle(bestDihedral);
        }
    }
    return;
}

void RotatableDihedral::WiggleUsingAllRotamers(std::vector<cds::Atom*>& overlapAtomSet1,
                                               std::vector<cds::Atom*>& overlapAtomSet2, const int& angleIncrement)
{
    double bestDihedral        = this->CalculateDihedralAngle();
    unsigned int lowestOverlap = cds::CountOverlappingAtoms(overlapAtomSet1, overlapAtomSet2);
    ;
    for (auto& metadata : this->GetMetadata())
    {
        double lowerBound       = (metadata.default_angle_value_ - metadata.lower_deviation_);
        double upperBound       = (metadata.default_angle_value_ + metadata.upper_deviation_);
        unsigned int newOverlap = this->WiggleWithinRangesDistanceCheck(overlapAtomSet1, overlapAtomSet2,
                                                                        angleIncrement, lowerBound, upperBound);
        if (lowestOverlap >= (newOverlap + 1))
        {
            lowestOverlap = newOverlap;
            bestDihedral  = this->CalculateDihedralAngle();
            this->SetCurrentMetaData(metadata);
        }
        else
        {
            this->SetDihedralAngle(bestDihedral);
        }
    }
    return;
}

// User requested gg, this prevents flipping into gt like the above would do. i.e. cb won't want a flip, gp would.
void RotatableDihedral::WiggleWithinCurrentRotamer(std::vector<cds::Atom*>& overlapAtomSet1,
                                                   std::vector<cds::Atom*>& overlapAtomSet2, const int& angleIncrement)
{
    double lowerBound =
        (this->GetCurrentMetaData()->default_angle_value_ - this->GetCurrentMetaData()->lower_deviation_);
    double upperBound =
        (this->GetCurrentMetaData()->default_angle_value_ + this->GetCurrentMetaData()->upper_deviation_);
    // Set to lowest deviation, work through to highest. Set best value and return it for reference.
    this->WiggleWithinRangesDistanceCheck(overlapAtomSet1, overlapAtomSet2, angleIncrement, lowerBound, upperBound);
    return;
}

// User requested gg, this prevents flipping into gt like the above would do.
void RotatableDihedral::WiggleWithinCurrentRotamer(std::vector<cds::Residue*>& overlapResidueSet1,
                                                   std::vector<cds::Residue*>& overlapResidueSet2,
                                                   const int& angleIncrement)
{
    double lowerBound =
        (this->GetCurrentMetaData()->default_angle_value_ - this->GetCurrentMetaData()->lower_deviation_);
    double upperBound =
        (this->GetCurrentMetaData()->default_angle_value_ + this->GetCurrentMetaData()->upper_deviation_);
    // Set to lowest deviation, work through to highest. Set best value and return it for reference.
    this->WiggleWithinRangesDistanceCheck(overlapResidueSet1, overlapResidueSet2, angleIncrement, lowerBound,
                                          upperBound);
    return;
}

// double RotatableDihedral::WiggleWithinRanges(std::vector<cds::Atom*>& overlapAtomSet1,
//                                              std::vector<cds::Atom*>& overlapAtomSet2, const int& angleIncrement,
//                                              const double& lowerBound, const double& upperBound)
//{
//     this->SetDihedralAngle(lowerBound);
//     double currentDihedral = lowerBound;
//     double bestDihedral    = lowerBound;
//     unsigned int currentOverlap  = cds::CalculateAtomicOverlapsBetweenNonBondedAtoms(overlapAtomSet1,
//     overlapAtomSet2); unsigned int lowestOverlap   = currentOverlap;
//     //    std::cout << "Starting wiggle with overlap : " << currentOverlap << " . Current Angle: " << currentDihedral
//     <<
//     //    " lowerBound: " << lowerBound << " upperBound: " << upperBound << "\n";
//     while (currentDihedral < upperBound)
//     {
//         currentDihedral += angleIncrement; // increment
//         this->SetDihedralAngle(currentDihedral);
//         currentOverlap = cds::CalculateAtomicOverlapsBetweenNonBondedAtoms(overlapAtomSet1, overlapAtomSet2);
//         if (lowestOverlap >= (currentOverlap + 0.01)) // rounding errors
//         {
//             lowestOverlap = currentOverlap;
//             bestDihedral  = currentDihedral;
//         }
//         // Prefer angles closer to default if overlap is the same.
//         else if ((lowestOverlap == currentOverlap) &&
//                  (std::abs(this->GetCurrentMetaData()->default_angle_value_ - bestDihedral) >
//                   std::abs(this->GetCurrentMetaData()->default_angle_value_ - currentDihedral)))
//         {
//             bestDihedral = currentDihedral;
//         }
//     }
//     this->SetDihedralAngle(bestDihedral);
//     return lowestOverlap;
// }

unsigned int RotatableDihedral::WiggleWithinRangesDistanceCheck(std::vector<cds::Atom*>& overlapAtomSet1,
                                                                std::vector<cds::Atom*>& overlapAtomSet2,
                                                                const int& angleIncrement, const double& lowerBound,
                                                                const double& upperBound)
{
    this->SetDihedralAngle(lowerBound);
    double currentDihedral      = lowerBound;
    double bestDihedral         = lowerBound;
    unsigned int currentOverlap = cds::CountOverlappingAtoms(overlapAtomSet1, overlapAtomSet2);
    unsigned int lowestOverlap  = currentOverlap;
    //    std::cout << "Starting wiggle with overlap : " << currentOverlap << " . Current Angle: " << currentDihedral <<
    //    " lowerBound: " << lowerBound << " upperBound: " << upperBound << "\n";
    while (currentDihedral < upperBound)
    {
        currentDihedral += angleIncrement; // increment
        this->SetDihedralAngle(currentDihedral);
        currentOverlap = cds::CountOverlappingAtoms(overlapAtomSet1, overlapAtomSet2);
        if (lowestOverlap >= (currentOverlap + 1))
        {
            lowestOverlap = currentOverlap;
            bestDihedral  = currentDihedral;
        }
        // Prefer angles closer to default if overlap is the same.
        else if ((lowestOverlap == currentOverlap) &&
                 (std::abs(this->GetCurrentMetaData()->default_angle_value_ - bestDihedral) >
                  std::abs(this->GetCurrentMetaData()->default_angle_value_ - currentDihedral)))
        {
            bestDihedral = currentDihedral;
        }
    }
    this->SetDihedralAngle(bestDihedral);
    return lowestOverlap;
}

unsigned int RotatableDihedral::WiggleWithinRangesDistanceCheck(std::vector<cds::Residue*>& overlapResidueSet1,
                                                                std::vector<cds::Residue*>& overlapResidueSet2,
                                                                const int& angleIncrement, const double& lowerBound,
                                                                const double& upperBound)
{
    this->SetDihedralAngle(lowerBound);
    double currentDihedral      = lowerBound;
    double bestDihedral         = lowerBound;
    unsigned int currentOverlap = cds::CountOverlappingAtoms(overlapResidueSet1, overlapResidueSet2);
    unsigned int lowestOverlap  = currentOverlap;
    // std::stringstream ss;
    // ss << "Starting wiggle with overlap : " << currentOverlap << " . Current Angle: " << currentDihedral <<
    // " lowerBound: " << lowerBound << " upperBound: " << upperBound << "\n";
    // gmml::log(__LINE__,__FILE__,gmml::INF, ss.str());
    if (this->GetCurrentMetaData() == nullptr) // you don't need these checks if you have RAII OLIVER
    {
        std::string message = "Error: current metadata not set for RotatableDihedral " + this->GetName();
        gmml::log(__LINE__, __FILE__, gmml::ERR, message);
        throw std::runtime_error(message);
    }
    while (currentDihedral < upperBound)
    {
        currentDihedral += angleIncrement; // increment
        this->SetDihedralAngle(currentDihedral);
        currentOverlap = cds::CountOverlappingAtoms(overlapResidueSet1, overlapResidueSet2);
        if (lowestOverlap >= (currentOverlap + 1))
        {
            lowestOverlap = currentOverlap;
            bestDihedral  = currentDihedral;
        }
        // Prefer angles closer to default if overlap is the same.
        else if ((lowestOverlap == currentOverlap) &&
                 (std::abs(this->GetCurrentMetaData()->default_angle_value_ - bestDihedral) >
                  std::abs(this->GetCurrentMetaData()->default_angle_value_ - currentDihedral)))
        {
            bestDihedral = currentDihedral;
        }
    }
    this->SetDihedralAngle(bestDihedral);
    return lowestOverlap;
}

// Private

bool RotatableDihedral::IsThereHydrogenForPsiAngle()
{
    for (auto& neighbor : atom3_->getNeighbors())
    {
        if (neighbor->getName().at(0) == 'H')
        {
            // std::cout << "In "; this->Print(); std::cout << "Replaced atom4_ with " << neighbor->getId() << "\n";
            atom4_ = neighbor;
            return true;
        }
    }
    return false;
}

std::unique_ptr<cds::Atom> RotatableDihedral::CreateHydrogenAtomForPsiAngle()
{
    std::vector<Coordinate*> neighborsCoords = cds::getCoordinatesFromAtoms(atom3_->getNeighbors());
    Coordinate newCoord = cds::CreateCoordinateForCenterAwayFromNeighbors(atom3_->getCoordinate(), neighborsCoords);
    std::unique_ptr<cds::Atom> newAtom = std::make_unique<cds::Atom>("HHH", newCoord);
    atom3_->addBond(newAtom.get());
    return newAtom;
}

void RotatableDihedral::Initialize(std::vector<cds::Atom*> atoms, bool reverseAtomsThatMove)
{
    this->SetAtoms(atoms);
    this->SetIsAtomsThatMoveReversed(reverseAtomsThatMove);
    currentMetadata_ = nullptr; // When does this get set? This is bad.
    // this->DetermineAtomsThatMove(); // Will be done the first time SetDihedralAngle is called.
    return;
}

void RotatableDihedral::SetAtoms(std::vector<cds::Atom*> atoms)
{
    atom1_ = atoms.at(0);
    atom2_ = atoms.at(1);
    atom3_ = atoms.at(2);
    atom4_ = atoms.at(3);
    return;
}

std::string RotatableDihedral::Print() const
{
    std::stringstream ss;
    ss << atom1_->getName() << ", " << atom2_->getName() << ", " << atom3_->getName() << ", " << atom4_->getName()
       << ": " << this->CalculateDihedralAngle() << ".\n";
    gmml::log(__LINE__, __FILE__, gmml::INF, ss.str());
    return ss.str();
}
