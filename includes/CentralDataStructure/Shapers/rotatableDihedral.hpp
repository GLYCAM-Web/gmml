#ifndef GMML_INCLUDES_CENTRAL_DATA_STRUCTURE_SHAPERS_ROTATABLE_DIHEDRAL_HPP
#define GMML_INCLUDES_CENTRAL_DATA_STRUCTURE_SHAPERS_ROTATABLE_DIHEDRAL_HPP

#include "includes/MolecularMetadata/GLYCAM/dihedralangledata.hpp"
#include "includes/GeometryTopology/geometrytopology.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/External_Libraries/PCG/pcg_random.h"
#include "includes/CentralDataStructure/Overlaps/cdsOverlaps.hpp"

#include <random>

using gmml::MolecularMetadata::GLYCAM::DihedralAngleData;
using gmml::MolecularMetadata::GLYCAM::DihedralAngleDataVector;
using GeometryTopology::Coordinate;

// Seed with a real random value, if available
static pcg_extras::seed_seq_from<std::random_device> seed_source;
// Make a random number engine
static pcg32 rng(seed_source); // Eclipse complains about ambiguity, and yet it compiles...

// This class stores the four atoms that define a dihedral angle, the atoms that move when it is rotated
// and, if moved, the previous dihedral angle, which allows me to reset easily.
namespace cds
{
template<class atomT>
class RotatableDihedral
{
public:
    //////////////////////////////////////////////////////////
    //                       CONSTRUCTOR                    //
    //////////////////////////////////////////////////////////
    RotatableDihedral(atomT* atom1, atomT* atom2, atomT* atom3, atomT* atom4, bool reverseAtomsThatMove = true);
    RotatableDihedral(atomT* atom1, atomT* atom2, atomT* atom3, atomT* atom4, std::vector<atomT*> extraAtomsThatMove, bool reverseAtomsThatMove = true);
    //////////////////////////////////////////////////////////
    //                       ACCESSOR                       //
    //////////////////////////////////////////////////////////
    inline const DihedralAngleDataVector& GetMetadata() const {return assigned_metadata_;}
    DihedralAngleDataVector GetLikelyMetadata() const;
    int GetNumberOfRotamers( bool likelyShapesOnly = false) const;
    std::string GetName() const;
    double CalculateDihedralAngle(const std::string type = "default") const;
    //////////////////////////////////////////////////////////
    //                       MUTATOR                        //
    //////////////////////////////////////////////////////////
    void DetermineAtomsThatMove(); // Based on connectivities, this figures out which atoms will move when the dihedral is rotated.
    void SetDihedralAngle(double dihedral_angle); // Sets the dihedral angle by rotating the bond between atom2 and atom3, moving atom4 and connected.
    void SetDihedralAngleToPrevious(); // Sets the dihedral to previous dihedral angle
    double RandomizeDihedralAngle(); // Randomly sets dihedral angle values between 0 and 360
    void AddMetadata(DihedralAngleData metadata);
    inline void ClearMetadata() {assigned_metadata_.clear();}
    void SetRandomAngleEntryUsingMetadata(bool useRanges = true);
    void SetSpecificAngleEntryUsingMetadata(bool useRanges = false, int angleEntryNumber = 0);
    bool SetSpecificShape(std::string dihedralName, std::string selectedRotamer);
    void WiggleWithinCurrentRotamer(std::vector<atomT*> &overlapAtomSet1, std::vector<atomT*> &overlapAtomSet2, const int &angleIncrement);
    void WiggleUsingAllRotamers(std::vector<atomT*>& overlapAtomSet1, std::vector<atomT*> &overlapAtomSet2, const int &angleIncrement);
    //////////////////////////////////////////////////////////
    //                       DISPLAY FUNCTION               //
    //////////////////////////////////////////////////////////
    std::string Print() const;
private:
    //////////////////////////////////////////////////////////
    //                  PRIVATE ACCESSORS                   //
    //////////////////////////////////////////////////////////
    std::vector<atomT*> GetAtoms() const;
    inline std::vector<atomT*> GetAtomsThatMove() const {return atoms_that_move_;}
    inline bool GetIsAtomsThatMoveReversed() const {return isAtomsThatMoveReversed_;}
    inline double GetPreviousDihedralAngle() const {return previous_dihedral_angle_;}
    std::vector<double> GetAllPossibleAngleValues(const int interval = 5) const;
    //////////////////////////////////////////////////////////
    //                  PRIVATE MUTATORS                    //
    //////////////////////////////////////////////////////////
    void AddExtraAtomsThatMove(std::vector<atomT*> extraAtoms);
    double RandomizeDihedralAngleWithinRanges(std::vector<std::pair<double,double> > ranges); // Takes in a set of ranges, e.g. 10 to 30, 45-55 etc. Randomly selects a range and randomly sets value within that range.
    double RandomizeDihedralAngleWithinRange(double min, double max);// Randomly sets dihedral angle to a value within the given range. E.g. Between 25 and 30 degrees.
    //////////////////////////////////////////////////////////
    //                  PRIVATE FUNCTIONS                   //
    //////////////////////////////////////////////////////////
    void Initialize(std::vector<atomT*> atoms, bool reverseAtomsThatMove = true);
    void SetAtoms(std::vector<atomT*> atoms);
    inline void SetAtomsThatMove(std::vector<atomT*> atoms) {atoms_that_move_ = atoms;}
    inline void SetIsAtomsThatMoveReversed(bool b) {isAtomsThatMoveReversed_ = b;}
    inline void RecordPreviousDihedralAngle(double d) { previous_dihedral_angle_ = d;}
    inline void SetWasEverRotated(bool b) {wasEverRotated_ = b;}
    void UpdateAtomsIfPsi();
    atomT* CreateHydrogenAtomForPsi(atomT* centralAtom);
    inline bool CheckIfEverRotated() const {return wasEverRotated_;}
    inline void SetCurrentMetaData(const DihedralAngleData &d) {currentMetadata_ = &d;}
    inline const DihedralAngleData* GetCurrentMetaData() {return currentMetadata_;}
    double WiggleWithinRanges(std::vector<atomT*>& overlapAtomSet1, std::vector<atomT*> &overlapAtomSet2, const int &angleIncrement, const double& lowerBound, const double& upperBound);
    //////////////////////////////////////////////////////////
    //                       ATTRIBUTES                     //
    //////////////////////////////////////////////////////////
    // The four atoms that define the dihedral angle. The bond between atom2_ and atom3_ is what is rotated.
    atomT* atom1_;
    atomT* atom2_;
    atomT* atom3_;
    atomT* atom4_;
    // A vector of pointers to the atoms that are connected to atom2_ and atom3_, and will be rotated when that bond is rotated.
    std::vector<atomT*> atoms_that_move_;
    std::vector<atomT*> extra_atoms_that_move_;
    bool isAtomsThatMoveReversed_;
    double previous_dihedral_angle_; // I often want to reset a dihedral angle after rotating it, so recording the previous angle makes this easy.
    DihedralAngleDataVector assigned_metadata_;
    const DihedralAngleData* currentMetadata_;
    bool wasEverRotated_; // Need this, as it might add a H atom for psi
};

template<class atomT>
RotatableDihedral<atomT>::RotatableDihedral(atomT* atom1, atomT* atom2, atomT* atom3, atomT* atom4, bool reverseAtomsThatMove)
{
    std::vector<atomT*> atoms {atom1, atom2, atom3, atom4};
    this->Initialize(atoms, reverseAtomsThatMove);
}

template<class atomT>
RotatableDihedral<atomT>::RotatableDihedral(atomT* atom1, atomT* atom2, atomT* atom3, atomT* atom4, std::vector<atomT*> extraAtomsThatMove, bool reverseAtomsThatMove)
{
    std::vector<atomT*> atoms {atom1, atom2, atom3, atom4};
    this->Initialize(atoms, reverseAtomsThatMove);
    this->AddExtraAtomsThatMove(extraAtomsThatMove);
}

template<class atomT>
double RotatableDihedral<atomT>::CalculateDihedralAngle(const std::string type) const
{
    if (type == "glycamReport")
    {
        if ( (this->GetName() == "Phi") && (atom1_->getName() == "C1") )
        {
            //Coordinate* o5Coord = atom1_->GetResidue()->GetAtom("O5")->getCoordinate();
            Coordinate* o5Coord = atom1_->seekNeighbor("O5", atom1_->getResidueId());
            return GeometryTopology::CalculateDihedralAngle(o5Coord, atom2_->getCoordinate(), atom3_->getCoordinate(), atom4_->getCoordinate());
        }
    }
    return GeometryTopology::CalculateDihedralAngle(atom1_->getCoordinate(), atom2_->getCoordinate(), atom3_->getCoordinate(), atom4_->getCoordinate());
}

template<class atomT>
std::vector<atomT*> RotatableDihedral<atomT>::GetAtoms() const
{
    std::vector<atomT*> atoms = {atom1_, atom2_, atom3_, atom4_};
    return atoms;
}

template<class atomT>
int RotatableDihedral<atomT>::GetNumberOfRotamers(bool likelyShapesOnly) const
{
    if (this->GetMetadata().empty())
    {
        gmml::log(__LINE__,__FILE__,gmml::ERR, "Error in RotatableDihedral::GetNumberOfRotamers; no metadata has been set.\n");
        return 0;
    }
    else
    {
        int count = 0;
        for (auto & entry: this->GetMetadata())
        {
            if ( (entry.weight_ < 0.01) && (likelyShapesOnly) ) // I'm hardcoding it, I don't care.
            {}// Do nought.
            else
                count++;
        }
        return count;
    }
}

template<class atomT>
DihedralAngleDataVector RotatableDihedral<atomT>::GetLikelyMetadata() const
{
    DihedralAngleDataVector returningMetadata;
    for (auto &entry : this->GetMetadata())
    {
        if(entry.weight_ >= 0.01 ) // HARDCODE EVERYTHING.
        {
            //std::cout << "Likely entry for: " << entry.dihedral_angle_name_ << " weigh: " << entry.weight_ << "\n";
            returningMetadata.push_back(entry);
        }
        //else {std::cout << "UnLikely entry for: " << entry.dihedral_angle_name_ << " weigh: " << entry.weight_ << "\n";}
    }
    return returningMetadata;
}

template<class atomT>
std::vector<double> RotatableDihedral<atomT>::GetAllPossibleAngleValues(const int interval) const
{
    std::vector<double> allPossibleAngleValues;
    if (assigned_metadata_.empty())
    {
        std::stringstream ss;
        ss << "Error in RotatableDihedral::GetAllPossibleAngleValues; no metadata has been set for:\n"
                << atom1_->getId() << " " << atom2_->getId() << " " << atom3_->getId() << " " << atom4_->getId() << "\n";
        gmml::log(__LINE__,__FILE__,gmml::ERR, ss.str());
        throw std::runtime_error(ss.str());
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

template<class atomT>
std::string RotatableDihedral<atomT>::GetName() const
{
    if (this->GetLikelyMetadata().empty())
        return "Boo";
    else
        return this->GetLikelyMetadata().at(0).dihedral_angle_name_;
}

template<class atomT>
void RotatableDihedral<atomT>::DetermineAtomsThatMove()
{
    // In keeping with giving residues as GlcNAc1-4Gal, and wanting the moving atoms to be in the opposite direction (defaults):
    std::vector<atomT*> atoms_that_move;
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

template<class atomT>
void RotatableDihedral<atomT>::AddExtraAtomsThatMove(std::vector<atomT*> extraAtoms)
{
    extra_atoms_that_move_.insert(extra_atoms_that_move_.end(), extraAtoms.begin(), extraAtoms.end());
    return;
}

template<class atomT>
void RotatableDihedral<atomT>::SetDihedralAngle(const double dihedral_angle)
{
    if(!this->wasEverRotated_)
    {
        this->SetWasEverRotated(true);
        this->DetermineAtomsThatMove();
    }
    Coordinate* a1 = atom1_->getCoordinate();
    Coordinate* a2 = atom2_->getCoordinate();
    Coordinate* a3 = atom3_->getCoordinate();
    Coordinate* a4 = atom4_->getCoordinate();
    if(this->GetIsAtomsThatMoveReversed())
    { // A function that was agnostic of this would be great.
        a1 = atom4_->getCoordinate();
        a2 = atom3_->getCoordinate();
        a3 = atom2_->getCoordinate();
        a4 = atom1_->getCoordinate();
    }
    this->RecordPreviousDihedralAngle(this->CalculateDihedralAngle());
    std::stringstream ss;
    ss << "Setting dihedral for " << atom1_->getId() << ":"  << atom2_->getId() << ":"  << atom3_->getId() << ":"  << atom4_->getId() <<  ": " << dihedral_angle << "\n";
    gmml::log(__LINE__,__FILE__,gmml::INF, ss.str());
    std::vector<Coordinate*> coordinatesToMove;
    for(auto &atom : this->GetAtomsThatMove())
    {
        coordinatesToMove.push_back(atom->getCoordinate());
    }
    for(auto &atom : extra_atoms_that_move_)
    {
        coordinatesToMove.push_back(atom->getCoordinate());
    }
    GeometryTopology::SetDihedralAngle(a1, a2, a3, a4, dihedral_angle, coordinatesToMove);
    return;
}

template<class atomT>
void RotatableDihedral<atomT>::SetDihedralAngleToPrevious()
{
    this->SetDihedralAngle(this->GetPreviousDihedralAngle());
}

template<class atomT>
double RotatableDihedral<atomT>::RandomizeDihedralAngle()
{
    return RotatableDihedral::RandomizeDihedralAngleWithinRange(0.0, 360.0);
}

template<class atomT>
double RotatableDihedral<atomT>::RandomizeDihedralAngleWithinRange(double min, double max)
{
    std::uniform_real_distribution<> angle_distribution(min, max); // define the range
    double random_angle = angle_distribution(rng);
    /*******************************************/
    /*               IMPORTANT                 */
    /*******************************************/

    this->SetDihedralAngle(random_angle); // THIS IS IMPORTANT!!! THIS SHOULD BE SEPARATED?!?! The two other functions call this one. Seems fragile.

    /*******************************************/
    /*               IMPORTANT                 */
    /*******************************************/
    return random_angle;
}

template<class atomT>
double RotatableDihedral<atomT>::RandomizeDihedralAngleWithinRanges(std::vector<std::pair<double,double>> ranges)
{
    std::uniform_int_distribution<> distr(0, (ranges.size() - 1)); // define the range
    // Select one of the ranges
    int range_selection = distr(rng);
    // create an angle within the selected range
    return this->RandomizeDihedralAngleWithinRange(ranges.at(range_selection).first, ranges.at(range_selection).second);
}

template<class atomT>
void RotatableDihedral<atomT>::AddMetadata(DihedralAngleData metadata)
{
    assigned_metadata_.push_back(metadata);
    this->UpdateAtomsIfPsi();
    return;
}

template<class atomT>
void RotatableDihedral<atomT>::SetRandomAngleEntryUsingMetadata(bool useRanges)
{
    if (assigned_metadata_.empty())
    {
        std::stringstream ss;
        ss << "Error in RotatableDihedral::SetRandomAngleEntryUsingMetadata; no metadata has been set for:.\n"
                << atom1_->getId() << atom2_->getId() << atom3_->getId() << atom4_->getId();
        gmml::log(__LINE__,__FILE__,gmml::ERR, ss.str());
        throw std::runtime_error(ss.str());
    }
    else
    {   // first randomly pick one of the meta data entries. If there is only one entry, it will randomly always be it.
        std::uniform_int_distribution<> distr(0, (assigned_metadata_.size() - 1)); // define the range
        DihedralAngleData& entry = assigned_metadata_.at(distr(rng));
        this->SetCurrentMetaData(entry);
        if (useRanges)
        {
            double lower = (entry.default_angle_value_ - entry.lower_deviation_) ;
            double upper = (entry.default_angle_value_ + entry.upper_deviation_) ;
            this->RandomizeDihedralAngleWithinRange(lower, upper);
        }
        else
        { // Just set it to the entries default value
            this->RandomizeDihedralAngleWithinRange(entry.default_angle_value_, entry.default_angle_value_);
        }
    }
    return;
}

template<class atomT>
void RotatableDihedral<atomT>::SetSpecificAngleEntryUsingMetadata(bool useRanges, int angleEntryNumber)
{
    if (assigned_metadata_.empty())
    {
        std::stringstream ss;
        ss << "Error in RotatableDihedral::SetSpecificAngleUsingMetadata; no metadata has been set for:\n"
                << atom1_->getId() << " " << atom2_->getId() << " " << atom3_->getId() << " " << atom4_->getId() << "\n";
        gmml::log(__LINE__,__FILE__,gmml::ERR, ss.str());
        throw std::runtime_error(ss.str());
    }
    else if (assigned_metadata_.size() <= angleEntryNumber)
    {
        std::stringstream ss;
        ss << "Error in RotatableDihedral::SetSpecificAngleUsingMetadata; angleEntryNumber of " << angleEntryNumber << " is too large as metadata.size() is " << assigned_metadata_.size() << ".\n";
        gmml::log(__LINE__,__FILE__,gmml::ERR, ss.str());
        throw std::runtime_error(ss.str());
    }
    else
    {
        DihedralAngleData& entry = assigned_metadata_.at(angleEntryNumber);
        this->SetCurrentMetaData(entry);
        if (useRanges)
        {
            double lower = (entry.default_angle_value_ - entry.lower_deviation_) ;
            double upper = (entry.default_angle_value_ + entry.upper_deviation_) ;
            this->RandomizeDihedralAngleWithinRange(lower, upper);
            std::stringstream ss;
            ss << entry.dihedral_angle_name_ << " was set to " <<  entry.default_angle_value_ << "+/-" << entry.lower_deviation_ << ". Atoms: "
                    << atom1_->getId() << " " << atom2_->getId() << " " << atom3_->getId() << " " << atom4_->getId() << "\n";
            gmml::log(__LINE__,__FILE__,gmml::INF, ss.str());
        }
        else
        {
            this->SetDihedralAngle(entry.default_angle_value_);
            std::stringstream ss;
            ss << entry.dihedral_angle_name_ << " was set to " <<  entry.default_angle_value_ << ". Atoms: "
                    << atom1_->getId() << " " << atom2_->getId() << " " << atom3_->getId() << " " << atom4_->getId() << "\n";
            gmml::log(__LINE__,__FILE__,gmml::INF, ss.str());
        }
    }
    return;
}

template<class atomT>
bool RotatableDihedral<atomT>::SetSpecificShape(std::string dihedralName, std::string selectedRotamer)
{
    if (assigned_metadata_.empty())
    {
        gmml::log(__LINE__,__FILE__,gmml::ERR, "No metadata set");
    }
    else
    {
        gmml::log(__LINE__,__FILE__,gmml::INF, "Made it here with " + dihedralName + " and " + selectedRotamer);
        if (dihedralName == this->GetMetadata().at(0).dihedral_angle_name_)
        {
            for(auto &metadata : this->GetMetadata())
            {
                std::stringstream ss;
                ss << dihedralName << ": " <<  selectedRotamer <<  " vs metadata " <<  metadata.dihedral_angle_name_ << ": " << metadata.rotamer_name_;
                gmml::log(__LINE__,__FILE__,gmml::INF, ss.str());
                if (metadata.rotamer_name_ == selectedRotamer)
                {
                    this->SetDihedralAngle(metadata.default_angle_value_);
                    this->SetCurrentMetaData(metadata);
                    ss << "Bingo! Setting " << dihedralName << " to " << metadata.default_angle_value_ << " for " << atom1_->getId() << " " << atom2_->getId() << " " << atom3_->getId() << " " << atom4_->getId() << "\n";
                    gmml::log(__LINE__,__FILE__,gmml::INF, ss.str());
                    return true;
                }
            }
        }
    }
    return false;
}

template<class atomT>
void RotatableDihedral<atomT>::WiggleUsingAllRotamers(std::vector<atomT*>& overlapAtomSet1, std::vector<atomT*> &overlapAtomSet2, const int &angleIncrement)
{
    double bestDihedral = this->CalculateDihedralAngle();
    double lowestOverlap = cds::CalculateAtomicOverlapsBetweenNonBondedAtoms(overlapAtomSet1, overlapAtomSet2);;
    for(auto &metadata : this->GetMetadata())
    {
        double lowerBound = (metadata.default_angle_value_ - metadata.lower_deviation_);
        double upperBound = (metadata.default_angle_value_ + metadata.upper_deviation_);
        double newOverlap = this->WiggleWithinRanges(overlapAtomSet1, overlapAtomSet2, angleIncrement, lowerBound, upperBound);
        if (lowestOverlap >= (newOverlap + 0.01))
        {
            lowestOverlap = newOverlap;
            bestDihedral = this->CalculateDihedralAngle();
            this->SetCurrentMetaData(metadata);
        }
        else
        {
            this->SetDihedralAngle(bestDihedral);
        }
    }
    return;
}

template<class atomT>
void RotatableDihedral<atomT>::WiggleWithinCurrentRotamer(std::vector<atomT*>& overlapAtomSet1, std::vector<atomT*> &overlapAtomSet2, const int &angleIncrement)
{
    double lowerBound = (this->GetCurrentMetaData()->default_angle_value_ - this->GetCurrentMetaData()->lower_deviation_);
    double upperBound = (this->GetCurrentMetaData()->default_angle_value_ + this->GetCurrentMetaData()->upper_deviation_);
    // Set to lowest deviation, work through to highest. Set best value and return it for reference.
    this->WiggleWithinRanges(overlapAtomSet1, overlapAtomSet2, angleIncrement, lowerBound, upperBound);
    return;
}

template<class atomT>
double RotatableDihedral<atomT>::WiggleWithinRanges(std::vector<atomT*>& overlapAtomSet1, std::vector<atomT*> &overlapAtomSet2, const int &angleIncrement, const double& lowerBound, const double& upperBound)
{
    this->SetDihedralAngle(lowerBound);
    double currentDihedral = lowerBound;
    double bestDihedral = lowerBound;
    double currentOverlap = cds::CalculateAtomicOverlapsBetweenNonBondedAtoms(overlapAtomSet1, overlapAtomSet2);
    double lowestOverlap = currentOverlap;
//    std::cout << "Starting wiggle with overlap : " << currentOverlap << " . Current Angle: " << currentDihedral << " lowerBound: " << lowerBound << " upperBound: " << upperBound << "\n";
    while(currentDihedral < upperBound)
    {
        currentDihedral += angleIncrement; // increment
        this->SetDihedralAngle(currentDihedral);
        currentOverlap = cds::CalculateAtomicOverlapsBetweenNonBondedAtoms(overlapAtomSet1, overlapAtomSet2);
        if (lowestOverlap >= (currentOverlap + 0.01)) // rounding errors
        {
            lowestOverlap = currentOverlap;
            bestDihedral = currentDihedral;
        }
        // Prefer angles closer to default if overlap is the same.
        else if ( (lowestOverlap == currentOverlap) && (std::abs(this->GetCurrentMetaData()->default_angle_value_ - bestDihedral)
                > std::abs(this->GetCurrentMetaData()->default_angle_value_ - currentDihedral)) )
        {
            bestDihedral = currentDihedral;
        }
    }
    this->SetDihedralAngle(bestDihedral);
    return lowestOverlap;
}

// Private

template<class atomT>
void RotatableDihedral<atomT>::Initialize(std::vector<atomT*> atoms, bool reverseAtomsThatMove)
{
    this->SetWasEverRotated(false);
    this->SetAtoms(atoms);
    this->SetIsAtomsThatMoveReversed(reverseAtomsThatMove);
    currentMetadata_ = nullptr;
    // this->DetermineAtomsThatMove(); // Will be done the first time SetDihedralAngle is called.
    return;
}

template<class atomT>
void RotatableDihedral<atomT>::SetAtoms(std::vector<atomT*> atoms)
{
    atom1_ = atoms.at(0);
    atom2_ = atoms.at(1);
    atom3_ = atoms.at(2);
    atom4_ = atoms.at(3);
    return;
}

template<class atomT>
void RotatableDihedral<atomT>::UpdateAtomsIfPsi() // Move up to residue linkage so that we have residue info for adding H atom into.
{
    for(auto &entry : assigned_metadata_)
    {
        // If it's a psi angle and is supposed to be defined by a H...
        if ((entry.dihedral_angle_name_.compare("Psi")==0) && (entry.atom4_.at(0)=='H'))
        {// Find the neighbor of current atom3 which is a hydrogen, and set it to be the new atom4;
            bool foundHydrogen = false;
            for(auto &neighbor : atom3_->getNeighbors())
            {
                if(neighbor->getName().at(0)=='H')
                {
//                    std::cout << "In ";
//                    this->Print();
//                    std::cout << "Replaced atom4_ with " << neighbor->getId() << "\n";
                    atom4_ = neighbor;
                    foundHydrogen = true;
                }
            }
            int numberOfNeighbors = atom3_->getNeighbors().size();
            if (!foundHydrogen && numberOfNeighbors <= 3)
            {
                //std::cout << "Creating HHH atom for " << atom3_->getId() << " as it has " << numberOfNeighbors << " neighbors\n";
                atom4_ = RotatableDihedral::CreateHydrogenAtomForPsi(atom3_);
            }
        }
    }
    return;
}

template<class atomT>
atomT* RotatableDihedral<atomT>::CreateHydrogenAtomForPsi(atomT* centralAtom)
{
    if(centralAtom->getNeighbors().size() != 4)
    {
        std::stringstream ss;
        ss << "Error in CreateHydrogenAtomForPsi. centralAtom neighbors size = " <<
                centralAtom->getNeighbors().size() << " for " << centralAtom->getId();
        gmml::log(__LINE__,__FILE__,gmml::ERR, ss.str());
        throw std::runtime_error(ss.str());
    }
    std::vector<Coordinate*> threeNeighborCoords;
    for (auto &neighbor : centralAtom->getNeighbors())
    {
        threeNeighborCoords.push_back(neighbor->getCoordinate());
    }
    Coordinate newCoord = GeometryTopology::CreateMissingCoordinateForTetrahedralAtom(centralAtom->getCoordinate(), threeNeighborCoords);
    atomT* newAtom = new atomT(centralAtom->GetResidue(), "HHH", newCoord);
    centralAtom->GetResidue()->AddAtom(newAtom);
    centralAtom->GetNode()->AddNodeNeighbor(newAtom);
    // Have to create a new node for this new atom.
    MolecularModeling::AtomNode *leakyNode = new MolecularModeling::AtomNode();
    newAtom->SetNode(leakyNode); // A jest!
    newAtom->GetNode()->AddNodeNeighbor(centralAtom);
    return newAtom;
}

template<class atomT>
std::string RotatableDihedral<atomT>::Print() const
{
   std::stringstream ss;
   ss << atom1_->getName() << ", " << atom2_->getName() << ", " << atom3_->getName() << ", " << atom4_->getName() << ": " << this->CalculateDihedralAngle()  << ".\n";
   gmml::log(__LINE__, __FILE__, gmml::INF, ss.str());
   return ss.str();
}
}//namespace
#endif
