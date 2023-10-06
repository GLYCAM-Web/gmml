#ifndef GMML_INCLUDES_CENTRAL_DATA_STRUCTURE_SHAPERS_ROTATABLE_DIHEDRAL_HPP
#define GMML_INCLUDES_CENTRAL_DATA_STRUCTURE_SHAPERS_ROTATABLE_DIHEDRAL_HPP

#include "includes/CentralDataStructure/atom.hpp"
#include "includes/CentralDataStructure/coordinate.hpp"
#include "includes/MolecularMetadata/GLYCAM/dihedralangledata.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/External_Libraries/PCG/pcg_random.h"
#include "includes/CentralDataStructure/Overlaps/cdsOverlaps.hpp"

#include <random>

using gmml::MolecularMetadata::GLYCAM::DihedralAngleData;
using gmml::MolecularMetadata::GLYCAM::DihedralAngleDataVector;
// Seed with a real random value, if available
static pcg_extras::seed_seq_from<std::random_device> seed_source;
// Make a random number engine
static pcg32
    rng(seed_source); // Eclipse complains about ambiguity, and yet it compiles... // @suppress("Ambiguous problem")

// This class stores the four atoms that define a dihedral angle, the atoms that move when it is rotated
// and, if moved, the previous dihedral angle, which allows me to reset easily.
namespace cds
{
    class RotatableDihedral
    {
      public:
        //////////////////////////////////////////////////////////
        //                       CONSTRUCTOR                    //
        //////////////////////////////////////////////////////////
        RotatableDihedral(cds::Atom* atom1, cds::Atom* atom2, cds::Atom* atom3, cds::Atom* atom4,
                          bool reverseAtomsThatMove = true);
        RotatableDihedral(cds::Atom* atom1, cds::Atom* atom2, cds::Atom* atom3, cds::Atom* atom4,
                          std::vector<cds::Atom*> extraAtomsThatMove, bool reverseAtomsThatMove = true);

        //////////////////////////////////////////////////////////
        //                       ACCESSOR                       //
        //////////////////////////////////////////////////////////
        inline const DihedralAngleDataVector& GetMetadata() const
        {
            return assigned_metadata_;
        }

        DihedralAngleDataVector GetLikelyMetadata() const;
        int GetNumberOfRotamers(bool likelyShapesOnly = false) const;
        std::string GetName() const;
        double CalculateDihedralAngle(const std::string type = "default") const;
        //////////////////////////////////////////////////////////
        //                       MUTATOR                        //
        //////////////////////////////////////////////////////////
        void DetermineAtomsThatMove(); // Based on connectivities, this figures out which atoms will move when the
                                       // dihedral is rotated.
        void SetDihedralAngle(double dihedral_angle); // Sets the dihedral angle by rotating the bond between atom2 and
                                                      // atom3, moving atom4 and connected.
        void SetDihedralAngleToPrevious();            // Sets the dihedral to previous dihedral angle
        double RandomizeDihedralAngle();              // Randomly sets dihedral angle values between 0 and 360
        void AddMetadata(DihedralAngleData metadata);

        inline void ClearMetadata()
        {
            assigned_metadata_.clear();
        }

        void SetRandomAngleEntryUsingMetadata(bool useRanges = true);
        void SetSpecificAngleEntryUsingMetadata(bool useRanges = false, long unsigned int angleEntryNumber = 0);
        bool SetSpecificShape(std::string dihedralName, std::string selectedRotamer);
        void WiggleWithinCurrentRotamer(std::vector<cds::Atom*>& overlapAtomSet1,
                                        std::vector<cds::Atom*>& overlapAtomSet2, const int& angleIncrement);
        void WiggleWithinCurrentRotamer(std::vector<cds::Residue*>& overlapResidueSet1,
                                        std::vector<cds::Residue*>& overlapResidueSet2, const int& angleIncrement);
        void WiggleUsingAllRotamers(std::vector<cds::Atom*>& overlapAtomSet1, std::vector<cds::Atom*>& overlapAtomSet2,
                                    const int& angleIncrement);
        bool IsThereHydrogenForPsiAngle(); // ToDo which of these functions are only public for ResidueLinkage? How
                                           // about friend function huh?
        std::unique_ptr<cds::Atom> CreateHydrogenAtomForPsiAngle();
        //////////////////////////////////////////////////////////
        //                       DISPLAY FUNCTION               //
        //////////////////////////////////////////////////////////
        std::string Print() const;

      private:
        //////////////////////////////////////////////////////////
        //                  PRIVATE ACCESSORS                   //
        //////////////////////////////////////////////////////////
        std::vector<cds::Atom*> GetAtoms() const;

        inline bool GetIsAtomsThatMoveReversed() const
        {
            return isAtomsThatMoveReversed_;
        }

        inline double GetPreviousDihedralAngle() const
        {
            return previous_dihedral_angle_;
        }

        std::vector<double> GetAllPossibleAngleValues(const int interval = 5) const;
        //////////////////////////////////////////////////////////
        //                  PRIVATE MUTATORS                    //
        //////////////////////////////////////////////////////////
        void AddExtraAtomsThatMove(std::vector<cds::Atom*>& extraAtoms);
        void SetAtomsThatMove(std::vector<cds::Atom*>& atoms);
        double
        RandomizeDihedralAngleWithinRanges(std::vector<std::pair<double, double>>
                                               ranges); // Takes in a set of ranges, e.g. 10 to 30, 45-55 etc. Randomly
                                                        // selects a range and randomly sets value within that range.
        double RandomizeDihedralAngleWithinRange(double min,
                                                 double max); // Randomly sets dihedral angle to a value within the
                                                              // given range. E.g. Between 25 and 30 degrees.
        //////////////////////////////////////////////////////////
        //                  PRIVATE FUNCTIONS                   //
        //////////////////////////////////////////////////////////
        void Initialize(std::vector<cds::Atom*> atoms, bool reverseAtomsThatMove = true);
        void SetAtoms(std::vector<cds::Atom*> atoms);

        inline void SetIsAtomsThatMoveReversed(bool b)
        {
            isAtomsThatMoveReversed_ = b;
        }

        inline void RecordPreviousDihedralAngle(double d)
        {
            previous_dihedral_angle_ = d;
        }

        inline void SetWasEverRotated(bool b)
        {
            wasEverRotated_ = b;
        }

        inline bool CheckIfEverRotated() const
        {
            return wasEverRotated_;
        }

        inline void SetCurrentMetaData(const DihedralAngleData& d)
        {
            currentMetadata_ = &d;
        }

        inline const DihedralAngleData* GetCurrentMetaData()
        {
            return currentMetadata_;
        }

        // double WiggleWithinRanges(std::vector<cds::Atom*>& overlapAtomSet1, std::vector<cds::Atom*>& overlapAtomSet2,
        //                                  const int& angleIncrement, const double& lowerBound, const double&
        //                                  upperBound);
        unsigned int WiggleWithinRangesDistanceCheck(std::vector<cds::Atom*>& overlapAtomSet1,
                                                     std::vector<cds::Atom*>& overlapAtomSet2,
                                                     const int& angleIncrement, const double& lowerBound,
                                                     const double& upperBound);
        unsigned int WiggleWithinRangesDistanceCheck(std::vector<cds::Residue*>& overlapResidueSet1,
                                                     std::vector<cds::Residue*>& overlapResidueSet2,
                                                     const int& angleIncrement, const double& lowerBound,
                                                     const double& upperBound);

        inline std::vector<cds::Coordinate*>& GetCoordinatesThatMove()
        {
            return coordinatesThatMove_;
        }

        //////////////////////////////////////////////////////////
        //                       ATTRIBUTES                     //
        //////////////////////////////////////////////////////////
        // The four atoms that define the dihedral angle. The bond between atom2_ and atom3_ is what is rotated.
        cds::Atom* atom1_;
        cds::Atom* atom2_;
        cds::Atom* atom3_;
        cds::Atom* atom4_;
        // A vector of pointers to the atoms that are connected to atom2_ and atom3_, and will be rotated when that bond
        // is rotated.
        //    std::vector<cds::Atom*> atoms_that_move_;
        //    std::vector<cds::Atom*> extra_atoms_that_move_;
        std::vector<cds::Coordinate*> coordinatesThatMove_;
        bool isAtomsThatMoveReversed_ = true;
        double previous_dihedral_angle_; // I often want to reset a dihedral angle after rotating it, so recording the
                                         // previous angle makes this easy.
        DihedralAngleDataVector assigned_metadata_;
        const DihedralAngleData* currentMetadata_;
        bool wasEverRotated_ = false; // Need this, as it might add a H atom for psi
    };

} // namespace cds
#endif
