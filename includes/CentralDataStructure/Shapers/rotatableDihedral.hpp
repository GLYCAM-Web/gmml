#ifndef GMML_INCLUDES_GEOMETRYTOPOLOGY_RESIDUELINKAGES_ROTATABLE_DIHEDRAL_HPP
#define GMML_INCLUDES_GEOMETRYTOPOLOGY_RESIDUELINKAGES_ROTATABLE_DIHEDRAL_HPP
/*
 * This class stores the four atoms that define a dihedral angle, the atoms that move when it is rotated
 * and, if moved, the previous dihedral angle, which allows me to reset easily.
 */

#include "includes/CentralDataStructure/atom.hpp"
#include "../../MolecularModeling/residue.hpp"
#include "../../MolecularMetadata/GLYCAM/dihedralangledata.hpp"

using cds::Atom;
using gmml::MolecularMetadata::GLYCAM::DihedralAngleData;
using gmml::MolecularMetadata::GLYCAM::DihedralAngleDataVector;

class Rotatable_dihedral
{
public:
    //////////////////////////////////////////////////////////
    //                       CONSTRUCTOR                    //
    //////////////////////////////////////////////////////////
    Rotatable_dihedral(Atom *atom1, Atom *atom2, Atom *atom3, Atom *atom4, bool reverseAtomsThatMove = true);
    Rotatable_dihedral(Atom *atom1, Atom *atom2, Atom *atom3, Atom *atom4, AtomVector extraAtomsThatMove, bool reverseAtomsThatMove = true);
    //////////////////////////////////////////////////////////
    //                       ACCESSOR                       //
    //////////////////////////////////////////////////////////
    const DihedralAngleDataVector& GetMetadata() const;
    DihedralAngleDataVector GetLikelyMetadata() const;
    int GetNumberOfRotamers( bool likelyShapesOnly = false) const;
    std::string GetName() const;
    double CalculateDihedralAngle(const std::string type = "default") const;
    //////////////////////////////////////////////////////////
    //                       MUTATOR                        //
    //////////////////////////////////////////////////////////
    // Based on connectivities, this figures out which atoms will move when the dihedral is rotated.
    void DetermineAtomsThatMove();
    // Sets the dihedral angle by rotating the bond between atom2 and atom3, moving atom4 and connected.
    void SetDihedralAngle(double dihedral_angle);
    // Sets the dihedral to previous dihedral angle
    void SetDihedralAngleToPrevious();
    // Randomly sets dihedral angle values between 0 and 360
    double RandomizeDihedralAngle();
    void AddMetadata(DihedralAngleData metadata);
    void ClearMetadata();
    void SetRandomAngleEntryUsingMetadata(bool useRanges = true);
    void SetSpecificAngleEntryUsingMetadata(bool useRanges = false, int angleEntryNumber = 0);
    bool SetSpecificShape(std::string dihedralName, std::string selectedRotamer);
    void WiggleWithinCurrentRotamer(std::vector<Atom*> &overlapAtomSet1, std::vector<Atom*> &overlapAtomSet2, const int &angleIncrement);
    void WiggleUsingAllRotamers(std::vector<Atom*>& overlapAtomSet1, std::vector<Atom*> &overlapAtomSet2, const int &angleIncrement);
    //////////////////////////////////////////////////////////
    //                       DISPLAY FUNCTION               //
    //////////////////////////////////////////////////////////
    std::string Print() const;
private:
    //////////////////////////////////////////////////////////
    //                  PRIVATE ACCESSORS                   //
    //////////////////////////////////////////////////////////
    std::vector<Atom*> GetAtoms() const;
    std::vector<Atom*> GetAtomsThatMove() const;
    bool GetIsAtomsThatMoveReversed() const;
    double GetPreviousDihedralAngle() const;
    std::vector<double> GetAllPossibleAngleValues(const int interval = 5) const;
    //////////////////////////////////////////////////////////
    //                  PRIVATE MUTATORS                    //
    //////////////////////////////////////////////////////////
    void AddExtraAtomsThatMove(std::vector<Atom*> extraAtoms);
    // Takes in a set of ranges, e.g. 10 to 30, 45-55 etc. Randomly selects a range and randomly sets value within that range.
    double RandomizeDihedralAngleWithinRanges(std::vector<std::pair<double,double> > ranges);
    // Randomly sets dihedral angle to a value within the given range. E.g. Between 25 and 30 degrees.
    void SetMetadata(DihedralAngleDataVector metadataVector);
    double RandomizeDihedralAngleWithinRange(double min, double max);
    //////////////////////////////////////////////////////////
    //                  PRIVATE FUNCTIONS                   //
    //////////////////////////////////////////////////////////
    void Initialize(std::vector<Atom*> atoms, bool reverseAtomsThatMove = true);
    void SetAtoms(std::vector<Atom*> atoms);
    void SetAtomsThatMove(std::vector<Atom*> atoms);
    void SetIsAtomsThatMoveReversed(bool isAtomsThatMoveReversed);
    void RecordPreviousDihedralAngle(double dihedral_angle);
    void UpdateAtomsIfPsi();
    Atom* CreateHydrogenAtomForPsi(Atom *centralAtom);
    void SetWasEverRotated(bool wasEverRotated);
    bool CheckIfEverRotated() const;
    inline void SetCurrentMetaData(const DihedralAngleData &d) {currentMetadata_ = &d;}
    inline const DihedralAngleData* GetCurrentMetaData() {return currentMetadata_;}
    double WiggleWithinRanges(std::vector<Atom*>& overlapAtomSet1, std::vector<Atom*> &overlapAtomSet2, const int &angleIncrement, const double& lowerBound, const double& upperBound);
    //////////////////////////////////////////////////////////
    //                       ATTRIBUTES                     //
    //////////////////////////////////////////////////////////
    // The four atoms that define the dihedral angle. The bond between atom2_ and atom3_ is what is rotated.
    Atom *atom1_;
    Atom *atom2_;
    Atom *atom3_;
    Atom *atom4_;
    // A vector of pointers to the atoms that are connected to atom2_ and atom3_, and will be rotated when that bond is rotated.
    std::vector<Atom*> atoms_that_move_;
    std::vector<Atom*> extra_atoms_that_move_;
    bool isAtomsThatMoveReversed_;
    // I often want to reset a dihedral angle after rotating it, so recording the previous angle makes this easy.
    double previous_dihedral_angle_;
    DihedralAngleDataVector assigned_metadata_;
    const DihedralAngleData* currentMetadata_;
    bool wasEverRotated_; // Need this, as it might add a H atom for psi
};
#endif // ROTATABLE_DIHEDRAL_H
