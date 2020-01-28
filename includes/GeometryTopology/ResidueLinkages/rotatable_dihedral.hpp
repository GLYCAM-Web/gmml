#ifndef ROTATABLE_DIHEDRAL_H
#define ROTATABLE_DIHEDRAL_H
/*
 * This class stores the four atoms that define a dihedral angle, the atoms that move when it is rotated
 * and, if moved, the previous dihedral angle, which allows me to reset easily.
 */

#include "../../MolecularModeling/atom.hpp"
#include "../../MolecularModeling/residue.hpp"
#include "../../MolecularMetadata/GLYCAM/dihedralangledata.hpp"

using MolecularModeling::Atom;
using MolecularModeling::AtomVector;
using gmml::MolecularMetadata::GLYCAM::DihedralAngleData;
using gmml::MolecularMetadata::GLYCAM::DihedralAngleDataVector;

class Rotatable_dihedral;
typedef std::vector<Rotatable_dihedral> RotatableDihedralVector;
class Rotatable_dihedral
{
public:
    //////////////////////////////////////////////////////////
    //                    TYPE DEFINITION                   //
    //////////////////////////////////////////////////////////


    //////////////////////////////////////////////////////////
    //                       CONSTRUCTOR                    //
    //////////////////////////////////////////////////////////

    Rotatable_dihedral(Atom *atom1, Atom *atom2, Atom *atom3, Atom *atom4, bool reverseAtomsThatMove = true);
    Rotatable_dihedral(Atom *atom1, Atom *atom2, Atom *atom3, Atom *atom4, AtomVector extraAtomsThatMove, bool reverseAtomsThatMove = true);

//    Rotatable_dihedral(AtomVector atoms, bool reverseAtomsThatMove = true);
//    Rotatable_dihedral(AtomVector atoms, AtomVector atoms_that_move);


    //////////////////////////////////////////////////////////
    //                       ACCESSOR                       //
    //////////////////////////////////////////////////////////

    double CalculateDihedralAngle() const;
    AtomVector GetAtoms() const;
    AtomVector GetAtomsThatMove();
    bool GetIsAtomsThatMoveReversed();
    double GetPreviousDihedralAngle();
    DihedralAngleDataVector GetMetadata();
    int GetNumberOfRotamers();
    std::vector<double> GetAllPossibleAngleValues(int interval = 5);

    //////////////////////////////////////////////////////////
    //                       FUNCTIONS                      //
    //////////////////////////////////////////////////////////


    //////////////////////////////////////////////////////////
    //                       MUTATOR                        //
    //////////////////////////////////////////////////////////

    // Based on connectivities, this figures out which atoms will move when the dihedral is rotated.
    void DetermineAtomsThatMove();
    void AddExtraAtomsThatMove(AtomVector extraAtoms);
    // Sets the dihedral angle by rotating the bond between atom2 and atom3, moving atom4 and connected.
    void SetDihedralAngle(double dihedral_angle);
    // Sets the dihedral to previous dihedral angle
    void SetDihedralAngleToPrevious();
    // Randomly sets dihedral angle values between 0 and 360
    double RandomizeDihedralAngle();
    // Takes in a set of ranges, e.g. 10 to 30, 45-55 etc. Randomly selects a range and randomly sets value within that range.
    double RandomizeDihedralAngleWithinRanges(std::vector<std::pair<double,double> > ranges);
    // Randomly sets dihedral angle to a value within the given range. E.g. Between 25 and 30 degrees.
    double RandomizeDihedralAngleWithinRange(double min, double max);

    // ALTER CONSTRUCTOR SO THESE  next two ARE PRIVATE?

    // A residue-residue linkage will have metadata for each rotatable_dihedral. Multiple rotamers means multiple entries.
    void SetMetadata(DihedralAngleDataVector metadataVector);
    void AddMetadata(DihedralAngleData metadata);

    void SetRandomAngleEntryUsingMetadata(bool useRanges = true);
    void SetSpecificAngleEntryUsingMetadata(bool useRanges = false, int angleEntryNumber = 0);


    //////////////////////////////////////////////////////////
    //                       DISPLAY FUNCTION               //
    //////////////////////////////////////////////////////////

    void Print();

    //////////////////////////////////////////////////////////
    //                       OPERATORS                      //
    //////////////////////////////////////////////////////////

private:

    //////////////////////////////////////////////////////////
    //                  PRIVATE FUNCTIONS                   //
    //////////////////////////////////////////////////////////

    void Initialize(AtomVector atoms, bool reverseAtomsThatMove = true);
    void SetAtoms(AtomVector atoms);
    void SetAtomsThatMove(AtomVector atoms);
    void SetIsAtomsThatMoveReversed(bool isAtomsThatMoveReversed);
    void RecordPreviousDihedralAngle(double dihedral_angle);
    void UpdateAtomsIfPsi();
    Atom* CreateHydrogenAtomForPsi(Atom *centralAtom);
    void SetWasEverRotated(bool wasEverRotated);
    bool CheckIfEverRotated();

    //////////////////////////////////////////////////////////
    //                       ATTRIBUTES                     //
    //////////////////////////////////////////////////////////

    // The four atoms that define the dihedral angle. The bond between atom2_ and atom3_ is what is rotated.
    Atom *atom1_;
    Atom *atom2_;
    Atom *atom3_;
    Atom *atom4_;
    // A vector of pointers to the atoms that are connected to atom2_ and atom3_, and will be rotated when that bond is rotated.
    AtomVector atoms_that_move_;
    AtomVector extra_atoms_that_move_;
    bool isAtomsThatMoveReversed_;
    // I often want to reset a dihedral angle after rotating it, so recording the previous angle makes this easy.
    double previous_dihedral_angle_;
    DihedralAngleDataVector assigned_metadata_;
    bool wasEverRotated_; // Need this, as it might add a H atom for psi

};

std::ostream& operator<<(std::ostream& os, const Rotatable_dihedral&);

#endif // ROTATABLE_DIHEDRAL_H
