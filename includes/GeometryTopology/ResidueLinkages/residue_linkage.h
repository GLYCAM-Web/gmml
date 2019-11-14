#ifndef RESIDUE_LINKAGE_H
#define RESIDUE_LINKAGE_H
/*
 * This class figures out the rotatable bonds between two residues
 * Starts/ends at the CA atoms in proteins. Looks for cycles (as they aren't rotatable).
 * Stores each rotatable bond as a rotatable_dihedral object.
 */
//#include "gmml.hpp"
#include "../../MolecularModeling/atom.hpp"
#include "../../MolecularModeling/residue.hpp"
#include "../../MolecularModeling/Selections/selections.hpp"
#include "rotatable_dihedral.h"

using MolecularModeling::Residue;
using MolecularModeling::ResidueVector;


class Residue_linkage;
typedef std::vector<Residue_linkage> ResidueLinkageVector;
class Residue_linkage
{
public:
    //////////////////////////////////////////////////////////
    //                    TYPE DEFINITION                   //
    //////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////
    //                       CONSTRUCTOR                    //
    //////////////////////////////////////////////////////////

    Residue_linkage();
    Residue_linkage(Residue *residue1, Residue *residue2);

    //////////////////////////////////////////////////////////
    //                       ACCESSOR                       //
    //////////////////////////////////////////////////////////

    ResidueVector GetResidues();
    RotatableDihedralVector GetRotatableDihedrals() const;
    int GetNumberOfRotatableDihedrals();

    //////////////////////////////////////////////////////////
    //                       MUTATOR                        //
    //////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////
    //                       FUNCTIONS                      //
    //////////////////////////////////////////////////////////

    void SetDefaultDihedralAnglesUsingMetadata();
    void SetRandomDihedralAnglesUsingMetadata();
    void SetCustomDihedralAngles(std::vector <double> dihedral_angles);
    void SetDihedralAnglesToPrevious();
    void SetRandomDihedralAngles();
    void DetermineAtomsThatMove();
    // Simple meaning you only check each rotatable_dihedral in series, not every combination.
    void SimpleWiggle(AtomVector overlapAtomSet1, AtomVector overlapAtomSet2, double overlapTolerance = 0.01, int interval = 5);

    //////////////////////////////////////////////////////////
    //                       DISPLAY FUNCTION               //
    //////////////////////////////////////////////////////////

    void Print();

    //////////////////////////////////////////////////////////
    //                       OPERATORS                      //
    //////////////////////////////////////////////////////////

private:

    //////////////////////////////////////////////////////////
    //                    PRIVATE FUNCTIONS                 //
    //////////////////////////////////////////////////////////

    void InitializeClass(Residue *from_this_residue1, Residue *to_this_residue2);
    bool CheckIfViableLinkage();
    RotatableDihedralVector FindRotatableDihedralsConnectingResidues(Atom *from_this_connection_atom1, Atom *to_this_connection_atom2);
    // Previous function generates a list of linearly connected atoms that define the rotatable bonds
    // This function splits that list into groups of 4 and creates rotatable_dihedral objects
    RotatableDihedralVector SplitAtomVectorIntoRotatableDihedrals(AtomVector atoms);
    gmml::MolecularMetadata::GLYCAM::DihedralAngleDataVector FindMetadata(Atom *from_this_connection_atom1, Atom *to_this_connection_atom2);
    void AddMetadataToRotatableDihedrals(gmml::MolecularMetadata::GLYCAM::DihedralAngleDataVector metadata);
    void SetResidues(Residue *residue1, Residue *residue2);
    void SetConnectionAtoms(Residue *residue1, Residue *residue2);

    //////////////////////////////////////////////////////////
    //                       ATTRIBUTES                     //
    //////////////////////////////////////////////////////////

    Residue* from_this_residue1_;
    Residue* to_this_residue2_;
    Atom* from_this_connection_atom1_;
    Atom* to_this_connection_atom2_;
    RotatableDihedralVector rotatable_dihedrals_;
    //gmml::MolecularMetadata::GLYCAM::DihedralAngleDataVector metadata_;
};

std::ostream& operator<<(std::ostream& os, const Residue_linkage&);

#endif // RESIDUE_LINKAGE_H
