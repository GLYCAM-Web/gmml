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
#include "rotatable_dihedral.hpp"



using MolecularModeling::Residue;
using MolecularModeling::ResidueVector;


class Residue_linkage;
typedef std::vector<Residue_linkage> ResidueLinkageVector;
class Residue_linkage
{
public:
    //////////////////////////////////////////////////////////
    //                       CONSTRUCTOR                    //
    //////////////////////////////////////////////////////////
    Residue_linkage(Residue *nonReducingResidue1, Residue *reducingResidue2, bool reverseAtomsThatMove = true);
    Residue_linkage(Residue *nonReducingResidue1, Residue *reducingResidue2, AtomVector alsoMovingAtoms, bool reverseAtomsThatMove = true);
    //////////////////////////////////////////////////////////
    //                       ACCESSOR                       //
    //////////////////////////////////////////////////////////
    ResidueVector GetResidues();
    bool GetIfReversedAtomsThatMove();
    std::vector<Rotatable_dihedral> GetRotatableDihedrals() const;
    std::vector<Rotatable_dihedral> GetRotatableDihedralsWithMultipleRotamers();
    int GetNumberOfRotatableDihedrals();
    int GetNumberOfShapes(bool likelyShapesOnly = false);
    Residue* GetFromThisResidue1();
    Residue* GetToThisResidue2();
    Atom* GetFromThisConnectionAtom1();
    Atom* GetToThisConnectionAtom2();
    bool CheckIfConformer();
    bool GetIfExtraAtoms();
    AtomVector GetExtraAtoms();
    void AddExtraAtomsThatMove(AtomVector extraAtoms);
    unsigned long long GetIndex();
    std::string GetName();
    //////////////////////////////////////////////////////////
    //                       MUTATOR                        //
    //////////////////////////////////////////////////////////
    void SetRotatableDihedrals(std::vector<Rotatable_dihedral> rotatableDihedrals);
    void SetIfReversedAtomsThatMove(bool reversedAtomsThatMove);
    //////////////////////////////////////////////////////////
    //                       FUNCTIONS                      //
    //////////////////////////////////////////////////////////
    void GenerateAllShapesUsingMetadata();
    void SetDefaultShapeUsingMetadata();
    void SetRandomShapeUsingMetadata(bool useRanges = true);
    void SetSpecificShapeUsingMetadata(int shapeNumber, bool useRanges = false);
    void SetSpecificShape(std::string dihedralName, std::string selectedRotamer);
    void SetCustomDihedralAngles(std::vector <double> dihedral_angles);
    void SetShapeToPrevious();
    void SetRandomDihedralAngles();
    void DetermineAtomsThatMove();
    // Simple meaning you only check each rotatable_dihedral in series, not every combination.
    void SimpleWiggle(AtomVector overlapAtomSet1, AtomVector overlapAtomSet2, double overlapTolerance = 0.01, int angleIncrement = 5);
    void SetIndex(unsigned long long index);
    void SetName(std::string name);
    //////////////////////////////////////////////////////////
    //                       DISPLAY FUNCTION               //
    //////////////////////////////////////////////////////////
    std::string Print();
private:
    //////////////////////////////////////////////////////////
    //                    PRIVATE FUNCTIONS                 //
    //////////////////////////////////////////////////////////
    void InitializeClass(Residue *from_this_residue1, Residue *to_this_residue2, bool reverseAtomsThatMove);
    bool CheckIfViableLinkage();
    std::vector<Rotatable_dihedral> FindRotatableDihedralsConnectingResidues(Atom *from_this_connection_atom1, Atom *to_this_connection_atom2);
    //AtomVector DealWithBranchesFromLinkages(AtomVector linearLinkageAtoms, Atom *cycle_point1, Atom *cycle_point2);
    // Previous function generates a list of linearly connected atoms that define the rotatable bonds
    // This function splits that list into groups of 4 and creates rotatable_dihedral objects
    std::vector<Rotatable_dihedral> SplitAtomVectorIntoRotatableDihedrals(AtomVector atoms);
    gmml::MolecularMetadata::GLYCAM::DihedralAngleDataVector FindMetadata(Atom *from_this_connection_atom1, Atom *to_this_connection_atom2);
    void AddMetadataToRotatableDihedrals(gmml::MolecularMetadata::GLYCAM::DihedralAngleDataVector metadata);
    void SetResidues(Residue *residue1, Residue *residue2);
    void SetConnectionAtoms(Residue *residue1, Residue *residue2);
    void SetConformerUsingMetadata(bool useRanges = false, int conformerNumber = 0);
    unsigned long long GenerateIndex();
    std::string DetermineLinkageNameFromResidueNames();
    //////////////////////////////////////////////////////////
    //                       ATTRIBUTES                     //
    //////////////////////////////////////////////////////////
    Residue* from_this_residue1_;
    Residue* to_this_residue2_;
    Atom* from_this_connection_atom1_;
    Atom* to_this_connection_atom2_;
    std::vector<Rotatable_dihedral> rotatable_dihedrals_;
    bool reverseAtomsThatMove_;
    AtomVector extraAtomsThatMove_;
    bool isExtraAtoms_ = true;
    unsigned long long index_;
    std::string name_; //e.g. "DGalpb1-6DGlcpNAc"
};
#endif // RESIDUE_LINKAGE_H
