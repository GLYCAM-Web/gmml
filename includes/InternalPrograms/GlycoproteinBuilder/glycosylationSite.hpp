#ifndef GMML_INCLUDES_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_GLYCOSYLATIONSITE_HPP
#define GMML_INCLUDES_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_GLYCOSYLATIONSITE_HPP
#include <iomanip> // For setting precision and formating in std::cout 
#include <algorithm> //  std::erase, std::remove
#include "includes/gmml.hpp"

enum OverlapType
{
    BEAD,
    ATOMIC,
};
enum MoleculeType
{
    PROTEIN,
    GLYCAN,
    ALL
};
using namespace MolecularModeling;
class GlycosylationSite
{
public:
    //////////////////////////////////////////////////////////
    //                       CONSTRUCTOR                    //
    //////////////////////////////////////////////////////////
    GlycosylationSite(Assembly* glycoprotein, std::string residueNumber, std::string glycanInputType, std::string glycanName, std::string prepFileLocation);
    //////////////////////////////////////////////////////////
    //                       ACCESSOR                       //
    //////////////////////////////////////////////////////////
    inline double GetGlycanOverlap() 							{return glycan_overlap_;}
    inline double GetProteinOverlap() 							{return protein_overlap_;}
    inline std::string GetResidueNumber() 						{return residue_number_;}
    inline std::string GetResidueId()							{return this->GetResidue()->GetId();}
    inline double GetOverlap()									{return (glycan_overlap_ + protein_overlap_);}
//    ResidueLinkageVector GetFirstAnd1_6Linkages(); // This should be a selection in a separate place.
//    ResidueLinkageVector GetFirstAnd2_XLinkages(); // This should be a selection in a separate place.
    //////////////////////////////////////////////////////////
    //                       FUNCTIONS                      //
    //////////////////////////////////////////////////////////
    void AddBeads(AtomVector proteinBeads);
    void Remove(Assembly *glycoproteinAssembly);
    double CalculateOverlaps(OverlapType overlapType = BEAD, MoleculeType moleculeType = ALL, bool recordOverlap = true, bool printOverlap = false);
    void UpdateAtomsThatMoveInLinkages();
//    void StashCoordinates();
//    void SetStashedCoordinates();
    void Wiggle(OverlapType overlapType = BEAD, bool firstLinkageOnly = false, double tolerance = 0.1, int interval = 5);
    //std::vector<GlycosylationSite> GetXClosestSitesWithinOverlapDistanceY(std::vector<GlycosylationSite> &glycosites, int maxNumberOfSitesToConsider);
    void SetRandomDihedralAnglesUsingMetadata();
    void SetRandomDihedralAnglesUsingMetadataForNthLinkage(int linkage_number);
    void SetDefaultDihedralAnglesUsingMetadata();
    void ResetDihedralAngles();
    void SetOtherGlycosites(std::vector<GlycosylationSite> &glycosites);
    inline AtomVector GetSelfGlycanBeads() 						{return self_glycan_beads_;}
    inline void SetOtherGlycanBeads(AtomVector *beads)			{other_glycan_beads_ = *beads;}
    //////////////////////////////////////////////////////////
    //                       DISPLAY FUNCTION               //
    //////////////////////////////////////////////////////////
    void PrintOverlaps();
    void Print(std::string type = "All");
    //////////////////////////////////////////////////////////
    //                       OPERATORS                      //
    //////////////////////////////////////////////////////////
    inline bool operator==(const GlycosylationSite &rhs) const
    {
        return rhs.residue_->GetId() == residue_->GetId();
    }
    inline bool operator!=(const GlycosylationSite &rhs) const
    {
        return residue_->GetId() == rhs.residue_->GetId();
    }
    inline bool operator<(const GlycosylationSite &rhs) const
    {
        return residue_->GetId() < rhs.residue_->GetId();
    }
    inline bool operator>(const GlycosylationSite &rhs) const
    {
        return residue_->GetId() > rhs.residue_->GetId();
    }
private:
    //////////////////////////////////////////////////////////
    //                  PRIVATE ACCESSOR                    //
    //////////////////////////////////////////////////////////
    inline std::string GetGlycanName() 							{return glycan_name_;}
    inline ResidueLinkageVector GetRotatableBonds() 			{return all_residue_linkages_;}
    inline std::vector<GlycosylationSite*> GetOtherGlycosites() {return other_glycosites_;}
    inline int GetInternalBondCount() 							{return internalBondCount_;}
    inline AtomVector GetProteinBeads() 						{return protein_beads_;}
    inline AtomVector GetOtherGlycanBeads() 					{return other_glycan_beads_;}
    inline Assembly* GetAttachedGlycan()			 			{return &glycan_;}
    inline  Residue* GetResidue() 								{return residue_;}
    inline std::vector<Atom*> GetProteinAtoms() 				{return proteinAtoms_;}
    double GetWeightedOverlap(double glycan_weight, double protein_weight);
    //////////////////////////////////////////////////////////
    //                  PRIVATE FRIENDS                     //
    //////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////
    //                  PRIVATE MUTATOR                     //
    //////////////////////////////////////////////////////////
    inline void SetProteinAtoms(std::vector<Atom*> atoms) 		{proteinAtoms_ = atoms;}
    inline void SetGlycanName(std::string s) 					{glycan_name_ = s;}
    inline void SetResidueNumber(std::string s) 				{residue_number_ = s;}
    inline void SetResidue(Residue* r)							{residue_ = r;}
    inline void SetGlycan(Assembly g)							{glycan_ = g;}
    inline void SetGlycanOverlap(double d)						{glycan_overlap_ = d;}
    inline void SetProteinOverlap(double d)						{protein_overlap_ = d;}
    inline void SetSelfGlycanBeads(AtomVector *beads)			{self_glycan_beads_ = *beads;}
    void Rename_Protein_Residue_From_GLYCAM_To_Standard();
    inline void SetInternalBondCount(int i) {internalBondCount_ = i;}
    void SetOverlap(MoleculeType moleculeType, double overlap);
    void SetProteinBeads(AtomVector *beads);
    //////////////////////////////////////////////////////////
    //               INITIALIZATION FUNCTIONS               //
    //////////////////////////////////////////////////////////
    void FindSetProteinResidue(std::string residueName, ResidueVector residues);
    void DetermineOverlapProteinAtoms(Assembly* glycoprotein, Residue* attachmentResidue);
    void BuildGlycan(std::string glycanInputType, std::string glycanName, std::string prepFileLocation);
    void AttachGlycan(Assembly* glycoprotein);
    //////////////////////////////////////////////////////////
    //                  PRIVATE FUNCTIONS                   //
    //////////////////////////////////////////////////////////
    void Prepare_Glycans_For_Superimposition_To_Particular_Residue(std::string amino_acid_name);
    void Superimpose_Glycan_To_Glycosite(Residue *glycosite_residue);
    void Rename_Protein_Residue_To_GLYCAM_Nomenclature();
    void FigureOutResidueLinkagesInGlycan(Residue *from_this_residue1, Residue *to_this_residue2, ResidueLinkageVector *residue_linkages);
    void RecursivelyGetAllNeighboringResidues(Atom* current_atom, ResidueVector* neighbors);
    Atom* GetConnectingProteinAtom(std::string residue_name);
    void WiggleOneLinkage(Residue_linkage &linkage, OverlapType overlapType = BEAD, double tolerance = 0.1, int interval = 5);
    double CalculateBeadOverlaps(MoleculeType moleculeType = ALL);
    //double Calculate_and_print_bead_overlaps();
    double CalculateAtomicOverlaps(MoleculeType moleculeType = ALL, bool print = false);
    double CalculateBeadOverlaps(AtomVector &atomsA, AtomVector &atomsB);
    bool NoNewInternalCloseContacts();
    //////////////////////////////////////////////////////////
    //                       ATTRIBUTES                     //
    //////////////////////////////////////////////////////////
    std::string glycan_name_;
    std::string residue_number_; // Don't need to be saved.
    Residue* residue_;                                  /*!< A pointer back to the residue for this glycosite >*/
    Assembly glycan_;
    AtomVector superimposition_atoms_;  // Don't need to be saved.             /*!< The 3 atoms used for superimposition of glycan to sidechain >*/
    double glycan_overlap_;
    double protein_overlap_;
    ResidueLinkageVector all_residue_linkages_;
    AtomVector self_glycan_beads_;
    AtomVector other_glycan_beads_;
    AtomVector protein_beads_;
    std::vector<Atom*> proteinAtoms_;
    std::vector<GlycosylationSite*> other_glycosites_;
    int internalBondCount_; // For checking for internal clashes within the glycan
};
#endif // GMML_INCLUDES_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_GLYCOSYLATIONSITE_HPP
