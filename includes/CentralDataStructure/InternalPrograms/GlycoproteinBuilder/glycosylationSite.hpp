#ifndef GMML_INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_GLYCOSYLATIONSITE_HPP
#define GMML_INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_GLYCOSYLATIONSITE_HPP

#include "includes/CentralDataStructure/CondensedSequence/carbohydrate.hpp"
#include "includes/CentralDataStructure/assembly.hpp"
#include "includes/CentralDataStructure/Shapers/residueLinkage.hpp"

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
using cds::Atom;
using cds::Residue;
using cds::ResidueLinkage;
using cds::RotatableDihedral;
//using cds::Assembly;
using cdsCondensedSequence::Carbohydrate;

class GlycosylationSite
{
public:
    //////////////////////////////////////////////////////////
    //                       CONSTRUCTOR                    //
    //////////////////////////////////////////////////////////
//    GlycosylationSite(Assembly* glycoprotein, std::string residueNumber, std::string glycanInputType, std::string glycanName, std::string prepFileLocation);
    GlycosylationSite(Residue* residue, std::vector<Residue*> otherProteinResidues, std::string glycanInputString, std::string prepFileLocation);
   // ToDo: want this or with && and move the carb: GlycosylationSite(Residue* residue, std::vector<Residue*> otherProteinResidues, Carbohydrate carbohydrate);
    //////////////////////////////////////////////////////////
    //                       ACCESSOR                       //
    //////////////////////////////////////////////////////////
    inline double GetGlycanOverlap() 							{return glycan_overlap_;}
    inline double GetProteinOverlap() 							{return protein_overlap_;}
    inline std::string GetResidueId()							{return this->GetResidue()->getId();}
    inline double GetOverlap()									{return (glycan_overlap_ + protein_overlap_);}
//    std::vector<ResidueLinkage> GetFirstAnd1_6Linkages(); // This should be a selection in a separate place.
//    std::vector<ResidueLinkage> GetFirstAnd2_XLinkages(); // This should be a selection in a separate place.
    //////////////////////////////////////////////////////////
    //                       MUTATOR                        //
    //////////////////////////////////////////////////////////
    inline void SetProteinAtoms(std::vector<Atom*> proteinAtoms) {proteinAtoms_ = proteinAtoms;}
    void SetOtherGlycosites(std::vector<GlycosylationSite> &glycosites);
    inline std::vector<Atom*> GetSelfGlycanBeads()                      {return self_glycan_beads_;}
    inline void SetOtherGlycanBeads(std::vector<Atom*> *beads)          {other_glycan_beads_ = *beads;}
    //////////////////////////////////////////////////////////
    //                       FUNCTIONS                      //
    //////////////////////////////////////////////////////////
    void AddBeads(std::vector<Atom*> proteinBeads);
    double CalculateOverlaps(OverlapType overlapType = BEAD, MoleculeType moleculeType = ALL, bool recordOverlap = true, bool printOverlap = false);
    void UpdateAtomsThatMoveInLinkages();
//    void StashCoordinates();
//    void SetStashedCoordinates();
    void Wiggle(OverlapType overlapType = BEAD, bool firstLinkageOnly = false, double tolerance = 0.1, int interval = 5);
    //std::vector<GlycosylationSite> GetXClosestSitesWithinOverlapDistanceY(std::vector<GlycosylationSite> &glycosites, int maxNumberOfSitesToConsider);
    void SetRandomDihedralAnglesUsingMetadata();
    void SetRandomDihedralAnglesUsingMetadataForNthLinkage(long unsigned int linkage_number);
    void SetDefaultDihedralAnglesUsingMetadata();
    void ResetDihedralAngles();
    void Rename_Protein_Residue_From_GLYCAM_To_Standard();
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
        return rhs.residue_->getId() == residue_->getId();
    }
    inline bool operator!=(const GlycosylationSite &rhs) const
    {
        return residue_->getId() == rhs.residue_->getId();
    }
    inline bool operator<(const GlycosylationSite &rhs) const
    {
        return residue_->getId() < rhs.residue_->getId();
    }
    inline bool operator>(const GlycosylationSite &rhs) const
    {
        return residue_->getId() > rhs.residue_->getId();
    }
private:
    //////////////////////////////////////////////////////////
    //                  PRIVATE ACCESSOR                    //
    //////////////////////////////////////////////////////////
    inline std::vector<cds::ResidueLinkage> GetRotatableBonds() 			{return all_residue_linkages_;}
    inline std::vector<GlycosylationSite*> GetOtherGlycosites() {return other_glycosites_;}
    inline unsigned long int GetInternalBondCount() 							{return internalBondCount_;}
    inline std::vector<Atom*> GetProteinBeads() 						{return protein_beads_;}
    inline std::vector<Atom*> GetProteinAtoms()                         {return proteinAtoms_;}
    inline std::vector<Atom*> GetOtherGlycanBeads() 					{return other_glycan_beads_;}
    inline Carbohydrate* GetAttachedGlycan()			 			{return &glycan_;}
    inline  Residue* GetResidue() 								{return residue_;}
    double GetWeightedOverlap(double glycan_weight, double protein_weight);
    //////////////////////////////////////////////////////////
    //                  PRIVATE FRIENDS                     //
    //////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////
    //                  PRIVATE MUTATOR                     //
    //////////////////////////////////////////////////////////
    inline void SetGlycanOverlap(double d)						{glycan_overlap_ = d;}
    inline void SetProteinOverlap(double d)						{protein_overlap_ = d;}
    inline void SetSelfGlycanBeads(std::vector<Atom*> *beads)			{self_glycan_beads_ = *beads;}
    inline void SetInternalBondCount(int i) {internalBondCount_ = i;}
    void SetOverlap(MoleculeType moleculeType, double overlap);
    void SetProteinBeads(std::vector<Atom*> *beads);
    //////////////////////////////////////////////////////////
    //               INITIALIZATION FUNCTIONS               //
    //////////////////////////////////////////////////////////
    void FindSetProteinResidue(std::string residueName, std::vector<Residue*> residues);
    void BuildGlycan(std::string glycanInputType, std::string glycanName, std::string prepFileLocation);
    void AttachGlycan();
    //////////////////////////////////////////////////////////
    //                  PRIVATE FUNCTIONS                   //
    //////////////////////////////////////////////////////////
    void Prepare_Glycans_For_Superimposition_To_Particular_Residue(std::string amino_acid_name);
    void Superimpose_Glycan_To_Glycosite(Residue *glycosite_residue);
    void Rename_Protein_Residue_To_GLYCAM_Nomenclature();
    Atom* GetConnectingProteinAtom(std::string residue_name);
    void WiggleOneLinkage(ResidueLinkage &linkage, OverlapType overlapType = BEAD, double tolerance = 0.1, int interval = 5);
    double CalculateBeadOverlaps(MoleculeType moleculeType = ALL);
    //double Calculate_and_print_bead_overlaps();
    double CalculateAtomicOverlaps(MoleculeType moleculeType = ALL, bool print = false);
    double CalculateBeadOverlaps(std::vector<Atom*> &atomsA, std::vector<Atom*> &atomsB);
    bool NoNewInternalCloseContacts();
    //////////////////////////////////////////////////////////
    //                       ATTRIBUTES                     //
    //////////////////////////////////////////////////////////

    Residue* residue_;                                  /*!< A pointer back to the residue for this glycosite >*/
    Carbohydrate glycan_;
    double glycan_overlap_ = 0.0;
    double protein_overlap_ = 0.0;
    std::vector<cds::ResidueLinkage> all_residue_linkages_;
    std::vector<Atom*> self_glycan_beads_;
    std::vector<Atom*> other_glycan_beads_;
    std::vector<Atom*> protein_beads_;
    std::vector<Residue*> otherProteinResidues_;
    std::vector<GlycosylationSite*> other_glycosites_;
    std::vector<Atom*> proteinAtoms_;
    unsigned long int internalBondCount_; // For checking for internal clashes within the glycan
};
#endif // GMML_INCLUDES_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_GLYCOSYLATIONSITE_HPP
