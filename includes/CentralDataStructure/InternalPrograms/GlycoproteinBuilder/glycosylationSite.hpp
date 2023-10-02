#ifndef GMML_INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_GLYCOSYLATIONSITE_HPP
#define GMML_INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_GLYCOSYLATIONSITE_HPP

#include "includes/CentralDataStructure/CondensedSequence/carbohydrate.hpp"
#include "includes/CentralDataStructure/assembly.hpp"
#include "includes/CentralDataStructure/Shapers/residueLinkage.hpp"

enum Resolution
{
    RESIDUE,
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
// using cds::Assembly;
using cdsCondensedSequence::Carbohydrate;

class GlycosylationSite
{
  public:
    //////////////////////////////////////////////////////////
    //                       CONSTRUCTOR                    //
    //////////////////////////////////////////////////////////
    //    GlycosylationSite(Assembly* glycoprotein, std::string residueNumber, std::string glycanInputType, std::string
    //    glycanName, std::string prepFileLocation);
    GlycosylationSite(Residue* residue, Carbohydrate* carbohydrate, std::vector<Residue*> otherProteinResidues,
                      unsigned int glycanStartResidueNumber);

    // ToDo: want this or with && and move the carb: GlycosylationSite(Residue* residue, std::vector<Residue*>
    // otherProteinResidues, Carbohydrate carbohydrate);
    //////////////////////////////////////////////////////////
    //                       ACCESSOR                       //
    //////////////////////////////////////////////////////////
    //    inline double GetGlycanOverlap() 							{return
    //    glycan_overlap_;} inline double GetProteinOverlap() 							{return
    //    protein_overlap_;}
    inline std::string GetResidueId()
    {
        return this->GetResidue()->getStringId();
    }

    //    inline double GetOverlap()									{return
    //    (glycan_overlap_
    //    + protein_overlap_);} std::vector<ResidueLinkage> GetFirstAnd1_6Linkages(); // This should be a selection in a
    //    separate place. std::vector<ResidueLinkage> GetFirstAnd2_XLinkages(); // This should be a selection in a
    //    separate place.
    //////////////////////////////////////////////////////////
    //                       MUTATOR                        //
    //////////////////////////////////////////////////////////
    // inline void SetProteinAtoms(std::vector<Atom*> proteinAtoms) {proteinAtoms_ = proteinAtoms;}
    void SetOtherGlycosites(std::vector<GlycosylationSite>& glycosites);
    //////////////////////////////////////////////////////////
    //                       FUNCTIONS                      //
    //////////////////////////////////////////////////////////
    int CalculateOverlaps(Resolution overlapType = RESIDUE, MoleculeType moleculeType = ALL);
    //    void StashCoordinates();
    //    void SetStashedCoordinates();
    void Wiggle(bool firstLinkageOnly = false, int interval = 5);
    // std::vector<GlycosylationSite> GetXClosestSitesWithinOverlapDistanceY(std::vector<GlycosylationSite> &glycosites,
    // int maxNumberOfSitesToConsider);
    void SetRandomDihedralAnglesUsingMetadata();
    // void SetDefaultDihedralAnglesUsingMetadata();
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
    inline bool operator==(const GlycosylationSite& rhs) const
    {
        return rhs.residue_->getStringId() == residue_->getStringId();
    }

    inline bool operator!=(const GlycosylationSite& rhs) const
    {
        return residue_->getStringId() != rhs.residue_->getStringId();
    }

    inline bool operator<(const GlycosylationSite& rhs) const
    {
        return residue_->getStringId() < rhs.residue_->getStringId();
    }

    inline bool operator>(const GlycosylationSite& rhs) const
    {
        return residue_->getStringId() > rhs.residue_->getStringId();
    }

  private:
    //////////////////////////////////////////////////////////
    //                  PRIVATE ACCESSOR                    //
    //////////////////////////////////////////////////////////
    inline std::vector<GlycosylationSite*> GetOtherGlycosites()
    {
        return other_glycosites_;
    }

    inline std::vector<Residue*> GetOtherProteinResidues()
    {
        return otherProteinResidues_;
    }

    inline Carbohydrate* GetGlycan()
    {
        return glycan_;
    }

    inline Residue* GetResidue() const
    {
        return residue_;
    }

    inline ResidueLinkage& GetProteinGlycanLinkage()
    {
        return this->GetGlycan()->GetGlycosidicLinkages().front();
    }

    inline unsigned long int GetInternalBondCount()
    {
        return internalBondCount_;
    }

    //////////////////////////////////////////////////////////
    //                  PRIVATE FRIENDS                     //
    //////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////
    //                  PRIVATE MUTATOR                     //
    //////////////////////////////////////////////////////////
    //    inline void SetGlycanOverlap(double d)						{glycan_overlap_ = d;}
    //    inline void SetProteinOverlap(double d)						{protein_overlap_ = d;}
    // inline void SetSelfGlycanBeads(std::vector<Atom*> *beads)			{self_glycan_beads_ = *beads;}
    inline void SetInternalBondCount(unsigned long int i)
    {
        internalBondCount_ = i;
    }

    // void SetOverlap(MoleculeType moleculeType, double overlap);
    // void SetProteinBeads(std::vector<Atom*> beads);
    //////////////////////////////////////////////////////////
    //               INITIALIZATION FUNCTIONS               //
    //////////////////////////////////////////////////////////
    // void FindSetProteinResidue(std::string residueName, std::vector<Residue*> residues);
    // void BuildGlycan(std::string glycanInputType, std::string glycanName, std::string prepFileLocation);
    void AttachGlycan(unsigned int glycanResidueStartNumber);
    //////////////////////////////////////////////////////////
    //                  PRIVATE FUNCTIONS                   //
    //////////////////////////////////////////////////////////
    void RenumberGlycanToMatch(unsigned int startNumber);
    void Prepare_Glycans_For_Superimposition_To_Particular_Residue(std::string amino_acid_name);
    void Superimpose_Glycan_To_Glycosite(Residue* glycosite_residue);
    void Rename_Protein_Residue_To_GLYCAM_Nomenclature();
    Atom* GetConnectingProteinAtom(const std::string residue_name) const;
    void WiggleOneLinkage(ResidueLinkage& linkage, int interval = 5);
    // double Calculate_and_print_bead_overlaps();
    int CalculateOverlaps(Resolution resolutionLevel, const std::vector<Residue*> residuesA,
                          const std::vector<Residue*> residuesB);
    bool NoNewInternalCloseContacts();
    //////////////////////////////////////////////////////////
    //                       ATTRIBUTES                     //
    //////////////////////////////////////////////////////////
    Residue* residue_; /*!< A pointer back to the residue for this glycosite >*/
    Carbohydrate* glycan_;
    // ResidueLinkage residueGlycanLinkage_;
    //    glycan_overlap_ = 0.0;
    //    double protein_overlap_ = 0.0;
    std::vector<GlycosylationSite*> other_glycosites_;
    std::vector<Residue*> otherProteinResidues_;
    unsigned long int internalBondCount_ = 0;
    //  std::vector<Atom*> proteinAtoms_;
};
#endif // GMML_INCLUDES_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_GLYCOSYLATIONSITE_HPP
