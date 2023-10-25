#ifndef GMML_INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_GLYCOSYLATIONSITE_HPP
#define GMML_INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_GLYCOSYLATIONSITE_HPP

#include "includes/CentralDataStructure/CondensedSequence/carbohydrate.hpp"
#include "includes/CentralDataStructure/assembly.hpp"
#include "includes/CentralDataStructure/Shapers/residueLinkage.hpp"

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
using cdsCondensedSequence::Carbohydrate;

class GlycosylationSite
{
  public:
    //////////////////////////////////////////////////////////
    //                       CONSTRUCTOR                    //
    //////////////////////////////////////////////////////////
    GlycosylationSite(Residue* residue, Carbohydrate* carbohydrate, std::vector<Residue*> otherProteinResidues,
                      unsigned int glycanStartResidueNumber);

    //////////////////////////////////////////////////////////
    //                       ACCESSOR                       //
    //////////////////////////////////////////////////////////
    inline std::string GetResidueId()
    {
        return this->GetResidue()->getStringId();
    }

    //////////////////////////////////////////////////////////
    //                       MUTATOR                        //
    //////////////////////////////////////////////////////////
    void SetOtherGlycosites(std::vector<GlycosylationSite>& glycosites);
    //////////////////////////////////////////////////////////
    //                       FUNCTIONS                      //
    //////////////////////////////////////////////////////////
    unsigned int CountOverlaps(MoleculeType moleculeType = MoleculeType::ALL);
    unsigned int CountOverlapsFast();
    //    void StashCoordinates();
    //    void SetStashedCoordinates();
    void Wiggle(bool firstLinkageOnly = false, int interval = 5, bool useAllResiduesForOverlap = false);
    void SetRandomDihedralAnglesUsingMetadata();
    void ResetDihedralAngles();
    void Rename_Protein_Residue_From_GLYCAM_To_Standard();
    void AddOtherGlycositesToLinkageOverlapAtoms();
    void UpdateOverlapAtomsInLinkages(unsigned int maxProteinResidues = 20);
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
    inline std::vector<GlycosylationSite*>& GetOtherGlycosites()
    {
        return other_glycosites_;
    }

    inline std::vector<Residue*>& GetOtherProteinResidues()
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

    //////////////////////////////////////////////////////////
    //                  PRIVATE FUNCTIONS                   //
    //////////////////////////////////////////////////////////
    void AttachGlycan(unsigned int glycanResidueStartNumber);
    void RenumberGlycanToMatch(unsigned int startNumber);
    void Prepare_Glycans_For_Superimposition_To_Particular_Residue(std::string amino_acid_name);
    void Superimpose_Glycan_To_Glycosite(Residue* glycosite_residue);
    void Rename_Protein_Residue_To_GLYCAM_Nomenclature();
    Atom* GetConnectingProteinAtom(const std::string residue_name) const;
    void WiggleOneLinkage(ResidueLinkage& linkage, int interval = 5, bool useAllResiduesForOverlap = false);
    unsigned int CountOverlaps(const std::vector<Residue*>& residuesA, const std::vector<Residue*>& residuesB);
    //////////////////////////////////////////////////////////
    //                       ATTRIBUTES                     //
    //////////////////////////////////////////////////////////
    Residue* residue_; /*!< A pointer back to the residue for this glycosite >*/
    Carbohydrate* glycan_;
    std::vector<GlycosylationSite*> other_glycosites_;
    std::vector<Residue*> otherProteinResidues_;
};
#endif // GMML_INCLUDES_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_GLYCOSYLATIONSITE_HPP
