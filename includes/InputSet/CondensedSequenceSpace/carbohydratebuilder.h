#ifndef CARBOHYDRATEBUILDER_H
#define CARBOHYDRATEBUILDER_H
#include "condensedsequence.hpp"
#include "../../MolecularModeling/assembly.hpp"

namespace CondensedSequenceSpace
{
class carbohydrateBuilder
{
public:

    //////////////////////////////////////////////////////////
    //                       CONSTRUCTORS                   //
    //////////////////////////////////////////////////////////

  //  carbohydrateBuilder();
    carbohydrateBuilder(std::string selectedBuildType = "build", std::string condensedSequence = "DManp[2S,3Me]a1-6DManpa1-6[DGlcpNAcb1-4][DNeu5Aca2-6DGalpb1-4DGlcpNAc[3S]b1-2DManpa1-3]DManpb1-4DGlcpNAc[6Me]b1-4DGlcpNAcb1-OH", std::string prepFilePath = "/home/oliver/Programs/GLYCAM_Dev_Env/V_2/Web_Programs/gems/gmml/dat/prep/GLYCAM_06j-1.prep");

    //////////////////////////////////////////////////////////
    //                       ACCESSORS                      //
    //////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////
    //                       MUTATOR                        //
    //////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////
    //                        FUNCTIONS                     //
    //////////////////////////////////////////////////////////
    void SetDefaultDihedralAngleGeometryWithMetadata();
    void ResolveOverlaps();


private:
    void FigureOutResidueLinkagesInGlycan(MolecularModeling::Residue *from_this_residue1, MolecularModeling::Residue *to_this_residue2, ResidueLinkageVector *residue_linkages);
    void InitializeClass(std::string selectedBuildType = "build", std::string condensedSequence = "DManp[2S,3Me]a1-6DManpa1-6[DGlcpNAcb1-4][DNeu5Aca2-6DGalpb1-4DGlcpNAc[3S]b1-2DManpa1-3]DManpb1-4DGlcpNAc[6Me]b1-4DGlcpNAcb1-OH", std::string prepFilePath = "/home/oliver/Programs/GLYCAM_Dev_Env/V_2/Web_Programs/gems/gmml/dat/prep/GLYCAM_06j-1.prep");

    //////////////////////////////////////////////////////////
    //                       ATTRIBUTES                     //
    //////////////////////////////////////////////////////////

    MolecularModeling::Assembly assembly_;
    CondensedSequence sequence_;
    //GlycanMetadataContainer metadataInformation_; // need this class
    ResidueLinkageVector glycosidicLinkages_;

};
}
#endif // CARBOHYDRATEBUILDER_H



// Proposed usage examples:

//carbohydrateBuilder builder(std::string sequenceString, std::string prepFilePath, enum jobTypeOption);

// Scrap
