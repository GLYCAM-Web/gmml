#ifndef GMML_INCLUDES_INTERNALPROGRAMS_CARBOHYDRATEBUILDER_CARBOHYDRATEBUILDER_HPP
#define GMML_INCLUDES_INTERNALPROGRAMS_CARBOHYDRATEBUILDER_CARBOHYDRATEBUILDER_HPP
#include <string>
#include <vector>
#include "includes/MolecularModeling/assembly.hpp"
#include "includes/Abstract/absBuilder.hpp"

namespace CondensedSequence
{ // For specifying a specific shape to be built with GenerateSpecific3DStructure
struct SingleRotamerInfo
{
    std::string linkageIndex; // What Dan is calling linkageLabel. Internal index determined at C++ level and given to frontend to track.
    std::string linkageName; // Can be whatever the user wants it to be, default to same as index.
    std::string dihedralName; // omg / phi / psi / chi1 / chi2
    std::string selectedRotamer; // gg / tg / g- etc
    std::string numericValue; // user entered 64 degrees. Could be a v2 feature.
};
typedef std::vector<SingleRotamerInfo> SingleRotamerInfoVector;

struct DihedralOptions
{   // CONSTRUCTOR
    DihedralOptions () {}
    DihedralOptions(std::string name, std::vector<std::string> rotamers) : dihedralName_ (name), rotamers_ (rotamers) {}
    // DATA
    std::string dihedralName_; // omg / phi / psi / chi1 / chi2
    std::vector<std::string> rotamers_; // gg / tg / g- etc
};
typedef std::vector<DihedralOptions> DihedralOptionsVector;

struct LinkageOptions
{   // CONSTRUCTOR
    LinkageOptions () {}
    LinkageOptions(std::string name, std::string index, std::string res1, std::string res2, DihedralOptionsVector likely, DihedralOptionsVector possible)
                    : linkageName_ (name), indexOrderedLabel_ (index), firstResidueNumber_ (res1), secondResidueNumber_ (res2),
                      likelyRotamers_ (likely), possibleRotamers_ (possible) {}
    // DATA
    std::string linkageName_;
    std::string indexOrderedLabel_;
    std::string firstResidueNumber_;
    std::string secondResidueNumber_;
    DihedralOptionsVector likelyRotamers_;
    DihedralOptionsVector possibleRotamers_;
};
typedef std::vector<LinkageOptions> LinkageOptionsVector;

class carbohydrateBuilder : public Abstract::absBuilder
{
public:
    //////////////////////////////////////////////////////////
    //                       CONSTRUCTORS                   //
    //////////////////////////////////////////////////////////
    carbohydrateBuilder(std::string condensedSequence = "DManp[2S,3Me]a1-6DManpa1-6[DGlcpNAcb1-4][DNeup5Aca2-6DGalpb1-4DGlcpNAc[3S]b1-2DManpa1-3]DManpb1-4DGlcpNAc[6Me]b1-4DGlcpNAcb1-OH", std::string prepFilePath = "../dat/prep/GLYCAM_06j-1_GAGS.prep");
    //////////////////////////////////////////////////////////
    //                       ACCESSORS                      //
    //////////////////////////////////////////////////////////
    std::string GetInputSequenceString();
    MolecularModeling::Assembly* GetAssembly();
    ResidueLinkageVector* GetGlycosidicLinkages(); // Dodgy. Why do I do this?
    //////////////////////////////////////////////////////////
    //                       MUTATOR                        //
    //////////////////////////////////////////////////////////
    // These two being public is funky.
    // For evaluation I don't want to bother with this, so they only trigger when 3D structures are generated,
    // but for gp builder, I don't write 3D structure until the end and so it does it's count internal bonds while the shape isn't good.
    // Making public while I come up with a better idea.
    // Probably carb builder needs an eval constructor, or that might go into sequence.
    // Actually I think these need to be done by default to be safe.
    void SetDefaultShapeUsingMetadata();
    void ResolveOverlaps();
    //////////////////////////////////////////////////////////
    //                        FUNCTIONS                     //
    //////////////////////////////////////////////////////////
    void GenerateSingle3DStructureDefaultFiles(std::string fileOutputDirectory = "unspecified", std::string outputFileNaming = "structure");
    void GenerateSingle3DStructureSingleFile(std::string fileOutputDirectory = "unspecified", std::string fileType = "PDB", std::string outputFileNaming = "structure");
    LinkageOptionsVector GenerateUserOptionsDataStruct();
    void GenerateSpecific3DStructure(SingleRotamerInfoVector conformerInfo, std::string fileOutputDirectory = "unspecified");
    // The following generates no files
    void GenerateUpToNRotamers(int maxRotamers = 32); // Will not be used by gems, but leaving the functionality as could be useful.
    int GetNumberOfShapes(bool likelyShapesOnly = false);
    std::string Print();
private:
    void Write3DStructureFile(std::string fileOutputDirectory = "unspecified", std::string fileType = "PDB", std::string filename = "output");
    void SetInputSequenceString(std::string sequence);
    void FigureOutResidueLinkagesInGlycan(MolecularModeling::Residue *from_this_residue1, MolecularModeling::Residue *to_this_residue2, ResidueLinkageVector *residue_linkages);
    void InitializeClass(std::string condensedSequence);
    // This does not belong in this class:
    ResidueLinkageVector SplitLinkagesIntoPermutants(ResidueLinkageVector &inputLinkages);
    void generateLinkagePermutationsRecursively(ResidueLinkageVector::iterator linkage, ResidueLinkageVector::iterator end, int maxRotamers = 32, int rotamerCount = 0);
    Residue_linkage* selectLinkageWithIndex(ResidueLinkageVector &inputLinkages, int indexQuery);
    void resetLinkageIDsToStartFromZero(ResidueLinkageVector &inputLinkages);
    std::string convertIncomingRotamerNamesToStandard(std::string incomingName);
    //////////////////////////////////////////////////////////
    //                       ATTRIBUTES                     //
    //////////////////////////////////////////////////////////
    MolecularModeling::Assembly assembly_;
    std::string inputSequenceString_;
    ResidueLinkageVector glycosidicLinkages_;
};
}
#endif // GMML_INCLUDES_INTERNALPROGRAMS_CARBOHYDRATEBUILDER_CARBOHYDRATEBUILDER_HPP
