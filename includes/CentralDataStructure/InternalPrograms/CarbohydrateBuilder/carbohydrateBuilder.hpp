#ifndef GMML_INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_CARBOHYDRATEBUILDER_CARBOHYDRATEBUILDER_HPP
#define GMML_INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_CARBOHYDRATEBUILDER_CARBOHYDRATEBUILDER_HPP

// This is becoming an interface to the carbohydrate class in gmml for Gems.
#include "includes/CentralDataStructure/CondensedSequence/carbohydrate.hpp"
#include <string>
#include <vector>

namespace cdsCondensedSequence
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
    DihedralOptions(const std::string name, const std::vector<std::string> rotamers) : dihedralName_ (name), rotamers_ (rotamers) {}
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

class carbohydrateBuilder : public cdsCondensedSequence::Carbohydrate
{
public:
    //////////////////////////////////////////////////////////
    //                       CONSTRUCTORS                   //
    //////////////////////////////////////////////////////////
    carbohydrateBuilder(std::string condensedSequence = "DManp[2S,3Me]a1-6DManpa1-6[DGlcpNAcb1-4][DNeup5Aca2-6DGalpb1-4DGlcpNAc[3S]b1-2DManpa1-3]DManpb1-4DGlcpNAc[6Me]b1-4DGlcpNAcb1-OH", std::string prepFilePath = "../dat/prep/GLYCAM_06j-1_GAGS_KDN.prep") : cdsCondensedSequence::Carbohydrate(condensedSequence, prepFilePath) {}

    //////////////////////////////////////////////////////////
    //                       ACCESSORS                      //
    //////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////
    //                       MUTATOR                        //
    //////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////
    //                        FUNCTIONS                     //
    //////////////////////////////////////////////////////////
    void GenerateSingle3DStructureDefaultFiles(std::string fileOutputDirectory = "unspecified", std::string outputFileNaming = "structure");
    LinkageOptionsVector GenerateUserOptionsDataStruct();
    void GenerateSpecific3DStructure(SingleRotamerInfoVector conformerInfo, std::string fileOutputDirectory = "unspecified");
    void GenerateUpToNRotamers(int maxRotamers = 32); // Will not be used by gems, but leaving the functionality as could be useful.
  //  std::string GetNumberOfShapes(bool likelyShapesOnly = false);
private:
    void generateLinkagePermutationsRecursively(std::vector<cds::ResidueLinkage>::iterator linkage, std::vector<cds::ResidueLinkage>::iterator end, int maxRotamers = 32, int rotamerCount = 0);
//  void resetLinkageIDsToStartFromZero(ResidueLinkageVector &inputLinkages);
    std::string convertIncomingRotamerNamesToStandard(std::string incomingName);
    //////////////////////////////////////////////////////////
    //                       ATTRIBUTES                     //
    //////////////////////////////////////////////////////////
    cdsCondensedSequence::Carbohydrate carbohydrate_;
};
}
#endif // GMML_INCLUDES_INTERNALPROGRAMS_CARBOHYDRATEBUILDER_CARBOHYDRATEBUILDER_HPP
