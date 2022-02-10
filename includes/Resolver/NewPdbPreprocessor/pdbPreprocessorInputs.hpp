#ifndef GMML_INCLUDES_RESOLVER_NEWPDBPREPROCSSOR_INPUTS_HPP
#define GMML_INCLUDES_RESOLVER_NEWPDBPREPROCSSOR_INPUTS_HPP

#include <string>
#include <vector>
#include "includes/utils.hpp" // gmml::Split
#include "includes/CodeUtils/logging.hpp"

// These structs are all for interacting with the website via gems.

namespace pdb
{
struct PreprocessorOptions // This is inputs into the preprocessor function. Defaults are fine, can contain user options instead/aswell.
{
    //Constructors
    PreprocessorOptions() :
        chainNTermination_("NH3+"), // aka zwitterionic
        chainCTermination_("CO2-"), // aka zwitterionic
        gapNTermination_("COCH3"), // aka ACE
        gapCTermination_("NHCH3") {} // aka NME
    PreprocessorOptions(
            std::vector<std::pair<std::string,std::string>> hisSelections,
            std::string chainNTermination = "NH3+",
            std::string chainCTermination = "CO2-",
            std::string gapNTermination = "COCH3",
            std::string gapCTermination = "NHCH3"
    ) :
        chainNTermination_(chainNTermination),
        chainCTermination_(chainCTermination),
        gapNTermination_(gapNTermination),
        gapCTermination_(gapCTermination),
        hisSelections_(hisSelections) {}
    // Members
    std::string chainNTermination_;
    std::string chainCTermination_;
    std::string gapNTermination_;
    std::string gapCTermination_;
    std::vector<std::pair<std::string,std::string>> hisSelections_; // e.g. pair: residue id like this <"HIS_20_?_A_1", "HID">
};

struct ResidueId
{
    //Constructor
    ResidueId(std::string inputId)
    {
        std::vector<std::string> tokens = gmml::Split(inputId, "_");
        if (tokens.size() == 5)
        {
            inputId_ = inputId;
            name_ = tokens.at(0);
            number_ = tokens.at(1);
            insertionCode_ = tokens.at(2);
            chainId_ = tokens.at(3);
            model_ = tokens.at(4);
        }
        else
        {
            gmml::log(__LINE__,__FILE__,gmml::ERR, "Could not parse residue id info from this id: " + inputId);
        }
    }
    // Members
    std::string inputId_;
    std::string name_;
    std::string number_;
    std::string insertionCode_;
    std::string chainId_;
    std::string model_;
};

struct DisulphideBond // The original spelling is phabulous.
{
    //Constructor
    DisulphideBond(std::string res1Id, std::string res2Id, double distance) : residue1_(res1Id), residue2_(res2Id), distance_(distance) {}
    //Members
    ResidueId residue1_;
    ResidueId residue2_;
    double distance_;
};

struct GapInAminoAcidChain
{
    //Constructor
    GapInAminoAcidChain(std::string chain, std::string resBefore, std::string resAfter, std::string cTerm, std::string nterm) : chainId_(chain), residueBeforeGap_(resBefore), residueAfterGap_(resAfter), terminationBeforeGap_(cTerm), terminationAfterGap_(nterm) {}
    //Members
    std::string chainId_;
    std::string residueBeforeGap_;
    std::string residueAfterGap_;
    std::string terminationBeforeGap_;
    std::string terminationAfterGap_;
};

//struct AtomId
//{
//    //Constructor
//    AtomId(std::string inputId)
//    {
//        std::vector<std::string> tokens = gmml::Split(inputId, "_");
//        if (tokens.size() == 7)
//        {
//            inputId_ = inputId;
//            name_ = tokens.at(0);
//            number_ = tokens.at(1);
//            resName_ = tokens.at(2);
//            resNumber_ = tokens.at(3);
//            resInsertionCode_ = tokens.at(4);
//            chainId_ = tokens.at(5);
//            model_ = tokens.at(6);
//        }
//        else
//        {
//            gmml::log(__LINE__,__FILE__,gmml::ERR, "Could not parse atom id info from this id: " + inputId);
//        }
//    }
//    // Members
//    std::string inputId_;
//    std::string name_;
//    std::string number_;
//    std::string resName_;
//    std::string resNumber_;
//    std::string resInsertionCode_;
//    std::string chainId_;
//    std::string model_;
//};

struct AtomInfo
{
    // Constructor
    AtomInfo(std::string atomName, ResidueId residue) : name_(atomName), residue_(residue) {}
    // Members
    std::string name_;
    ResidueId residue_;
};


struct ChainTerminal
{
    // Constructor
    ChainTerminal(std::string chainId, std::string startIndex, std::string endIndex, std::string nTerm, std::string cTerm) : chainId_(chainId), startIndex_(startIndex), endIndex_(endIndex), nTermination_(nTerm), cTermination_(cTerm) {}
    // Members
    std::string chainId_;
    std::string startIndex_;
    std::string endIndex_;
    std::string nTermination_;
    std::string cTermination_;
};

struct PreprocessorInformation
{
    std::vector<AtomInfo> unrecognizedAtoms_;
    std::vector<AtomInfo> missingHeavyAtoms_;
    std::vector<ResidueId> unrecognizedResidues_;
    std::vector<GapInAminoAcidChain> missingResidues_;
    std::vector<ChainTerminal> chainTerminals_;
    std::vector<ResidueId> hisResidues_;
    std::vector<DisulphideBond> cysBondResidues_;
};

}
#endif // GMML_INCLUDES_RESOLVER_NEWPDBPREPROCSSORINPUTS_HPP
