#ifndef GMML_INCLUDES_RESOLVER_NEWPDBPREPROCSSOR_INPUTS_HPP
#define GMML_INCLUDES_RESOLVER_NEWPDBPREPROCSSOR_INPUTS_HPP

#include <string>
#include <vector>

#include "includes/CentralDataStructure/Readers/Pdb/pdbResidueId.hpp" // residueId
// #include "includes/utils.hpp" // gmml::Split
#include "includes/CodeUtils/logging.hpp"

// These structs are all for interacting with the website via gems.

namespace pdb
{
    struct PreprocessorOptions // This is inputs into the preprocessor function. Defaults are fine, can contain user
                               // options instead/aswell.
    {
        // Constructors
        PreprocessorOptions()
        {}

        PreprocessorOptions(std::vector<std::pair<std::string, std::string>> hisSelections,
                            std::string chainNTermination = "NH3+", std::string chainCTermination = "CO2-",
                            std::string gapNTermination = "COCH3", std::string gapCTermination = "NHCH3")
            : chainNTermination_(chainNTermination), chainCTermination_(chainCTermination),
              gapNTermination_(gapNTermination), gapCTermination_(gapCTermination), hisSelections_(hisSelections)
        {}

        // Members
        std::string chainNTermination_ = "NH3+";  // aka zwitterionic
        std::string chainCTermination_ = "CO2-";  // aka zwitterionic
        std::string gapNTermination_   = "COCH3"; // aka ACE
        std::string gapCTermination_   = "NHCH3"; // aka NME
        std::vector<std::pair<std::string, std::string>>
            hisSelections_; // e.g. pair: residue id like this <"HIS_20_?_A_1", "HID">
    };

    // struct pdb::ResidueId
    //{
    //     //Constructor
    //     pdb::ResidueId(std::string inputId)
    //     {
    //         std::vector<std::string> tokens = gmml::Split(inputId, "_");
    //         if (tokens.size() == 4)
    //         {
    //             inputId_ = inputId;
    //             name_ = tokens.at(0);
    //             number_ = tokens.at(1);
    //             insertionCode_ = tokens.at(2);
    //             chainId_ = tokens.at(3);
    ////            model_ = tokens.at(4);
    //        }
    //        else
    //        {
    //            gmml::log(__LINE__,__FILE__,gmml::ERR, "Could not parse residue id info from this id: " + inputId);
    //        }
    //    }
    //    // Members
    //    std::string inputId_;
    //    std::string name_;
    //    std::string number_;
    //    std::string insertionCode_;
    //    std::string chainId_;
    ////    std::string model_;
    //};

    struct DisulphideBond // The original spelling is phabulous.
    {
        // Constructor
        DisulphideBond()
        {} // SWIG needs this

        DisulphideBond(const pdb::ResidueId& res1Id, const pdb::ResidueId& res2Id, const double& distance)
            : residue1_(res1Id), residue2_(res2Id), distance_(distance)
        {}

        // Members
        pdb::ResidueId residue1_;
        pdb::ResidueId residue2_;
        double distance_ = constants::dNotSet;
    };

    struct GapInAminoAcidChain
    {
        // Constructor
        GapInAminoAcidChain()
        {} // SWIG needs this

        GapInAminoAcidChain(const std::string& chain, const std::string& resBefore, const std::string& resAfter,
                            const std::string& cTerm, const std::string& nterm)
            : chainId_(chain), residueBeforeGap_(resBefore), residueAfterGap_(resAfter), terminationBeforeGap_(cTerm),
              terminationAfterGap_(nterm)
        {}

        // Members
        std::string chainId_;
        std::string residueBeforeGap_;
        std::string residueAfterGap_;
        std::string terminationBeforeGap_;
        std::string terminationAfterGap_;
    };

    // struct AtomId
    //{
    //     //Constructor
    //     AtomId(std::string inputId)
    //     {
    //         std::vector<std::string> tokens = gmml::Split(inputId, "_");
    //         if (tokens.size() == 7)
    //         {
    //             inputId_ = inputId;
    //             name_ = tokens.at(0);
    //             number_ = tokens.at(1);
    //             resName_ = tokens.at(2);
    //             resNumber_ = tokens.at(3);
    //             resInsertionCode_ = tokens.at(4);
    //             chainId_ = tokens.at(5);
    //             model_ = tokens.at(6);
    //         }
    //         else
    //         {
    //             gmml::log(__LINE__,__FILE__,gmml::ERR, "Could not parse atom id info from this id: " + inputId);
    //         }
    //     }
    //     // Members
    //     std::string inputId_;
    //     std::string name_;
    //     std::string number_;
    //     std::string resName_;
    //     std::string resNumber_;
    //     std::string resInsertionCode_;
    //     std::string chainId_;
    //     std::string model_;
    // };

    struct AtomInfo
    {
        // Constructor
        AtomInfo()
        {} // required for SWIG? Really SWIG?

        AtomInfo(const std::string& atomName, const pdb::ResidueId& residue) : name_(atomName), residue_(residue)
        {}

        // Members
        std::string name_;
        pdb::ResidueId residue_;
    };

    struct ChainTerminal
    {
        // Constructor
        ChainTerminal()
        {} // SWIG needs this

        ChainTerminal(const std::string& chainId, const std::string& startIndex, const std::string& endIndex,
                      const std::string& nTerm, const std::string& cTerm)
            : chainId_(chainId), startIndex_(startIndex), endIndex_(endIndex), nTermination_(nTerm),
              cTermination_(cTerm)
        {}

        // Members
        std::string chainId_;
        std::string startIndex_;
        std::string endIndex_;
        std::string nTermination_;
        std::string cTermination_;
    };

    struct NonNaturalProteinResidue
    {
        // Constructor
        NonNaturalProteinResidue()
        {}

        NonNaturalProteinResidue(const pdb::ResidueId& residueId) : residue_(residueId)
        {}

        // Members
        pdb::ResidueId residue_;
    };

    struct PreprocessorInformation
    {
        std::vector<AtomInfo> unrecognizedAtoms_;
        std::vector<AtomInfo> missingHeavyAtoms_;
        std::vector<pdb::ResidueId> unrecognizedResidues_;
        std::vector<GapInAminoAcidChain> missingResidues_;
        std::vector<ChainTerminal> chainTerminals_;
        std::vector<pdb::ResidueId> hisResidues_;
        std::vector<DisulphideBond> cysBondResidues_;
        std::vector<NonNaturalProteinResidue> nonNaturalProteinResidues_;
    };

} // namespace pdb
#endif // GMML_INCLUDES_RESOLVER_NEWPDBPREPROCSSORINPUTS_HPP
