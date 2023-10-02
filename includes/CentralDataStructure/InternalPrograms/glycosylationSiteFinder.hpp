#include <iostream>
#include <string>
#include <vector>
#include "includes/CentralDataStructure/residue.hpp"

// Structs
struct GlycosylationSiteInfo
{
    // Constructor
    GlycosylationSiteInfo(std::string chain, std::string residueNumber, std::string insertionCode,
                          std::string sequenceContext, std::vector<std::string> tags)
        : chain_(chain), residueNumber_(residueNumber), insertionCode_(insertionCode),
          sequenceContext_(sequenceContext), tags_(tags)
    {}

    // Data
    std::string chain_;
    std::string residueNumber_;
    std::string insertionCode_;
    std::string sequenceContext_;
    std::vector<std::string> tags_; // e.g. oLink, nLink, sequon, cysteineSequon, all, etc

    // Print
    inline std::string Print()
    {
        std::string output = "Chain: " + chain_ + "\nResidueNumber: " + residueNumber_ +
                             "\nInsertionCode: " + insertionCode_ + "\nSequenceContext: " + sequenceContext_ +
                             "\nTags:\n";
        for (auto& tag : tags_)
        {
            output += tag + "\n";
        }
        return output;
    }
};

namespace glycoproteinBuilder
{
    class GlycosylationSiteFinder
    {
      public:
        GlycosylationSiteFinder(std::vector<cds::Residue*> residues);

        inline std::vector<GlycosylationSiteInfo> GetTable()
        {
            return table_;
        }

        std::string PrintTable();

      private:
        std::vector<GlycosylationSiteInfo> table_;
    };
} // namespace glycoproteinBuilder
