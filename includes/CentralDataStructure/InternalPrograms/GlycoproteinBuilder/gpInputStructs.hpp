#ifndef GMML_INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_GPINPUTSTRUCTS_HPP
#define GMML_INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_GPINPUTSTRUCTS_HPP

#include <string>
#include <vector>
#include <iostream>

namespace glycoprotein
{
    struct GlycositeInput
    {
        // Constructor
        GlycositeInput(std::string proteinResidueId, std::string glycan)
            : proteinResidueId_(proteinResidueId), glycanInput_(glycan)
        {}

        std::string proteinResidueId_ = ""; // E.g. ?_20_A if no chain ID and residue number is 20 and insertion code is
                                            // A. C_20_? if chain id is C and there is no insertion code.
        std::string glycanInput_ =
            ""; // E.g. Man9 if "Library" glycanInputType. E.g. DGlcpNAcb1-4DGlcpNAcb1-OH if "Sequence".
    };

    struct GlycoproteinBuilderInputs
    {
        std::string substrateFileName_ = "Undefined"; // Program should throw if left as "Undefined".
        int number3DStructures_        = 1;           // ToDo Implement this.
        int maxThreads_                = 1;           // ToDo Implement this.
        int persistCycles_             = 5;
        int overlapTolerance_          = 1;
        bool isDeterministic_          = false;
        std::vector<GlycositeInput> glycositesInputVector_; // No default, program will throw if uninitialized.
    };

    GlycoproteinBuilderInputs readGPInputFile(std::string inputFileName);
} // namespace glycoprotein
#endif // GMML_INCLUDES_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_GPINPUTSTRUCTS_HPP
