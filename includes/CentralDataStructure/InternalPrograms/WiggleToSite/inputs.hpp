#ifndef GMML_INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_WIGGLETOSITE_INPUTS_HPP
#define GMML_INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_WIGGLETOSITE_INPUTS_HPP
#include "includes/CentralDataStructure/Readers/Pdb/pdbResidueId.hpp"
#include <string>

namespace gmmlPrograms
{
    struct WiggleToSiteInputs
    {
        // ctor
        WiggleToSiteInputs(std::string inputFileName);
        // Members
        std::string carbohydrateSequence_       = "";
        int carbohydrateSuperimpositionResidue_ = 0;
        int carbohydrateWigglingResidue_        = 0;
        std::string substrateFile_              = "";
        pdb::ResidueId superimpositionTargetResidue_;
        pdb::ResidueId wigglingTargetResidue_;
        int substrateModelNumber_ = 1;
        int persistCycles_        = 5;
        bool isDeterministic_     = false;
        std::string Print();
    };
} // namespace gmmlPrograms
#endif // GMML_INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_WIGGLETOSITE_INPUTS_HPP
