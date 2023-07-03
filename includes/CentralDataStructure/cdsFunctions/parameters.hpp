#ifndef INCLUDES_CENTRALDATASTRUCTURE_CDSFUNCTIONS_PARAMETERS_HPP_
#define INCLUDES_CENTRALDATASTRUCTURE_CDSFUNCTIONS_PARAMETERS_HPP_

#include "includes/CentralDataStructure/residue.hpp"
#include "includes/CentralDataStructure/Readers/Prep/prepFile.hpp"
#include "includes/CentralDataStructure/Readers/Lib/LibraryFile.hpp"
#include "includes/CodeUtils/logging.hpp"
#include <vector>
#include <unordered_map>

namespace cdsParameters
{
    class ParameterManager
    {
      public:
        // Constructor
        ParameterManager(); // Take in nothing for now, just use defaults
        // Functions
        cds::Residue* getParameterResidue(const std::string name) const;
        bool setAtomCharges(cds::Residue* queryResidue);
        void setAtomCharges(std::vector<cds::Residue*> queryResidues);

      private:
        // Functions
        void InitializeResidueMap(std::vector<cds::Residue*> incomingResidues);
        // Attributes
        std::unordered_map<std::string, cds::Residue*> parameterResidueMap_;
        std::vector<prep::PrepFile> prepFiles_;
        std::vector<lib::LibraryFile> libFiles_;
    };

    // Separate cause it can be
    bool setChargeForAtom(cds::Atom* queryAtom, std::vector<cds::Atom*> referenceAtoms);
    // bool setCharges(cds::Residue* queryResidue, std::vector<cds::Residue*> referenceResidues);
    void setCharges(std::vector<cds::Residue*> queryResidues);
} // namespace cdsParameters
#endif /* INCLUDES_CENTRALDATASTRUCTURE_CDSFUNCTIONS_PARAMETERS_HPP_ */
