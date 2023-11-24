#ifndef GLYCAM06RESIDUECODES_HPP
#define GLYCAM06RESIDUECODES_HPP
#include <string>
#include <map>
#include <vector>

namespace gmml
{
    namespace MolecularMetadata
    {
        namespace GLYCAM
        {
            struct ResidueNamesCodesTypes
            {
                std::string residueName_;
                std::string glycamCode_;
                std::string residueType_;
            };

            typedef std::vector<ResidueNamesCodesTypes> ResidueNamesCodesTypesVector;

            class Glycam06ResidueNamesToCodesLookupContainer
            {
              public:
                //////////////////////////////////////////////////////////
                //                       CONSTRUCTOR                    //
                //////////////////////////////////////////////////////////
                Glycam06ResidueNamesToCodesLookupContainer();
                //////////////////////////////////////////////////////////
                //                         TYPEDEFS                     //
                //////////////////////////////////////////////////////////

                //////////////////////////////////////////////////////////
                //                      QUERY FUNCTIONS                 //
                //////////////////////////////////////////////////////////
                std::string GetResidueForCode(std::string query) const;
                std::string GetCodeForResidue(std::string query) const;

              private:
                ResidueNamesCodesTypesVector ResidueNamesCodesTypesLookupTable_;
            };
        } // namespace GLYCAM
    }     // namespace MolecularMetadata
} // namespace gmml
#endif // GLYCAM06RESIDUECODES_HPP
