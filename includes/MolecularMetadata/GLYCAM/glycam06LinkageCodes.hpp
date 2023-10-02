#ifndef GLYCAM06_LINKAGE_CODES_HPP
#define GLYCAM06_LINKAGE_CODES_HPP
#include <string>
#include <map>
#include <vector>

namespace gmml
{
    namespace MolecularMetadata
    {
        namespace GLYCAM
        {
            struct LinkageCodes
            {
                std::string linkage_;
                std::string glycamLinkageCode_;
            };

            class Glycam06LinkageCodesLookupContainer
            {
              public:
                //////////////////////////////////////////////////////////
                //                       CONSTRUCTOR                    //
                //////////////////////////////////////////////////////////
                /*! \fn
                 * Default constructor
                 */
                Glycam06LinkageCodesLookupContainer();

                //////////////////////////////////////////////////////////
                //                         TYPEDEFS                     //
                //////////////////////////////////////////////////////////

                //////////////////////////////////////////////////////////
                //                      QUERY FUNCTIONS                 //
                //////////////////////////////////////////////////////////

                std::string GetCodeForLinkages(std::string query);

              private:
                std::vector<LinkageCodes> Glycam06LinkageCodesLookup_;
            };
        } // namespace GLYCAM
    }     // namespace MolecularMetadata
} // namespace gmml
#endif // GLYCAM06_LINKAGE_CODES_HPP
