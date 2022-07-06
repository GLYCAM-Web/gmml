#ifndef GMML_INCLUDES_MOLECULARMETADATA_GLYCAM_GLYCAM06LINKAGECODES_HPP
#define GMML_INCLUDES_MOLECULARMETADATA_GLYCAM_GLYCAM06LINKAGECODES_HPP
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
    std::string glycamLinkageCode_ ;
} ;

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
} // close namespace
} // close namespace
} // close namespace
#endif // GMML_INCLUDES_MOLECULARMETADATA_GLYCAM_GLYCAM06LINKAGECODES_HPP
