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

    //The below function has no code in the .cc file, so GEMS breaks.  Commented out for now -DM
    // std::string GetCodeForResidue(std::string query);
    std::string GetCodeForLinkage(std::string query);

private:
    std::vector<LinkageCodes> Glycam06LinkageCodesLookup_;
};
} // close namespace
} // close namespace
} // close namespace
#endif // GLYCAM06_LINKAGE_CODES_HPP
