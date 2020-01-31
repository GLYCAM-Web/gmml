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
    std::string residueName_ ;
    std::string glycamCode_ ;
    std::string residueType_ ;
} ;

typedef std::vector<ResidueNamesCodesTypes> ResidueNamesCodesTypesVector;

class Glycam06ResidueNamesToCodesLookupContainer
{
public:

    //////////////////////////////////////////////////////////
    //                       CONSTRUCTOR                    //
    //////////////////////////////////////////////////////////
    /*! \fn
    * Default constructor
    */
    Glycam06ResidueNamesToCodesLookupContainer();

    //////////////////////////////////////////////////////////
    //                         TYPEDEFS                     //
    //////////////////////////////////////////////////////////


    //////////////////////////////////////////////////////////
    //                      QUERY FUNCTIONS                 //
    //////////////////////////////////////////////////////////

    //The below function has no code in the .cc file, so GEMS breaks.  Commented out for now -DM
    // std::string GetCodeForResidue(std::string query);
    std::string GetResidueForCode(std::string query);

private:
    ResidueNamesCodesTypesVector ResidueNamesCodesTypesLookupTable_;
};
} // close namespace
} // close namespace
} // close namespace
#endif // GLYCAM06RESIDUECODES_HPP
