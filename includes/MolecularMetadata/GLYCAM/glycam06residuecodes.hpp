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
/**
         *   Glycam06NamesToCodesLookupMap
         *
         *   The first string is the name-code for a residue.  It is typically
         *   three or six characters long.
         *
         *     Examples:  Gal, Neu5Ac
         *
         *   The second string is a code used by GLYCAM in PDB files for that residue.
         *
         *     Examples:  Gal is G
         *
         *
         */

class Glycam06ResidueNamesToCodesLookup
{
public:

    //////////////////////////////////////////////////////////
    //                       CONSTRUCTOR                    //
    //////////////////////////////////////////////////////////
    /*! \fn
    * Default constructor
    */
    Glycam06ResidueNamesToCodesLookup();

    //////////////////////////////////////////////////////////
    //                         TYPEDEFS                     //
    //////////////////////////////////////////////////////////


    //////////////////////////////////////////////////////////
    //                      QUERY FUNCTIONS                 //
    //////////////////////////////////////////////////////////

    std::string GetCodeForResidue(std::string query);
    std::string GetResidueForCode(std::string query);

private:
    std::multimap<std::string, std::string> glycam06ResidueNamesToCodesLookupMap_;
};
} // close namespace
} // close namespace
} // close namespace
#endif // GLYCAM06RESIDUECODES_HPP
