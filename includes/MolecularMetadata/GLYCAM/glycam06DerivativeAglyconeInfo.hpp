#ifndef GMML_INCLUDES_MOLECULARMETADATA_GLYCAM_GLYCAM06DERIVATIVEAGLYCONEINFO_HPP
#define GMML_INCLUDES_MOLECULARMETADATA_GLYCAM_GLYCAM06DERIVATIVEAGLYCONEINFO_HPP

#include <string>
#include <map>
#include <vector>

namespace gmml
{
namespace MolecularMetadata
{
namespace GLYCAM
{

class Glycam06DerivativeAglyconeConnectionAtomLookup
{
public:

    //////////////////////////////////////////////////////////
    //                       CONSTRUCTOR                    //
    //////////////////////////////////////////////////////////
    /*! \fn
    * Default constructor
    */
    Glycam06DerivativeAglyconeConnectionAtomLookup();

    //////////////////////////////////////////////////////////
    //                         TYPEDEFS                     //
    //////////////////////////////////////////////////////////


    //////////////////////////////////////////////////////////
    //                      QUERY FUNCTIONS                 //
    //////////////////////////////////////////////////////////

    inline std::string GetConnectionAtomForResidue(std::string query)
    {
        for (auto &elem : glycam06DerivativeAglyconeConnectionAtomLookup_)
        {
            if (elem.first == query)
            {
                return elem.second;
            }
        }
        return "Derivative or aglycone residue is not currently supported by GLYCAM.";
    }
private:
    std::multimap<std::string, std::string> glycam06DerivativeAglyconeConnectionAtomLookup_;
};
} // close namespace
} // close namespace
} // close namespace

#endif // GMML_INCLUDES_MOLECULARMETADATA_GLYCAM_GLYCAM06DERIVATIVEAGLYCONEINFO_HPP
