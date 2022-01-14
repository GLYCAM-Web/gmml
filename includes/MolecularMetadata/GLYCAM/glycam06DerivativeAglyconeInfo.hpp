#ifndef GLYCAM06_DERIVATIVE_AGLYCONE_CONNECTION_ATOMS_HPP
#define GLYCAM06_DERIVATIVE_AGLYCONE_CONNECTION_ATOMS_HPP

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

#endif // GLYCAM06_DERIVATIVE_AGLYCONE_CONNECTION_ATOMS_HPP
