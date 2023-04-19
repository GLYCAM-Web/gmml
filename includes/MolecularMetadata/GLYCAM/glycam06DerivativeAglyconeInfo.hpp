#ifndef GLYCAM06_DERIVATIVE_AGLYCONE_CONNECTION_ATOMS_HPP
#define GLYCAM06_DERIVATIVE_AGLYCONE_CONNECTION_ATOMS_HPP

#include "includes/CodeUtils/logging.hpp"
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

    inline std::string GetConnectionAtomForResidue(const std::string query) const
    {
        for (auto &elem : glycam06DerivativeAglyconeConnectionAtomLookup_)
        {
            if (elem.first == query)
            {
                return elem.second;
            }
        }
        std::string message = "The selected derivative or aglycone residue is not currently supported by GLYCAM: " + query;
        gmml::log(__LINE__,__FILE__, gmml::ERR, message);
        throw std::runtime_error(message);
        return "";
    }
private:
    std::multimap<std::string, std::string> glycam06DerivativeAglyconeConnectionAtomLookup_;
};
} // close namespace
} // close namespace
} // close namespace

#endif // GLYCAM06_DERIVATIVE_AGLYCONE_CONNECTION_ATOMS_HPP
