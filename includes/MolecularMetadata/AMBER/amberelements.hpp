#ifndef AMBER_ELEMENTS_META_HPP
#define AMBER_ELEMENTS_META_HPP

/* File amberelements.hpp begun on 16 June 2018 by BLFoley */

#include <string>
#include <vector>
//#include <set>

namespace gmml
{
namespace MolecularMetadata
{
namespace AMBER
{
struct AmberElement
{
    std::string element_;  // The element symbol
    double mass_;     // The mass, in amu, used by most/all of Amber
};

class AmberElementContainer
{
public:

    //////////////////////////////////////////////////////////
    //                       CONSTRUCTOR                    //
    //////////////////////////////////////////////////////////
    /*! \fn
            * Default constructor
            */
    AmberElementContainer(); // Calls an initializer?

    //////////////////////////////////////////////////////////
    //                         TYPEDEFS                     //
    //////////////////////////////////////////////////////////

    typedef std::vector<AmberElement> AmberElementVector;

    //////////////////////////////////////////////////////////
    //                      QUERY FUNCTIONS                 //
    //////////////////////////////////////////////////////////

    inline AmberElement GetEntryOfElement(std::string query)
    {
        AmberElement matching_entry;
        // can overload the == operator for the struct to compare type_, and do a oneliner std::find, or this:
        for (const auto& entry : amberElementVector_)
        {
            if (entry.element_.compare(query)==0)
            {
                matching_entry = entry;
            }
        }
        return matching_entry;
    }
    inline double GetMassOfElement(std::string query)
    {
        double return_mass;
        // can overload the == operator for the struct to compare type_, and do a oneliner std::find, or this:
        for (const auto& entry : amberElementVector_)
        {
            if (entry.element_.compare(query)==0)
            {
                return_mass = entry.mass_;
            }
        }
        return return_mass;
    }
private:
    AmberElementVector amberElementVector_;
};
}
}
}
#endif // AMBER_ELEMENTS_META_HPP
