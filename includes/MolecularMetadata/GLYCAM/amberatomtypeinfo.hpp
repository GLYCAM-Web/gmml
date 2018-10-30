#ifndef AMBER_ATOM_TYPE_INFO_HPP
#define AMBER_ATOM_TYPE_INFO_HPP

#include <string>
#include <vector>
#include <set>

namespace gmml
{
namespace MolecularMetadata
{
namespace GLYCAM
{
// Move this to an amber space once the proper design is apparent
struct AmberAtomTypeInfo
{
    std::string type_ ;                 // The atom type
    std::string element_ ;              // The element symbol
    std::string amber_hybridization_ ;  // Hybridization expected by AMBER
    std::string hybridization_ ;        // Hybridization a chemist would specify
    double vdw_radius_ ;                // Default non-bonded atom radius
    std::string note_;                  // Note about this atom type
};

class AmberAtomTypeInfoContainer
{
public:

    //////////////////////////////////////////////////////////
    //                       CONSTRUCTOR                    //
    //////////////////////////////////////////////////////////
    /*! \fn
              * Default constructor
              */
    AmberAtomTypeInfoContainer();

    //////////////////////////////////////////////////////////
    //                         TYPEDEFS                     //
    //////////////////////////////////////////////////////////

    typedef std::vector<AmberAtomTypeInfo> AmberAtomTypeInfoVector;

    //////////////////////////////////////////////////////////
    //                      QUERY FUNCTIONS                 //
    //////////////////////////////////////////////////////////

    inline AmberAtomTypeInfo GetEntryWithAtomType(std::string query)
    {
        AmberAtomTypeInfo matching_entry;
        // can overload the == operator for the struct to compare type_, and do a oneliner std::find, or this:
        for (const auto& entry : amberAtomTypeInfoVector_)
        {
            if (entry.type_.compare(query)==0)
            {
                matching_entry = entry;
            }
        }
        return matching_entry;
    }
    inline std::string GetElementForAtomType(std::string query)
    {
        std::string element_of_matching_entry = "";
        // can overload the == operator for the struct to compare type_, and do a oneliner std::find, or this:
        for (const auto& entry : amberAtomTypeInfoVector_)
        {
            if (entry.type_.compare(query)==0)
            {
                element_of_matching_entry = entry.element_;
            }
        }
        return element_of_matching_entry;
    }
private:
    AmberAtomTypeInfoVector amberAtomTypeInfoVector_;
};
}
}
}
#endif // GLYCAM06META_HPP
