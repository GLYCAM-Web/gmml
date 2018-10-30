#ifndef BOND_LENGTH_BY_TYPE_PAIR_HPP
#define BOND_LENGTH_BY_TYPE_PAIR_HPP

#include <string>
#include <vector>
#include "../../common.hpp"
//#include "common.hpp"
namespace gmml
{
namespace MolecularMetadata
{
namespace GLYCAM
{
struct BondLengthByTypePair {
    std::string type1_;  // One of the atom types
    std::string type2_;  // The other atom type
    double length_;      // Default length in Angstroms
    std::string note_;   // Note about this bond.
} ;

class BondLengthByTypePairContainer
{
public:

    //////////////////////////////////////////////////////////
    //                       CONSTRUCTOR                    //
    //////////////////////////////////////////////////////////
    /*! \fn
            * Default constructor
            */
    BondLengthByTypePairContainer(); // Calls an initializer?

    //////////////////////////////////////////////////////////
    //                         TYPEDEFS                     //
    //////////////////////////////////////////////////////////

    typedef std::vector<BondLengthByTypePair> BondLengthByTypePairVector;

    //////////////////////////////////////////////////////////
    //                      QUERY FUNCTIONS                 //
    //////////////////////////////////////////////////////////

    inline double GetBondLengthForAtomTypes(std::string query1, std::string query2)
    {
        double matching_entry_bond_length = gmml::BOND_LENGTH;
        // can overload the == operator for the struct to compare type_, and do a oneliner std::find, or this:
        for (const auto& entry : bondLengthByTypePairVector_)
        {   //Search bidirectionally e.g Cg-Os, Os-Cg
            if ( ( (entry.type1_.compare(query1)==0) && (entry.type2_.compare(query2)==0) ) ||
                 ( (entry.type1_.compare(query2)==0) && (entry.type2_.compare(query1)==0) )
                 )
            {
                matching_entry_bond_length = entry.length_;
            }
        }
        return matching_entry_bond_length;
    }
private:
    BondLengthByTypePairVector bondLengthByTypePairVector_;
};


}
}
}
#endif // BOND_LENGTH_BY_TYPE_PAIR_HPP
