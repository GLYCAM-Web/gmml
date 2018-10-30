#ifndef CARBOHYDRATE_DIHEDRAL_ANGLE_PREFERENCE_HPP
#define CARBOHYDRATE_DIHEDRAL_ANGLE_PREFERENCE_HPP

#include <string>
#include <vector>
#include <regex>
// Query function includes:
#include "../../MolecularModeling/atom.hpp"
#include "../../MolecularModeling/residue.hpp"
#include "glycam06residueinfo.hpp"

namespace gmml
{
namespace MolecularMetadata
{
namespace GLYCAM
{

struct DihedralAngleData
{
    std::string linking_atom1_ ;
    std::string linking_atom2_ ;
    std::string dihedral_angle_name_ ;
    double default_angle_value_ ;
    double lower_range_ ;
    double upper_range_ ;
    std::string name_ ;
    int index_ ; // if two entries match the criteria, and have the same index, the later entry should overwrite the earlier.
    std::string residue1_condition_ ;
    std::string residue2_condition_ ;
    std::string atom4_ ; // I always want rotation to be in the direction of protein->glycan and within glycan from reducing terminal to non-reducing.
    std::string atom3_ ;
    std::string atom2_ ;
    std::string atom1_ ;

} ;

class DihedralAngleDataContainer
{
public:

    //////////////////////////////////////////////////////////
    //                       CONSTRUCTOR                    //
    //////////////////////////////////////////////////////////
    /*! \fn
    * Default constructor
    */
    DihedralAngleDataContainer(); // Calls an initializer?

    //////////////////////////////////////////////////////////
    //                         TYPEDEFS                     //
    //////////////////////////////////////////////////////////

    typedef std::vector<DihedralAngleData> DihedralAngleDataVector;

    //////////////////////////////////////////////////////////
    //                      QUERY FUNCTIONS                 //
    //////////////////////////////////////////////////////////
    // Pass in the two atoms on either side the residue-residue linkage
    inline DihedralAngleDataVector GetEntriesForLinkage( MolecularModeling::Atom* linking_atom1, MolecularModeling::Atom* linking_atom2)
    {
        DihedralAngleDataVector matching_entries;
        Glycam06NamesToTypesLookupContainer metadata_residueNamesToTypes;
        // Go through each entry in the metadata
        for (const auto& entry : dihedralAngleDataVector_)
        {
            // Create a regex of each entry's linking_atom1_ and 2_. These are regex queries.
            std::regex regex1(entry.linking_atom1_, std::regex_constants::ECMAScript);
            std::regex regex2(entry.linking_atom2_, std::regex_constants::ECMAScript);
            // If metadata entry matches (regex query) to the two linking atom names
            if ( (std::regex_search(linking_atom1->GetName(), regex1)) && (std::regex_search(linking_atom2->GetName(), regex2)) )
            {
                // Some entries have conditions for the residue, that they have certain tags. Make sure any conditions are met:
                std::vector<std::string> residue1_types = metadata_residueNamesToTypes.GetTypesForResidue(linking_atom1->GetResidue()->GetName());
                std::vector<std::string> residue2_types = metadata_residueNamesToTypes.GetTypesForResidue(linking_atom2->GetResidue()->GetName());
                if ( (checkIfResidueConditionsAreSatisfied(residue1_types, entry.residue1_condition_))
                  && (checkIfResidueConditionsAreSatisfied(residue2_types, entry.residue2_condition_)) )
                {
                    matching_entries.push_back(entry);
                }
            }
        }
        return matching_entries;
    }
private:
    //////////////////////////////////////////////////////////
    //                    PRIVATE FUNCTIONS                 //
    //////////////////////////////////////////////////////////
    inline bool checkIfResidueConditionsAreSatisfied(std::vector<std::string> residue_types, std::string entry_condition)
    {
        bool conditionSatisfied = false;

        if (entry_condition.compare("none")==0)
        {
            conditionSatisfied = true;
        }
        else if (std::find(residue_types.begin(), residue_types.end(), entry_condition) != residue_types.end())
        {
            conditionSatisfied = true;
        }
        return conditionSatisfied;
    }

    //////////////////////////////////////////////////////////
    //                       ATTRIBUTES                     //
    //////////////////////////////////////////////////////////
    DihedralAngleDataVector dihedralAngleDataVector_;
};
}
}
}
#endif // CARBOHYDRATE_DIHEDRAL_ANGLES_HPP
