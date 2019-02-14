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
    double lower_deviation_ ;
    double upper_deviation_ ;
    double weight_;
    std::string rotamer_type_ ; // permutation or conformer
    std::string rotamer_name_ ;
    int number_of_bonds_from_anomeric_carbon_;
    int index_ ; // Used to indicate whether multiple entries are meant to overwrite each other or generate an additional angle
    StringVector residue1_conditions_ ;
    StringVector residue2_conditions_ ;
//    std::string residue1_conditions_ ;
//    std::string residue2_conditions_ ;
    std::string atom1_ ;
    std::string atom2_ ;
    std::string atom3_ ;
    std::string atom4_ ;
} ;

typedef std::vector<DihedralAngleData> DihedralAngleDataVector;

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

   // typedef std::vector<DihedralAngleData> DihedralAngleDataVector;

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
            //std::cout << "Compare entry " << entry.linking_atom1_ << "-" << entry.linking_atom2_ << " : " << linking_atom1->GetName() << "-" << linking_atom2->GetName() <<"\n";
            std::regex regex1(entry.linking_atom1_, std::regex_constants::ECMAScript);
            std::regex regex2(entry.linking_atom2_, std::regex_constants::ECMAScript);
            // If metadata entry matches (regex query) to the two linking atom names
            if ( (std::regex_search(linking_atom1->GetName(), regex1)) && (std::regex_search(linking_atom2->GetName(), regex2)) )
            {
               // std::cout << "Checking for conditions: " << entry.residue1_condition_ << " + " << entry.residue2_condition_ << "\n";
                // Some entries have conditions for the residue, that they have certain tags. Make sure any conditions are met:
                std::vector<std::string> residue1_types = metadata_residueNamesToTypes.GetTypesForResidue(linking_atom1->GetResidue()->GetName());
                std::vector<std::string> residue2_types = metadata_residueNamesToTypes.GetTypesForResidue(linking_atom2->GetResidue()->GetName());
                if ( (checkIfResidueConditionsAreSatisfied(residue1_types, entry.residue1_conditions_))
                  && (checkIfResidueConditionsAreSatisfied(residue2_types, entry.residue2_conditions_)) )
                {
               //    std::cout << "Entry added: " << entry.linking_atom1_ << "-" << entry.linking_atom2_ << "\n";
                    matching_entries.push_back(entry);
                }
            }
//            std::string str("1231");
//            std::regex r("^(\\d)\\d"); // entire match will be 2 numbers
//            std::smatch m;
//            std::regex_search(str, m, r);
            //for(auto v: m) std::cout << v << std::endl;
//            std::smatch results1;
//            std::string string1 = linking_atom2->GetName();
//            bool well = std::regex_search(string1, results1, regex1);
//            std::cout << well << "\n";
//            for(auto v: results1) std::cout << v << std::endl;
        }
        // Not yet implemented: If two entries have same index number, delete the earlier entry.
        return matching_entries;
    }
private:
    //////////////////////////////////////////////////////////
    //                    PRIVATE FUNCTIONS                 //
    //////////////////////////////////////////////////////////
    // Some entries have conditions for the first or second residue to have a particular type (aka tag).
    // Most entries have "none" for condition. This checks first if condition is "none", and therefore satisfied.
    // Otherwise (else if) it checks if any of the residue_types match the condition for the entry, e.g. gauche_effect=galacto.
    inline bool checkIfResidueConditionsAreSatisfied(std::vector<std::string> residue_types, std::vector<std::string> entry_conditions)
    {
        for (const auto& entry_condition : entry_conditions)
        {
            if (entry_condition.compare("none")==0)
            {
                return true;
            }
            else if (!(std::find(residue_types.begin(), residue_types.end(), entry_condition) != residue_types.end()))
            {
                return false; //If any condition isn't satisified. return false.
            }
        }
        std::cout << "Logic error in dihedralangledata.hpp::checkIfResidueConditionsAreSatisfied" << std::endl;
    }

//    inline bool checkIfResidueConditionsAreSatisfied(std::vector<std::string> residue_types, std::string entry_condition)
//    {
//        bool conditionSatisfied = false;
//        if (entry_condition.compare("none")==0)
//        {
//            conditionSatisfied = true;
//        }
//        else if (std::find(residue_types.begin(), residue_types.end(), entry_condition) != residue_types.end())
//        {
//            conditionSatisfied = true;
//        }
//        return conditionSatisfied;
//    }

    //////////////////////////////////////////////////////////
    //                       ATTRIBUTES                     //
    //////////////////////////////////////////////////////////
    DihedralAngleDataVector dihedralAngleDataVector_;
};
}
}
}
#endif // CARBOHYDRATE_DIHEDRAL_ANGLES_HPP
