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
    std::string atom1_ ;
    std::string atom2_ ;
    std::string atom3_ ;
    std::string atom4_ ;
    inline bool operator== (const DihedralAngleData &other) {return (this->index_ == other.index_ && this->number_of_bonds_from_anomeric_carbon_ == other.number_of_bonds_from_anomeric_carbon_);}
    inline bool operator!= (const DihedralAngleData &other) {return (this->index_ != other.index_ || this->number_of_bonds_from_anomeric_carbon_ != other.number_of_bonds_from_anomeric_carbon_);}
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
    DihedralAngleDataVector GetEntriesForLinkage( MolecularModeling::Atom* linking_atom1, MolecularModeling::Atom* linking_atom2);

private:
    //////////////////////////////////////////////////////////
    //                    PRIVATE FUNCTIONS                 //
    //////////////////////////////////////////////////////////
    // Some entries have conditions for the first or second residue to have a particular type (aka tag).
    // Most entries have "none" for condition. This checks first if condition is "none", and therefore satisfied.
    // Otherwise (else if) it checks if any of the residue_types match the condition for the entry, e.g. gauche_effect=galacto.
    bool checkIfResidueConditionsAreSatisfied(std::vector<std::string> residue_types, std::vector<std::string> entry_conditions);


    //////////////////////////////////////////////////////////
    //                       ATTRIBUTES                     //
    //////////////////////////////////////////////////////////
    DihedralAngleDataVector dihedralAngleDataVector_;
};
}
}
}
#endif // CARBOHYDRATE_DIHEDRAL_ANGLES_HPP
