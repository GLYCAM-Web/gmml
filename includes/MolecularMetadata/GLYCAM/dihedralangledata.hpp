#ifndef CARBOHYDRATE_DIHEDRAL_ANGLE_PREFERENCE_HPP
#define CARBOHYDRATE_DIHEDRAL_ANGLE_PREFERENCE_HPP
#include "includes/MolecularMetadata/GLYCAM/glycam06residueinfo.hpp"
#include <string>
#include <vector>

namespace gmml
{
    namespace MolecularMetadata
    {
        namespace GLYCAM
        {
            struct DihedralAngleData
            {
                std::string linking_atom1_;
                std::string linking_atom2_;
                std::string dihedral_angle_name_;
                double default_angle_value_;
                double lower_deviation_;
                double upper_deviation_;
                double weight_;
                std::string rotamer_type_; // permutation or conformer
                std::string rotamer_name_;
                int number_of_bonds_from_anomeric_carbon_;
                int index_; // Used to indicate whether multiple entries are meant to overwrite each other or generate
                            // an additional angle
                std::vector<std::string> residue1_conditions_;
                std::vector<std::string> residue2_conditions_;
                std::string atom1_;
                std::string atom2_;
                std::string atom3_;
                std::string atom4_;

                inline std::string print()
                {
                    return atom1_ + "_" + atom2_ + "_" + atom3_ + "_" + atom4_ + " : " + rotamer_name_;
                }

                //////////////////////////////////////////////////////////
                //                       OPERATORS                      //
                //////////////////////////////////////////////////////////
                inline bool operator==(const DihedralAngleData& other)
                {
                    return (this->index_ == other.index_ &&
                            this->number_of_bonds_from_anomeric_carbon_ == other.number_of_bonds_from_anomeric_carbon_);
                }

                inline bool operator!=(const DihedralAngleData& other)
                {
                    return (this->index_ != other.index_ ||
                            this->number_of_bonds_from_anomeric_carbon_ != other.number_of_bonds_from_anomeric_carbon_);
                }
            };

            typedef std::vector<DihedralAngleData> DihedralAngleDataVector;

            class DihedralAngleDataContainer
            {
              public:
                //////////////////////////////////////////////////////////
                //                       CONSTRUCTOR                    //
                //////////////////////////////////////////////////////////
                DihedralAngleDataContainer();
                //////////////////////////////////////////////////////////
                //                      QUERY FUNCTIONS                 //
                //////////////////////////////////////////////////////////
                // Pass in the two atoms on either side the residue-residue linkage
                DihedralAngleDataVector GetEntriesForLinkage(const std::string atom1Name,
                                                             const std::string residue1Name,
                                                             const std::string atom2Name,
                                                             const std::string residue2Name) const;

              private:
                //////////////////////////////////////////////////////////
                //                    PRIVATE FUNCTIONS                 //
                //////////////////////////////////////////////////////////
                // Some entries have conditions for the first or second residue to have a particular type (aka tag).
                // Most entries have "none" for condition. This checks first if condition is "none", and therefore
                // satisfied. Otherwise (else if) it checks if any of the residue_types match the condition for the
                // entry, e.g. gauche_effect=galacto.
                bool checkIfResidueConditionsAreSatisfied(std::vector<std::string> residue_types,
                                                          std::vector<std::string> entry_conditions) const;
                //////////////////////////////////////////////////////////
                //                       ATTRIBUTES                     //
                //////////////////////////////////////////////////////////
                DihedralAngleDataVector dihedralAngleDataVector_;
            };
        } // namespace GLYCAM
    }     // namespace MolecularMetadata
} // namespace gmml
#endif // CARBOHYDRATE_DIHEDRAL_ANGLES_HPP
