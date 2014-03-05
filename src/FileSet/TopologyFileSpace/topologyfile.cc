
#include "../../../includes/FileSet/TopologyFileSpace/topologyfile.hpp"
#include "../../../includes/FileSet/TopologyFileSpace/topologyatomtype.hpp"
#include "../../../includes/FileSet/TopologyFileSpace/topologybondtype.hpp"
#include "../../../includes/FileSet/TopologyFileSpace/topologyangletype.hpp"
#include "../../../includes/FileSet/TopologyFileSpace/topologydihedraltype.hpp"

using namespace std;
using namespace TopologyFileSpace;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
TopologyFile::TopologyFile() {}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
std::string TopologyFile::GetTitle()
{
    return title_;
}
int TopologyFile::GetNumberOfAtoms()
{
    return number_of_atoms_;
}
int TopologyFile::GetNumberOfTypes()
{
    return number_of_types_;
}
int TopologyFile::GetNumberOfBondsIncludingHydrogen()
{
    return number_of_bonds_including_hydrogen_;
}
int TopologyFile::GetNumberOfAnglesIncludingHydrogen()
{
    return number_of_angles_including_hydrogen_;
}
int TopologyFile::GetNumberOfDihedralsIncludingHydrogen()
{
    return number_of_dihedrals_including_hydrogen_;
}
int TopologyFile::GetNumberOfBondsExcludingHydrogen()
{
    return number_of_bonds_excluding_hydrogen_;
}
int TopologyFile::GetNumberOfAnglesExcludingHydrogen()
{
    return number_of_angles_excluding_hydrogen_;
}
int TopologyFile::GetNumberOfDihedralsExcludingHydrogen()
{
    return number_of_dihedrals_excluding_hydrogen_;
}
int TopologyFile::GetNumberOfHydrogenParameters()
{
    return number_of_hydrogen_parameters_;
}
int TopologyFile::GetNumberOfParameters()
{
    return number_of_parameters_;
}
int TopologyFile::GetNumberOfExcludedAtoms()
{
    return number_of_excluded_atoms_;
}
int TopologyFile::GetNumberOfResidues()
{
    return number_of_residues_;
}
int TopologyFile::GetTotalNumberOfBonds()
{
    return total_number_of_bonds_;
}
int TopologyFile::GetTotalNumberOfAngles()
{
    return total_number_of_angles_;
}
int TopologyFile::GetTotalNumberOfDihedrals()
{
    return total_number_of_dihedrals_;
}
int TopologyFile::GetNumberOfBondTypes()
{
    return number_of_bond_types_;
}
int TopologyFile::GetNumberOfAnglesTypes()
{
    return number_of_angle_types_;
}
int TopologyFile::GetNumberOfDihedralTypes()
{
    return number_of_dihedral_types_;
}
int TopologyFile::GetNumberOfAtomTypesInParameterFile()
{
    return number_of_atom_types_in_parameter_file_;
}
int TopologyFile::GetNumberOfDistinctHydrogenBonds()
{
    return number_of_distinct_hydrogen_bonds_;
}
bool TopologyFile::GetHasPerturbation()
{
    return has_perturbation_;
}
int TopologyFile::GetNumberOfBondsPerturbed()
{
    return number_of_bonds_perturbed_;
}
int TopologyFile::GetNumberOfAnglesPerturbed()
{
    return number_of_angles_perturbed_;
}
int TopologyFile::GetNumberOfDihedralsPerturbed()
{
    return number_of_dihedrals_perturbed_;
}
int TopologyFile::GetNumberOfBondsGroupPerturbed()
{
    return number_of_bonds_group_perturbed_;
}
int TopologyFile::GetNumberOfAnglesGroupPerturbed()
{
    return number_of_angles_group_perturbed_;
}
bool TopologyFile::GetHasStandardPeriodicBox()
{
    return has_standard_periodic_box_;
}
int TopologyFile::GetNumberOfAtomsInLargestResidue()
{
    return number_of_atoms_in_largest_residue_;
}
bool TopologyFile::GetHasCapOption()
{
    return has_cap_option_;
}
int TopologyFile::GetNumberOfExtraPoints()
{
    return number_of_extra_points_;
}
int TopologyFile::GetNumberOfBeads()
{
    return number_of_beads_;
}
TopologyFile::TopologyAtomTypeMap TopologyFile::GetAtomTypes()
{
    return atom_types_;
}
TopologyFile::TopologyBondTypeMap TopologyFile::GetBondTypes()
{
    return bond_types_;
}
TopologyFile::TopologyAngleTypeMap TopologyFile::GetAngleTypes()
{
    return angle_types_;
}
TopologyFile::TopologyDihedralTypeMap TopologyFile::GetDihedralTypes()
{
    return dihedral_types_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void TopologyFile::SetTitle(const string title)
{
    title_ = title;
}
void TopologyFile::SetNumberOfAtoms(int number_of_atoms)
{
    number_of_atoms_ = number_of_atoms;
}
void TopologyFile::SetNumberOfTypes(int number_of_types)
{
    number_of_types_ = number_of_types;
}
void TopologyFile::SetNumberOfBondsIncludingHydrogen(int number_of_bonds_including_hydrogen)
{
    number_of_bonds_including_hydrogen_ = number_of_bonds_including_hydrogen;
}
void TopologyFile::SetNumberOfAnglesIncludingHydrogen(int number_of_angles_including_hydrogen)
{
    number_of_angles_including_hydrogen_ = number_of_angles_including_hydrogen;
}
void TopologyFile::SetNumberOfDihedralsIncludingHydrogen(int number_of_dihedrals_including_hydrogen)
{
    number_of_dihedrals_including_hydrogen_ = number_of_dihedrals_including_hydrogen;
}
void TopologyFile::SetNumberOfBondsExcludingHydrogen(int number_of_bonds_excluding_hydrogen)
{
    number_of_bonds_excluding_hydrogen_ = number_of_bonds_excluding_hydrogen;
}
void TopologyFile::SetNumberOfAnglesExcludingHydrogen(int number_of_angles_excluding_hydrogen)
{
    number_of_angles_excluding_hydrogen_ = number_of_angles_excluding_hydrogen;
}
void TopologyFile::SetNumberOfDihedralsExcludingHydrogen(int number_of_dihedrals_excluding_hydrogen)
{
    number_of_dihedrals_excluding_hydrogen_ = number_of_dihedrals_excluding_hydrogen;
}
void TopologyFile::SetNumberOfHydrogenParameters(int number_of_hydrogen_parameters)
{
    number_of_hydrogen_parameters_ = number_of_hydrogen_parameters;
}
void TopologyFile::SetNumberOfParameters(int number_of_parameters)
{
    number_of_parameters_ = number_of_parameters;
}
void TopologyFile::SetNumberOfExcludedAtoms(int number_of_excluded_atoms)
{
    number_of_excluded_atoms_ = number_of_excluded_atoms;
}
void TopologyFile::SetNumberOfResidues(int number_of_residues)
{
    number_of_residues_ = number_of_residues;
}
void TopologyFile::SetTotalNumberOfBonds(int total_number_of_bonds)
{
    total_number_of_bonds_ = total_number_of_bonds;
}
void TopologyFile::SetTotalNumberOfAngles(int total_number_of_angles)
{
    total_number_of_angles_ = total_number_of_angles;
}
void TopologyFile::SetTotalNumberOfDihedrals(int total_number_of_dihedrals)
{
    total_number_of_dihedrals_ = total_number_of_dihedrals;
}
void TopologyFile::SetNumberOfBondTypes(int number_of_bond_types)
{
    number_of_bond_types_ = number_of_bond_types;
}
void TopologyFile::SetNumberOfAngleTypes(int number_of_angle_types)
{
    number_of_angle_types_ = number_of_angle_types;
}
void TopologyFile::SetNumberOfDihedralTypes(int number_of_dihedral_types)
{
    number_of_dihedral_types_ = number_of_dihedral_types;
}
void TopologyFile::SetNumberOfAtomTypesInParameterFile(int number_of_atom_types_in_parameter_file)
{
    number_of_atom_types_in_parameter_file_ = number_of_atom_types_in_parameter_file;
}
void TopologyFile::SetNumberOfDistinctHydrogenBonds(int number_of_distinct_hydrogen_bonds)
{
    number_of_distinct_hydrogen_bonds_ = number_of_distinct_hydrogen_bonds;
}
void TopologyFile::SetHasPerturbation(bool has_perturbation)
{
    has_perturbation_ = has_perturbation;
}
void TopologyFile::SetNumberOfBondsPerturbed(int number_of_bonds_perturbed)
{
    number_of_bonds_perturbed_ = number_of_bonds_perturbed;
}
void TopologyFile::SetNumberOfAnglesPerturbed(int number_of_angles_perturbed)
{
    number_of_angles_perturbed_ = number_of_angles_perturbed;
}
void TopologyFile::SetNumberOfDihedralsPerturbed(int number_of_dihedrals_perturbed)
{
    number_of_dihedrals_perturbed_ = number_of_dihedrals_perturbed;
}
void TopologyFile::SetNumberOfBondsGroupPerturbed(int number_of_bonds_group_perturbed)
{
    number_of_bonds_group_perturbed_ = number_of_bonds_group_perturbed;
}
void TopologyFile::SetNumberOfAnglesGroupPerturbed(int number_of_angles_group_perturbed)
{
    number_of_angles_group_perturbed_ = number_of_angles_group_perturbed;
}
void TopologyFile::SetHasStandardPeriodicBox(bool has_standard_periodic_box)
{
    has_standard_periodic_box_ = has_standard_periodic_box;
}
void TopologyFile::SetNumberOfAtomsInLargestResidue(int number_of_atoms_in_largest_residue)
{
    number_of_atoms_in_largest_residue_ = number_of_atoms_in_largest_residue;
}
void TopologyFile::SetHasCapOption(bool has_cap_option)
{
    has_cap_option_ = has_cap_option;
}
void TopologyFile::SetNumberOfExtraPoints(int number_of_extra_points)
{
    number_of_extra_points_ = number_of_extra_points;
}
void TopologyFile::SetNumberOfBeads(int number_of_beads)
{
    number_of_beads_ = number_of_beads;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void TopologyFile::Print(ostream &out)
{
}




