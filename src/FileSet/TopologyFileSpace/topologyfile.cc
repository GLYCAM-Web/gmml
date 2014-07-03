#include <fstream>
#include <math.h>

#include "../../../includes/FileSet/TopologyFileSpace/topologyfile.hpp"
#include "../../../includes/FileSet/TopologyFileSpace/topologyatompair.hpp"
#include "../../../includes/FileSet/TopologyFileSpace/topologybondtype.hpp"
#include "../../../includes/FileSet/TopologyFileSpace/topologyangletype.hpp"
#include "../../../includes/FileSet/TopologyFileSpace/topologydihedraltype.hpp"
#include "../../../includes/FileSet/TopologyFileSpace/topologyassembly.hpp"
#include "../../../includes/FileSet/TopologyFileSpace/topologyatom.hpp"
#include "../../../includes/FileSet/TopologyFileSpace/topologybond.hpp"
#include "../../../includes/FileSet/TopologyFileSpace/topologyangle.hpp"
#include "../../../includes/FileSet/TopologyFileSpace/topologydihedral.hpp"
#include "../../../includes/FileSet/TopologyFileSpace/topologyresidue.hpp"
#include "../../../includes/FileSet/TopologyFileSpace/topologybond.hpp"
#include "../../../includes/FileSet/TopologyFileSpace/topologyangle.hpp"
#include "../../../includes/FileSet/TopologyFileSpace/topologydihedral.hpp"
#include "../../../includes/FileSet/TopologyFileSpace/topologyfileprocessingexception.hpp"
#include "../../../includes/utils.hpp"

using namespace std;
using namespace TopologyFileSpace;
using namespace gmml;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
TopologyFile::TopologyFile() {}

TopologyFile::TopologyFile(const string &top_file)
{
    ifstream in_file;
    try
    {
        in_file.open(top_file.c_str());
    }
    catch(exception &ex)
    {
        throw TopologyFileProcessingException(__LINE__, "File not found");
    }
    Read(in_file);
    in_file.close();            /// Close the parameter files
}

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
int TopologyFile::GetPerturbationOption()
{
    return perturbation_option_;
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
int TopologyFile::GetNumberOfDihedralsGroupPerturbed()
{
    return number_of_dihedrals_group_perturbed_;
}
int TopologyFile::GetStandardPeriodicBoxOption()
{
    return standard_periodic_box_option_;
}
int TopologyFile::GetNumberOfAtomsInLargestResidue()
{
    return number_of_atoms_in_largest_residue_;
}
int TopologyFile::GetCapOption()
{
    return cap_option_;
}
int TopologyFile::GetNumberOfExtraPoints()
{
    return number_of_extra_points_;
}
int TopologyFile::GetNumberOfBeads()
{
    return number_of_beads_;
}
TopologyFile::TopologyAtomPairMap TopologyFile::GetPairs()
{
    return pairs_;
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
TopologyAssembly* TopologyFile::GetAssembly()
{
    return assembly_;
}

TopologyFile::RadiusSet TopologyFile::GetRadiusSet()
{
    return radius_set_;
}
TopologyFile::TopologyBondMap TopologyFile::GetBonds()
{
    return bonds_;
}
TopologyFile::TopologyAngleMap TopologyFile::GetAngles()
{
    return angles_;
}
TopologyFile::TopologyDihedralMap TopologyFile::GetDihedrals()
{
    return dihedrals_;
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
void TopologyFile::SetPerturbationOption(int perturbation_option)
{
    perturbation_option_ = perturbation_option;
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
void TopologyFile::SetNumberOfDihedralsGroupPerturbed(int number_of_dihedrals_group_perturbed)
{
    number_of_dihedrals_group_perturbed_ = number_of_dihedrals_group_perturbed;
}
void TopologyFile::SetStandardPeriodicBoxOption(int standard_periodic_box_option)
{
    standard_periodic_box_option_ = standard_periodic_box_option;
}
void TopologyFile::SetNumberOfAtomsInLargestResidue(int number_of_atoms_in_largest_residue)
{
    number_of_atoms_in_largest_residue_ = number_of_atoms_in_largest_residue;
}
void TopologyFile::SetCapOption(int cap_option)
{
    cap_option_ = cap_option;
}
void TopologyFile::SetNumberOfExtraPoints(int number_of_extra_points)
{
    number_of_extra_points_ = number_of_extra_points;
}
void TopologyFile::SetNumberOfBeads(int number_of_beads)
{
    number_of_beads_ = number_of_beads;
}
void TopologyFile::SetAssembly(TopologyAssembly *assembly)
{
    assembly_ = assembly;
}

void TopologyFile::SetRadiusSet(RadiusSet radius_set)
{
    radius_set_.clear();
    for(RadiusSet::iterator it = radius_set.begin(); it != radius_set.end(); it++)
    {
        radius_set_.push_back(*it);
    }
}
void TopologyFile::SetBonds(TopologyBondMap bonds)
{
    bonds_.clear();
    for(TopologyBondMap::iterator it = bonds.begin(); it != bonds.end(); it++)
    {
        vector<string> atom_bond = (*it).first;
        TopologyBond* bond = (*it).second;
        bonds_[atom_bond] = bond;
    }
}
void TopologyFile::SetAngles(TopologyAngleMap angles)
{
    angles_.clear();
    for(TopologyAngleMap::iterator it = angles.begin(); it != angles.end(); it++)
    {
        vector<string> atom_angle = (*it).first;
        TopologyAngle* angle = (*it).second;
        angles_[atom_angle] = angle;
    }
}
void TopologyFile::SetDihedrals(TopologyDihedralMap dihedrals)
{
    dihedrals_.clear();
    for(TopologyDihedralMap::iterator it = dihedrals.begin(); it != dihedrals.end(); it++)
    {
        vector<string> atom_dihedral = (*it).first;
        TopologyDihedral* dihedral = (*it).second;
        dihedrals_[atom_dihedral] = dihedral;
    }
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////
void TopologyFile::Read(ifstream &in_file)
{
    this->ParseSections(in_file);
}

void TopologyFile::ParseSections(ifstream &in_stream)
{
    string line;
    /// Unable to read file
    if (!getline(in_stream, line))
    {
        throw TopologyFileProcessingException("Error reading file");
    }
    stringstream other;
    vector<string> atom_names = vector<string>();
    vector<double> charges = vector<double>();
    vector<int> atomic_numbers = vector<int>();
    vector<double> masses = vector<double>();
    vector<int> atom_type_indexes = vector<int>();
    vector<int> number_excluded_atoms = vector<int>();
    vector<int> nonbonded_parm_indexes = vector<int>();
    vector<string> residue_labels = vector<string>();
    vector<int> residue_pointers = vector<int>();
    vector<double> bond_force_constants = vector<double>();
    vector<double> bond_equil_values = vector<double>();
    vector<double> angle_force_constants = vector<double>();
    vector<double> angle_equil_values = vector<double>();
    vector<double> dihedral_force_constants = vector<double>();
    vector<double> dihedral_periodicities = vector<double>();
    vector<double> dihedral_phases = vector<double>();
    vector<double> scee_scale_factors = vector<double>();
    vector<double> scnb_scale_factors = vector<double>();
    vector<double> lennard_jones_acoefs = vector<double>();
    vector<double> lennard_jones_bcoefs = vector<double>();
    vector<int> bonds_inc_hydrogens = vector<int>();
    vector<int> bonds_without_hydrogens = vector<int>();
    vector<int> angles_inc_hydrogens = vector<int>();
    vector<int> angles_without_hydrogens = vector<int>();
    vector<int> dihedrals_inc_hydrogens = vector<int>();
    vector<int> dihedrals_without_hydrogens = vector<int>();
    vector<int> excluded_atoms_lists = vector<int>();
    vector<double> hbond_acoefs = vector<double>();
    vector<double> hbond_bcoefs = vector<double>();
    vector<double> hb_cuts = vector<double>();
    vector<string> amber_atom_types = vector<string>();
    vector<string> tree_chain_classifications = vector<string>();
    vector<string> radius_sets = vector<string>();
    vector<double> radiis = vector<double>();
    vector<double> screens = vector<double>();
    vector<int> ipols = vector<int>();
    while(!line.empty())
    {
        if(line.find("%FLAG") != string::npos)
        {
            stringstream section;
            PartitionSection(in_stream, line, section);
            string in_line;
//            cout << section.str() << endl;
            getline(section, in_line);
            if(in_line.find("%FLAG TITLE") != string::npos)
            {
                ParseTitlePartition(section);
            }
            else if(in_line.find("%FLAG POINTERS") != string::npos)
            {
                ParsePointersPartition(section);
            }
            else if(in_line.find("%FLAG ATOM_NAME") != string::npos)
            {
                atom_names = ParsePartition<string>(section);
            }
            else if(in_line.find("%FLAG CHARGE") != string::npos)
            {
                charges = ParsePartition<double>(section);
            }
            else if(in_line.find("%FLAG ATOMIC_NUMBER") != string::npos)
            {
                atomic_numbers = ParsePartition<int>(section);
            }
            else if(in_line.find("%FLAG MASS") != string::npos)
            {
                masses = ParsePartition<double>(section);
            }
            else if(in_line.find("%FLAG ATOM_TYPE_INDEX") != string::npos)
            {
                atom_type_indexes = ParsePartition<int>(section);
            }
            else if(in_line.find("%FLAG NUMBER_EXCLUDED_ATOMS") != string::npos)
            {
                number_excluded_atoms =ParsePartition<int>(section);
            }
            else if(in_line.find("%FLAG NONBONDED_PARM_INDEX") != string::npos)
            {
                nonbonded_parm_indexes = ParsePartition<int>(section);
            }
            else if(in_line.find("%FLAG RESIDUE_LABEL") != string::npos)
            {
                residue_labels = ParsePartition<string>(section);
            }
            else if(in_line.find("%FLAG RESIDUE_POINTER") != string::npos)
            {
                residue_pointers = ParsePartition<int>(section);
            }
            else if(in_line.find("%FLAG BOND_FORCE_CONSTANT") != string::npos)
            {
                bond_force_constants = ParsePartition<double>(section);
            }
            else if(in_line.find("%FLAG BOND_EQUIL_VALUE") != string::npos)
            {
                bond_equil_values = ParsePartition<double>(section);
            }
            else if(in_line.find("%FLAG ANGLE_FORCE_CONSTANT") != string::npos)
            {
                angle_force_constants = ParsePartition<double>(section);
            }
            else if(in_line.find("%FLAG ANGLE_EQUIL_VALUE") != string::npos)
            {
                angle_equil_values = ParsePartition<double>(section);
            }
            else if(in_line.find("%FLAG DIHEDRAL_FORCE_CONSTANT") != string::npos)
            {
                dihedral_force_constants = ParsePartition<double>(section);
            }
            else if(in_line.find("%FLAG DIHEDRAL_PERIODICITY") != string::npos)
            {
                dihedral_periodicities = ParsePartition<double>(section);
            }
            else if(in_line.find("%FLAG DIHEDRAL_PHASE") != string::npos)
            {
                dihedral_phases = ParsePartition<double>(section);
            }
            else if(in_line.find("%FLAG SCEE_SCALE_FACTOR") != string::npos)
            {
                scee_scale_factors = ParsePartition<double>(section);
            }
            else if(in_line.find("%FLAG SCNB_SCALE_FACTOR") != string::npos)
            {
                scnb_scale_factors = ParsePartition<double>(section);
            }
            else if(in_line.find("%FLAG SOLTY") != string::npos)
            {
            }
            else if(in_line.find("%FLAG LENNARD_JONES_ACOEF") != string::npos)
            {
                lennard_jones_acoefs = ParsePartition<double>(section);
            }
            else if(in_line.find("%FLAG LENNARD_JONES_BCOEF") != string::npos)
            {
                lennard_jones_bcoefs = ParsePartition<double>(section);
            }
            else if(in_line.find("%FLAG BONDS_INC_HYDROGEN") != string::npos)
            {
                bonds_inc_hydrogens = ParsePartition<int>(section);
            }
            else if(in_line.find("%FLAG BONDS_WITHOUT_HYDROGEN") != string::npos)
            {
                bonds_without_hydrogens = ParsePartition<int>(section);
            }
            else if(in_line.find("%FLAG ANGLES_INC_HYDROGEN") != string::npos)
            {
                angles_inc_hydrogens = ParsePartition<int>(section);
            }
            else if(in_line.find("%FLAG ANGLES_WITHOUT_HYDROGEN") != string::npos)
            {
                angles_without_hydrogens = ParsePartition<int>(section);
            }
            else if(in_line.find("%FLAG DIHEDRALS_INC_HYDROGEN") != string::npos)
            {
                dihedrals_inc_hydrogens = ParsePartition<int>(section);
            }
            else if(in_line.find("%FLAG DIHEDRALS_WITHOUT_HYDROGEN") != string::npos)
            {
                dihedrals_without_hydrogens = ParsePartition<int>(section);
            }
            else if(in_line.find("%FLAG EXCLUDED_ATOMS_LIST") != string::npos)
            {
                excluded_atoms_lists = ParsePartition<int>(section);
            }
            else if(in_line.find("%FLAG HBOND_ACOEF") != string::npos)
            {
                hbond_acoefs = ParsePartition<double>(section);
            }
            else if(in_line.find("%FLAG HBOND_BCOEF") != string::npos)
            {
                hbond_bcoefs = ParsePartition<double>(section);
            }
            else if(in_line.find("%FLAG HBCUT") != string::npos)
            {
                hb_cuts = ParsePartition<double>(section);
            }
            else if(in_line.find("%FLAG AMBER_ATOM_TYPE") != string::npos)
            {
                amber_atom_types = ParsePartition<string>(section);
            }
            else if(in_line.find("%FLAG TREE_CHAIN_CLASSIFICATION") != string::npos)
            {
                tree_chain_classifications = ParsePartition<string>(section);
            }
            else if(in_line.find("%FLAG JOIN_ARRAY") != string::npos)
            {
            }
            else if(in_line.find("%FLAG IROTAT") != string::npos)
            {
            }
            else if(in_line.find("%FLAG RADIUS_SET") != string::npos)
            {
                radius_sets = ParsePartition<string>(section);
            }
            else if(in_line.find("%FLAG RADII") != string::npos)
            {
                radiis = ParsePartition<double>(section);
            }
            else if(in_line.find("%FLAG SCREEN") != string::npos)
            {
                screens= ParsePartition<double>(section);
            }
            else if(in_line.find("%FLAG IPOL") != string::npos)
            {
                ipols = ParsePartition<int>(section);
            }
            else if(in_line.find("%FLAG SOLVENT_POINTERS") != string::npos)
            {
            }
            else if(in_line.find("%FLAG ATOMS_PER_MOLECULE") != string::npos)
            {
            }
            else if(in_line.find("%FLAG BOX_DIMENSIONS") != string::npos)
            {
            }
            else if(in_line.find("%FLAG CAP_INFO") != string::npos)
            {
            }
            else if(in_line.find("%FLAG CAP_INFO2") != string::npos)
            {
            }
            else if(in_line.find("%FLAG SOLVENT_POINTERS") != string::npos)
            {
            }
            else
            {
                // Ignore section
            }
        }
        else
        {
            other << line << endl;
            getline(in_stream, line);
        }
    }

    // Bond types
    for(int i = 0; i < number_of_bond_types_; i++)
    {
        bond_types_[i] = new TopologyBondType(i, bond_force_constants.at(i), bond_equil_values.at(i));
    }
    // Angle types
    for(int i = 0; i < number_of_angle_types_; i++)
    {
        angle_types_[i] = new TopologyAngleType(i, angle_force_constants.at(i), angle_equil_values.at(i));
    }
    // Dihedral types
    for(int i = 0; i < number_of_dihedral_types_; i++)
    {
        dihedral_types_[i] = new TopologyDihedralType(i, dihedral_force_constants.at(i), dihedral_periodicities.at(i), dihedral_phases.at(i),
                                                      scee_scale_factors.at(i), scnb_scale_factors.at(i));
    }

    // Radius set
    for(RadiusSet::iterator it = radius_sets.begin(); it != radius_sets.end(); it++)
    {
        radius_set_.push_back(*it);
    }
    // Lennard Jones coefficients for atom pairs
    for(int i = 0; i < number_of_atoms_; i++)
    {
        string atom_type_1 = amber_atom_types.at(i);
        TopologyAtomPair::TopologyCoefficientMap coefficient_a_map;
        TopologyAtomPair::TopologyCoefficientMap coefficient_b_map;
        for(int j = 0; j < number_of_atoms_; j++)
        {
            string atom_type_2 = amber_atom_types.at(j);
            double coefficient_a;
            double coefficient_b;
            int index = nonbonded_parm_indexes.at(number_of_types_ * (atom_type_indexes.at(i) - 1) + atom_type_indexes.at(j) - 1);
            if(index > 0)
            {
                coefficient_a = lennard_jones_acoefs.at(index - 1);
                coefficient_b = lennard_jones_bcoefs.at(index - 1);
            }
            else
            {
                coefficient_a = hbond_acoefs.at(index - 1);
                coefficient_b = hbond_bcoefs.at(index - 1);
            }
            coefficient_a_map[atom_type_2] = coefficient_a;
            coefficient_b_map[atom_type_2] = coefficient_b;
        }
        pairs_[atom_type_1] = new TopologyAtomPair(atom_type_1, coefficient_a_map, coefficient_b_map);

    }
    // Bonds, Angles, Dihedrals
    //Bonds in topology file
    for(int i = 0; i < number_of_bonds_including_hydrogen_; i++)
    {
        vector<string> bonds = vector<string>();

        int first_atom_index = (bonds_inc_hydrogens[i*3])/3;
        int second_atom_index = (bonds_inc_hydrogens[i*3+1])/3;

        bonds.push_back(atom_names[first_atom_index]);
        bonds.push_back(atom_names[second_atom_index]);

        TopologyBond* topology_bond = new TopologyBond(bonds);

        topology_bond->SetBondType(bond_types_[bonds_inc_hydrogens[i*3+2] - 1]);
        topology_bond->SetIncludingHydrogen(true);

        bonds_[bonds] = topology_bond;
    }
    for(int i = 0; i < number_of_bonds_excluding_hydrogen_; i++)
    {
        vector<string> bonds = vector<string>();

        int first_atom_index = (bonds_without_hydrogens[i*3])/3;
        int second_atom_index = (bonds_without_hydrogens[i*3+1])/3;

        bonds.push_back(atom_names[first_atom_index]);
        bonds.push_back(atom_names[second_atom_index]);

        TopologyBond* topology_bond = new TopologyBond(bonds);

        topology_bond->SetBondType(bond_types_[bonds_without_hydrogens[i*3+2] - 1]);
        topology_bond->SetIncludingHydrogen(false);

        bonds_[bonds] = topology_bond;
    }
    // Angles in topology file
    // Angles including hydrogen
    for(int i = 0; i < number_of_angles_including_hydrogen_; i++)
    {
        vector<string> angle_atoms = vector<string>();

        int atom_index_1 = angles_inc_hydrogens.at(i*4) / 3;
        int atom_index_2 = angles_inc_hydrogens.at(i*4+1) / 3;
        int atom_index_3 = angles_inc_hydrogens.at(i*4+2) / 3;

        angle_atoms.push_back(atom_names.at(atom_index_1));
        angle_atoms.push_back(atom_names.at(atom_index_2));
        angle_atoms.push_back(atom_names.at(atom_index_3));

        TopologyAngle* angle = new TopologyAngle(angle_atoms);

        angle->SetIncludingHydrogen(true);
        angle->SetAnlgeType(angle_types_[angles_inc_hydrogens.at(i*4+3) - 1]);

        angles_[angle_atoms] = angle;
    }
    // Angles excluding hydrogen
    for(int i = 0; i < number_of_angles_excluding_hydrogen_; i++)
    {
        vector<string> angle_atoms = vector<string>();

        int atom_index_1 = angles_without_hydrogens.at(i*4) / 3;
        int atom_index_2 = angles_without_hydrogens.at(i*4+1) / 3;
        int atom_index_3 = angles_without_hydrogens.at(i*4+2) / 3;

        angle_atoms.push_back(atom_names.at(atom_index_1));
        angle_atoms.push_back(atom_names.at(atom_index_2));
        angle_atoms.push_back(atom_names.at(atom_index_3));

        TopologyAngle* angle = new TopologyAngle(angle_atoms);

        angle->SetIncludingHydrogen(false);
        angle->SetAnlgeType(angle_types_[angles_without_hydrogens.at(i*4+3) - 1]);

        angles_[angle_atoms] = angle;
    }
    // Dihedrals in topology file
    // Dihedrals including hydrogen
    for(int i = 0; i < number_of_dihedrals_including_hydrogen_; i++)
    {
        vector<string> dihedral_atoms = vector<string>();

        int atom_index_1 = abs(dihedrals_inc_hydrogens.at(i*5)) / 3;
        int atom_index_2 = abs(dihedrals_inc_hydrogens.at(i*5+1)) / 3;
        int atom_index_3 = abs(dihedrals_inc_hydrogens.at(i*5+2)) / 3;
        int atom_index_4 = abs(dihedrals_inc_hydrogens.at(i*5+3)) / 3;

        dihedral_atoms.push_back(atom_names.at(atom_index_1));
        dihedral_atoms.push_back(atom_names.at(atom_index_2));
        dihedral_atoms.push_back(atom_names.at(atom_index_3));
        dihedral_atoms.push_back(atom_names.at(atom_index_4));

        TopologyDihedral* dihedral = new TopologyDihedral(dihedral_atoms);

        if(dihedrals_inc_hydrogens.at(i*5+2) < 0)
            dihedral->SetIgnoredGroupInteraction(true);
        else
            dihedral->SetIgnoredGroupInteraction(false);
        if(dihedrals_inc_hydrogens.at(i*5+3) < 0)
            dihedral->SetIsImproper(true);
        else
            dihedral->SetIsImproper(false);

        dihedral->SetIncludingHydrogen(true);
        dihedral->SetDihedralType(dihedral_types_[dihedrals_inc_hydrogens.at(i*5+4) - 1]);
        dihedrals_[dihedral_atoms] = dihedral;

    }
    // Dihedrals excluding hydrogen
    for(int i = 0; i < number_of_dihedrals_excluding_hydrogen_; i++)
    {
        vector<string> dihedral_atoms = vector<string>();

        int atom_index_1 = abs(dihedrals_without_hydrogens.at(i*5)) / 3;
        int atom_index_2 = abs(dihedrals_without_hydrogens.at(i*5+1)) / 3;
        int atom_index_3 = abs(dihedrals_without_hydrogens.at(i*5+2)) / 3;
        int atom_index_4 = abs(dihedrals_without_hydrogens.at(i*5+3)) / 3;

        dihedral_atoms.push_back(atom_names.at(atom_index_1));
        dihedral_atoms.push_back(atom_names.at(atom_index_2));
        dihedral_atoms.push_back(atom_names.at(atom_index_3));
        dihedral_atoms.push_back(atom_names.at(atom_index_4));

        TopologyDihedral* dihedral = new TopologyDihedral(dihedral_atoms);

        if(dihedrals_without_hydrogens.at(i*5+2) < 0)
            dihedral->SetIgnoredGroupInteraction(true);
        else
            dihedral->SetIgnoredGroupInteraction(false);
        if(dihedrals_without_hydrogens.at(i*5+3) < 0)
            dihedral->SetIsImproper(true);
        else
            dihedral->SetIsImproper(false);

        dihedral->SetIncludingHydrogen(false);
        dihedral->SetDihedralType(dihedral_types_[dihedrals_without_hydrogens.at(i*5+4) - 1]);

        dihedrals_[dihedral_atoms] = dihedral;
    }
    // Residues in topology file
    int start_index = 0;
    TopologyAssembly::TopologyResidueMap residues;
    for(vector<string>::iterator it = residue_labels.begin(); it != residue_labels.end(); it++)
    {
        string residue_name = *it;
        TopologyResidue::TopologyAtomMap atoms;
        int residue_index = distance(residue_labels.begin(), it) + 1;
        int starting_atom_index = residue_pointers.at(residue_index - 1);
        int ending_atom_index;
        if(residue_index < number_of_residues_)
        {
            ending_atom_index = residue_pointers.at(residue_index);
        }
        else
        {
            ending_atom_index = number_of_atoms_;
        }
        for(int i = starting_atom_index - 1; i < ending_atom_index - 1; i++)
        {
            TopologyAtom::ExcludedAtomNames excluded_atoms = TopologyAtom::ExcludedAtomNames();
            if(i > 0)
                start_index += number_excluded_atoms.at(i-1);
            for(int j = start_index; j < start_index + number_excluded_atoms.at(i); j++)
            {
                cout << excluded_atoms_lists.at(j) - 1 << endl;
                excluded_atoms.push_back(atom_names.at((excluded_atoms_lists.at(j) - 1 == -1) ? 0 : excluded_atoms_lists.at(j) - 1));
            }
            atoms[atom_names.at(i)] = new TopologyAtom(i + 1, atom_names.at(i), amber_atom_types.at(i), charges.at(i), atomic_numbers.at(i), masses.at(i), excluded_atoms,
                                                       number_excluded_atoms.at(i), radiis.at(i), screens.at(i), tree_chain_classifications.at(i));
        }
        residues[residue_name] = new TopologyResidue(residue_name, atoms, residue_index, starting_atom_index);
    }
    assembly_ = new TopologyAssembly();
    assembly_->SetAssemblyName(title_);
    assembly_->SetResidues(residues);
}

void TopologyFile::PartitionSection(ifstream &stream, string &line, stringstream& section)
{
    while(line[0] == '%')
    {
        if(line.find("%FLAG") != string::npos)
        {
            section << line << endl;
        }
        if(line.find("%FORMAT") != string::npos)
        {
            section << line << endl;
        }
        getline(stream, line);
    }
    while(line[0] != '%')
    {
        section << line << endl;
        getline(stream, line);
        if(line.empty())
            break;
    }
}

void TopologyFile::ParseTitlePartition(stringstream& stream)
{
    string line;
    // Section Header
    getline(stream, line);
    // Section Format
    getline(stream, line);
    while(!line.empty())
    {
        string format = "20a4";
        if(line[0] == '%')
        {
            if(line.find("%FORMAT") != string::npos)
            {
                format = line.substr(line.find_first_of('(', 0) + 1, line.find_last_of(')', line.length()) - line.find_first_of('(', 0) - 1);
            }
        }
        else
        {
            title_ += line;
        }
        getline(stream, line);
    }
    title_ = Trim(title_);
}

void TopologyFile::ParsePointersPartition(stringstream &stream)
{
    string line;
    vector<int> items = vector<int>();
    // Section Header
    getline(stream, line);
    // Section Format
    getline(stream, line);
    while(!line.empty())
    {
        string format = "10I8";
        if(line[0] == '%')
        {
            if(line.find("%FORMAT") != string::npos)
            {
                format = line.substr(line.find_first_of('(', 0) + 1, line.find_last_of(')', line.length()) - line.find_first_of('(', 0) - 1);
            }
        }
        else
        {
            vector<int> line_items = PartitionLine<int>(line, format);
            for(vector<int>::iterator it = line_items.begin(); it != line_items.end(); it++)
            {
                items.push_back(*it);
            }
        }
        getline(stream, line);
    }
    number_of_atoms_ = items.at(0);
    number_of_types_ = items.at(1);
    number_of_bonds_including_hydrogen_ = items.at(2);
    number_of_bonds_excluding_hydrogen_ = items.at(3);
    number_of_angles_including_hydrogen_ = items.at(4);
    number_of_angles_excluding_hydrogen_ = items.at(5);
    number_of_dihedrals_including_hydrogen_ = items.at(6);
    number_of_dihedrals_excluding_hydrogen_ = items.at(7);
    number_of_hydrogen_parameters_ = items.at(8);
    number_of_parameters_ = items.at(9);
    number_of_excluded_atoms_ = items.at(10);
    number_of_residues_ = items.at(11);
    total_number_of_bonds_ = items.at(12);
    total_number_of_angles_ = items.at(13);
    total_number_of_dihedrals_ = items.at(14);
    number_of_bond_types_ = items.at(15);
    number_of_angle_types_ = items.at(16);
    number_of_dihedral_types_ = items.at(17);
    number_of_atom_types_in_parameter_file_ = items.at(18);
    number_of_distinct_hydrogen_bonds_ = items.at(19);
    perturbation_option_ = items.at(20);
    number_of_bonds_perturbed_ = items.at(21);
    number_of_angles_perturbed_ = items.at(22);
    number_of_dihedrals_perturbed_ = items.at(23);
    number_of_bonds_group_perturbed_ = items.at(24);
    number_of_angles_group_perturbed_ = items.at(25);
    number_of_dihedrals_group_perturbed_ = items.at(26);
    standard_periodic_box_option_ = items.at(27);
    number_of_atoms_in_largest_residue_ = items.at(28);
    cap_option_ = items.at(29);
    number_of_extra_points_ = items.at(30);
    if(items.size() >= 31)
        number_of_beads_ = items.at(30);
}

template<typename T>
vector<T> TopologyFile::ParsePartition(stringstream &stream)
{
    string line;
    vector<T> section = vector<T>();
    getline(stream, line);
    string format;
    while(!line.empty())
    {
        if(line[0] == '%')
        {
            if(line.find("%FORMAT") != string::npos)
            {
                format = line.substr(line.find_first_of('(', 0) + 1, line.find_last_of(')', line.length()) - line.find_first_of('(', 0) - 1);
            }
        }
        else
        {
            vector<T> line_items = PartitionLine<T>(line, format);
            for(typename vector<T>::iterator it = line_items.begin(); it != line_items.end(); it++)
            {
                section.push_back(*it);
            }
        }
        getline(stream, line);
    }
    return section;
}


template<typename T>
vector<T> TopologyFile::PartitionLine(string line, string format)
{
    vector<T> items = vector<T>();
    if(format.compare("10I8") == 0)
    {
        int number_of_items = 10;
        int item_length = 8;
        for(int i = 0; i < number_of_items && item_length * (i+1) <= line.length(); i++)
        {
            string token = line.substr(i*item_length, item_length);
            token = Trim(token);
            items.push_back(ConvertString<T>(token));
        }
        return items;
    }
    if(format.compare("20a4") == 0)
    {
        int number_of_items = 20;
        int item_length = 4;
        for(int i = 0; i < number_of_items && item_length * (i+1) <= line.length(); i++)
        {
            string token = line.substr(i*item_length, item_length);
            token = Trim(token);
            items.push_back(ConvertString<T>(token));
        }
        return items;
    }
    if(format.compare("5E16.8") == 0)
    {
        int number_of_items = 5;
        int item_length = 16;
        for(int i = 0; i < number_of_items && item_length * (i+1) <= line.length(); i++)
        {
            string token = line.substr(i*item_length, item_length);
            token = Trim(token);
            double base = ConvertString<double>(Split(token, "E")[0]);
            int power = ConvertString<int>(Split(token, "E")[1]);
            stringstream temp, format;
            double result = base * pow10(power);
            format << setw(16) << fixed << setprecision(8);
            temp << format.str() << result;
            T item = ConvertString<T>(temp.str());
            items.push_back(item);
        }
        return items;
    }
    if(format.compare("1a80") == 0)
    {
        int number_of_items = 1;
        int item_length = 80;
        for(int i = 0; i < number_of_items && item_length * (i+1) <= line.length(); i++)
        {
            string token = line.substr(i*item_length, item_length);
            token = Trim(token);
            items.push_back(ConvertString<T>(token));
        }
        return items;
    }
    if(format.compare("1I8") == 0)
    {
        int number_of_items = 1;
        int item_length = 8;
        for(int i = 0; i < number_of_items && item_length * (i+1) <= line.length(); i++)
        {
            string token = line.substr(i*item_length, item_length);
            token = Trim(token);
            items.push_back(ConvertString<T>(token));
        }
        return items;
    }
}

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void TopologyFile::Print(ostream &out)
{
    out << "============================ " << title_ << " ===========================" << endl;
    assembly_->Print(out);
    out << "------------------------------ Bond Types ------------------------------" << endl;
    for(TopologyBondTypeMap::iterator it = bond_types_.begin(); it != bond_types_.end(); it++)
    {
        TopologyBondType* bond_type = (*it).second;
        bond_type->Print(out);
        out << endl;
    }
    out << "------------------------------ Angle Types ------------------------------" << endl;
    for(TopologyAngleTypeMap::iterator it = angle_types_.begin(); it != angle_types_.end(); it++)
    {
        TopologyAngleType* angle_type = (*it).second;
        angle_type->Print(out);
        out << endl;
    }
    out << "------------------------------ Dihedral Types ------------------------------" << endl;
    for(TopologyDihedralTypeMap::iterator it = dihedral_types_.begin(); it != dihedral_types_.end(); it++)
    {
        TopologyDihedralType* dihedral_type = (*it).second;
        dihedral_type->Print(out);
        out << endl;
    }
    out << "------------------------------ Bonds ------------------------------" << endl;
    for(TopologyBondMap::iterator it = bonds_.begin(); it != bonds_.end(); it++)
    {
        TopologyBond* bond = (*it).second;
        bond->Print(out);
    }
    out << "------------------------------ Angles ------------------------------" << endl;
    for(TopologyAngleMap::iterator it = angles_.begin(); it != angles_.end(); it++)
    {
        TopologyAngle* angle = (*it).second;
        angle->Print(out);
    }
    out << "------------------------------ Dihedrals ------------------------------" << endl;
    for(TopologyDihedralMap::iterator it = dihedrals_.begin(); it != dihedrals_.end(); it++)
    {
        TopologyDihedral* dihedral = (*it).second;
        dihedral->Print(out);
    }

}
