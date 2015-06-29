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
TopologyFile::TopologyFile()
{
    title_ = "";
    path_ = "";
    number_of_atoms_ = iNotSet;
    number_of_types_ = iNotSet;
    number_of_bonds_including_hydrogen_ = iNotSet;
    number_of_angles_including_hydrogen_ = iNotSet;
    number_of_dihedrals_including_hydrogen_ = iNotSet;
    number_of_bonds_excluding_hydrogen_ = iNotSet;
    number_of_angles_excluding_hydrogen_ = iNotSet;
    number_of_dihedrals_excluding_hydrogen_ = iNotSet;
    number_of_hydrogen_parameters_ = iNotSet;
    number_of_parameters_ = iNotSet;
    number_of_excluded_atoms_ = iNotSet;
    number_of_residues_ = iNotSet;
    total_number_of_bonds_ = iNotSet;
    total_number_of_angles_ = iNotSet;
    total_number_of_dihedrals_ = iNotSet;
    number_of_bond_types_ = iNotSet;
    number_of_angle_types_ = iNotSet;
    number_of_dihedral_types_ = iNotSet;
    number_of_atom_types_in_parameter_file_ = iNotSet;
    number_of_distinct_hydrogen_bonds_ = iNotSet;
    perturbation_option_ = iNotSet;
    number_of_bonds_perturbed_ = iNotSet;
    number_of_angles_perturbed_ = iNotSet;
    number_of_dihedrals_perturbed_ = iNotSet;
    number_of_bonds_group_perturbed_ = iNotSet;
    number_of_angles_group_perturbed_ = iNotSet;
    number_of_dihedrals_group_perturbed_ = iNotSet;
    standard_periodic_box_option_ = iNotSet;
    number_of_atoms_in_largest_residue_ = iNotSet;
    cap_option_ = iNotSet;
    number_of_extra_points_ = iNotSet;
    number_of_beads_ = iNotSet;
    pairs_ = TopologyAtomPairMap();
    bond_types_ = TopologyBondTypeMap();
    angle_types_ = TopologyAngleTypeMap();
    dihedral_types_ = TopologyDihedralTypeMap();
    assembly_ = new TopologyAssembly();
    radius_set_ = RadiusSet();
    bonds_ = TopologyBondMap();
    angles_ = TopologyAngleMap();
    dihedrals_ = TopologyDihedralMap();
}

TopologyFile::TopologyFile(const string &top_file)
{
    path_ = top_file;
    ifstream in_file;        
    if(std::ifstream(top_file.c_str()))
        in_file.open(top_file.c_str());
    else
    {
        throw TopologyFileProcessingException(__LINE__, "Topology file not found");
    }
    Read(in_file);
    in_file.close();            /// Close the parameter files
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
string TopologyFile::GetPath()
{
    return path_;
}
string TopologyFile::GetTitle()
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
TopologyFile::TopologyAtomTypeIndexMap TopologyFile::GetAtomTypeIndexMap()
{
    TopologyAtomTypeIndexMap atom_type_index_map;
    int index = 0;
    int count = 1;
    for(int i = 0; i < number_of_residues_; i++)
    {
        TopologyResidue* residue = this->assembly_->GetResidueByIndex(i+1);
        for(unsigned int j = 0; j < residue->GetAtoms().size(); j++)
        {
            TopologyAtom* atom = residue->GetAtomByIndex(index+1);
            index++;
            if(atom_type_index_map[atom->GetType()] == 0)
            {
                atom_type_index_map[atom->GetType()] = count;
                count++;
            }
        }
    }
    return atom_type_index_map;
}
TopologyBondType* TopologyFile::GetBondTypeByIndex(int index)
{
    for(TopologyBondTypeMap::iterator it = bond_types_.begin(); it != bond_types_.end(); it++)
    {
        TopologyBondType* bond_type = (*it).second;
        if(bond_type->GetIndex() == index)
            return bond_type;
    }
    return NULL;
}
TopologyAngleType* TopologyFile::GetAngleTypeByIndex(int index)
{
    for(TopologyAngleTypeMap::iterator it = angle_types_.begin(); it != angle_types_.end(); it++)
    {
        TopologyAngleType* angle_type = (*it).second;
        if(angle_type->GetIndex() == index)
            return angle_type;
    }
    return NULL;
}
TopologyDihedralType* TopologyFile::GetDihedralTypeByIndex(int index)
{
    for(TopologyDihedralTypeMap::iterator it = dihedral_types_.begin(); it != dihedral_types_.end(); it++)
    {
        TopologyDihedralType* dihedral_type = (*it).second;
        if(dihedral_type->GetIndex() == index)
        {
            return dihedral_type;
        }
    }
    return NULL;
}
TopologyAtomPair* TopologyFile::GetAtomPairByIndex(int index)
{
    for(TopologyAtomPairMap::iterator it = pairs_.begin(); it != pairs_.end(); it++)
    {
        TopologyAtomPair* atom_pair = (*it).second;
        if(atom_pair->GetIndex() == index)
            return atom_pair;
    }
    return NULL;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void TopologyFile::SetPath(string path)
{
    path_ = path;
}
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
        string atom_bond = (*it).first;
        TopologyBond* bond = (*it).second;
        bonds_[atom_bond] = bond;
    }
}
void TopologyFile::AddBond(TopologyBond *bond)
{
    vector<string> bonds = bond->GetBonds();
    vector<string> residue_names = bond->GetResidueNames();
    TopologyBondType* bond_type = bond->GetBondType();
    stringstream ss;
    ss << residue_names.at(0) << ":" << bonds.at(0) << "-" << residue_names.at(1) << ":" << bonds.at(1) << "_" << bond_type;
    bonds_[ss.str()] = bond;
}
void TopologyFile::SetBondTypes(TopologyBondTypeMap bond_types)
{
    bond_types_.clear();
    for(TopologyBondTypeMap::iterator it = bond_types.begin(); it != bond_types.end(); it++)
    {
        int index = (*it).first;
        TopologyBondType* bond_type = (*it).second;
        bond_types_[index] = bond_type;
    }
}
void TopologyFile::AddBondType(TopologyBondType *bond_type)
{
    int index = bond_type->GetIndex();
    bond_types_[index] = bond_type;
}
void TopologyFile::SetAngles(TopologyAngleMap angles)
{
    angles_.clear();
    for(TopologyAngleMap::iterator it = angles.begin(); it != angles.end(); it++)
    {
        string atom_angle = (*it).first;
        TopologyAngle* angle = (*it).second;
        angles_[atom_angle] = angle;
    }
}
void TopologyFile::AddAngle(TopologyAngle *angle)
{
    vector<string> angles = angle->GetAngles();
    vector<string> residue_names = angle->GetResidueNames();
    TopologyAngleType* angle_type = angle->GetAngleType();
    stringstream ss;
    ss << residue_names.at(0) << ":" << angles.at(0) << "-" << residue_names.at(1) << ":" << angles.at(1) << "-" << residue_names.at(2) << ":" << angles.at(2) << "_" << angle_type;
    angles_[ss.str()] = angle;
}
void TopologyFile::SetAngleTypes(TopologyAngleTypeMap angle_types)
{
    angle_types_.clear();
    for(TopologyAngleTypeMap::iterator it = angle_types.begin(); it != angle_types.end(); it++)
    {
        int index = (*it).first;
        TopologyAngleType* angle_type = (*it).second;
        angle_types_[index] = angle_type;
    }
}
void TopologyFile::AddAngleType(TopologyAngleType *angle_type)
{
    int index = angle_type->GetIndex();
    angle_types_[index] = angle_type;
}
void TopologyFile::SetDihedrals(TopologyDihedralMap dihedrals)
{
    dihedrals_.clear();
    for(TopologyDihedralMap::iterator it = dihedrals.begin(); it != dihedrals.end(); it++)
    {
        string atom_dihedral = (*it).first;
        TopologyDihedral* dihedral = (*it).second;
        dihedrals_[atom_dihedral] = dihedral;
    }
}
void TopologyFile::AddDihedral(TopologyDihedral *dihedral)
{
    vector<string> dihedrals = dihedral->GetDihedrals();
    vector<string> residue_names = dihedral->GetResidueNames();
    TopologyDihedralType* dihedral_type = dihedral->GetDihedralType();
    stringstream ss;
    ss << residue_names.at(0) << ":" << dihedrals.at(0) << "-" << residue_names.at(1) << ":" << dihedrals.at(1) << "-"
       << residue_names.at(2) << ":" << dihedrals.at(2) << "-" << residue_names.at(3) << ":" << dihedrals.at(3) << "_" << dihedral_type;
    dihedrals_[ss.str()] = dihedral;
}
void TopologyFile::SetDihedralTypes(TopologyDihedralTypeMap dihedral_types)
{
    dihedral_types_.clear();
    for(TopologyDihedralTypeMap::iterator it = dihedral_types.begin(); it != dihedral_types.end(); it++)
    {
        int index = (*it).first;
        TopologyDihedralType* dihedral_type = (*it).second;
        dihedral_types_[index] = dihedral_type;
    }
}
void TopologyFile::AddDihedralType(TopologyDihedralType *dihedral_type)
{
    int index = dihedral_type->GetIndex();
    dihedral_types_[index] = dihedral_type;
}
void TopologyFile::SetAtomPairs(TopologyAtomPairMap pairs)
{
    pairs_.clear();
    for(TopologyAtomPairMap::iterator it = pairs.begin(); it != pairs.end(); it++)
    {
        string atom_pair_types = (*it).first;
        TopologyAtomPair* pair = (*it).second;
        pairs_[atom_pair_types] = pair;
    }
}
void TopologyFile::AddAtomPair(TopologyAtomPair *pair)
{
    string atom_pair_types = pair->GetPairType();
    pairs_[atom_pair_types] = pair;
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
            else if(in_line.find("%FLAG IPOL") != string::npos)
            {
                ipols = ParsePartition<int>(section);
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
            if(pairs_.find(atom_type_1 + "-" + atom_type_2) != pairs_.end()|| pairs_.find(atom_type_2 + "-" + atom_type_1) != pairs_.end())
            {

            }
            else
            {
                pairs_[atom_type_1 + "-" + atom_type_2] = new TopologyAtomPair(atom_type_1 + "-" + atom_type_2, coefficient_a, coefficient_b, index);
//                cout << atom_type_1 + "-" + atom_type_2 << ":" << index << endl;
            }
        }

    }
    // Bonds, Angles, Dihedrals
    //Bonds in topology file
    for(int i = 0; i < number_of_bonds_including_hydrogen_; i++)
    {
        vector<string> bonds = vector<string>();
        vector<string> bond_atoms_residue_names = vector<string>();

        int first_atom_index = (bonds_inc_hydrogens[i*3])/3;
        int second_atom_index = (bonds_inc_hydrogens[i*3+1])/3;

        bonds.push_back(atom_names[first_atom_index]);
        bonds.push_back(atom_names[second_atom_index]);

        for(int j = 0; j < number_of_residues_; j++)
        {
            int start_index = residue_pointers.at(j) - 1;
            int end_index;
            if(j < number_of_residues_ - 1)
            {
                end_index = residue_pointers.at(j+1) - 1;
            }
            else
            {
                end_index = number_of_atoms_;
            }
            if(first_atom_index >= start_index && first_atom_index < end_index)
            {
                stringstream r;
                r << residue_labels.at(j) << "(" << (j+1) << ")";
                bond_atoms_residue_names.push_back(r.str());
                break;
            }
        }
        for(int j = 0; j < number_of_residues_; j++)
        {
            int start_index = residue_pointers.at(j) - 1;
            int end_index;
            if(j < number_of_residues_ - 1)
            {
                end_index = residue_pointers.at(j+1) - 1;
            }
            else
            {
                end_index = number_of_atoms_;
            }
            if(second_atom_index >= start_index && second_atom_index < end_index)
            {
                stringstream r;
                r << residue_labels.at(j) << "(" << (j+1) << ")";
                bond_atoms_residue_names.push_back(r.str());
                break;
            }
        }

        TopologyBond* topology_bond = new TopologyBond(bonds, bond_atoms_residue_names);

        topology_bond->SetBondType(bond_types_[bonds_inc_hydrogens[i*3+2] - 1]);
        topology_bond->SetIncludingHydrogen(true);

        stringstream key;
        key << bond_atoms_residue_names.at(0) << ":" << bonds.at(0) << "-" << bond_atoms_residue_names.at(1) << ":" << bonds.at(1) << "_" << topology_bond->GetBondType();
        bonds_[key.str()] = topology_bond;
    }
    for(int i = 0; i < number_of_bonds_excluding_hydrogen_; i++)
    {
        vector<string> bonds = vector<string>();
        vector<string> bond_atoms_residue_names = vector<string>();

        int first_atom_index = (bonds_without_hydrogens[i*3])/3;
        int second_atom_index = (bonds_without_hydrogens[i*3+1])/3;

        bonds.push_back(atom_names[first_atom_index]);
        bonds.push_back(atom_names[second_atom_index]);

        for(int j = 0; j < number_of_residues_; j++)
        {
            int start_index = residue_pointers.at(j) - 1;
            int end_index;
            if(j <number_of_residues_ - 1)
            {
                end_index = residue_pointers.at(j+1) - 1;
            }
            else
            {
                end_index = number_of_atoms_;
            }
            if(first_atom_index >= start_index && first_atom_index < end_index)
            {
                stringstream r;
                r << residue_labels.at(j) << "(" << (j+1) << ")";
                bond_atoms_residue_names.push_back(r.str());
                break;
            }
        }
        for(int j = 0; j < number_of_residues_; j++)
        {
            int end_index;
            int start_index = residue_pointers.at(j) - 1;
            if(j < number_of_residues_ - 1)
            {
                end_index = residue_pointers.at(j+1) - 1;
            }
            else
            {
                end_index = number_of_atoms_;
            }
            if(second_atom_index >= start_index && second_atom_index < end_index)
            {
                stringstream r;
                r << residue_labels.at(j) << "(" << (j+1) << ")";
                bond_atoms_residue_names.push_back(r.str());
                break;
            }
        }

        TopologyBond* topology_bond = new TopologyBond(bonds, bond_atoms_residue_names);

        topology_bond->SetBondType(bond_types_[bonds_without_hydrogens[i*3+2] - 1]);
        topology_bond->SetIncludingHydrogen(false);

        stringstream key;
        key << bond_atoms_residue_names.at(0) << ":" << bonds.at(0) << "-" << bond_atoms_residue_names.at(1) << ":" << bonds.at(1) << "_" << topology_bond->GetBondType();
        bonds_[key.str()] = topology_bond;
    }
    // Angles in topology file
    // Angles including hydrogen
    for(int i = 0; i < number_of_angles_including_hydrogen_; i++)
    {
        vector<string> angle_atoms = vector<string>();
        vector<string> angle_atoms_residue_names = vector<string>();

        int atom_index_1 = angles_inc_hydrogens.at(i*4) / 3;
        int atom_index_2 = angles_inc_hydrogens.at(i*4+1) / 3;
        int atom_index_3 = angles_inc_hydrogens.at(i*4+2) / 3;

        angle_atoms.push_back(atom_names.at(atom_index_1));
        angle_atoms.push_back(atom_names.at(atom_index_2));
        angle_atoms.push_back(atom_names.at(atom_index_3));

        for(int j = 0; j < number_of_residues_; j++)
        {
            int start_index = residue_pointers.at(j) - 1;
            int end_index;
            if(j < number_of_residues_ - 1)
            {
                end_index = residue_pointers.at(j+1) - 1;
            }
            else
            {
                end_index = number_of_atoms_;
            }
            if(atom_index_1 >= start_index && atom_index_1 < end_index)
            {
                stringstream r;
                r << residue_labels.at(j) << "(" << (j+1) << ")";
                angle_atoms_residue_names.push_back(r.str());
                break;
            }
        }
        for(int j = 0; j < number_of_residues_; j++)
        {
            int start_index = residue_pointers.at(j) - 1;
            int end_index;
            if(j < number_of_residues_ - 1)
            {
                end_index = residue_pointers.at(j+1) - 1;
            }
            else
            {
                end_index = number_of_atoms_;
            }
            if(atom_index_2 >= start_index && atom_index_2 < end_index)
            {
                stringstream r;
                r << residue_labels.at(j) << "(" << (j+1) << ")";
                angle_atoms_residue_names.push_back(r.str());
                break;
            }
        }
        for(int j = 0; j < number_of_residues_; j++)
        {
            int start_index = residue_pointers.at(j) - 1;
            int end_index;
            if(j < number_of_residues_ - 1)
            {
                end_index = residue_pointers.at(j+1) - 1;
            }
            else
            {
                end_index = number_of_atoms_;
            }
            if(atom_index_3 >= start_index && atom_index_3 < end_index)
            {
                stringstream r;
                r << residue_labels.at(j) << "(" << (j+1) << ")";
                angle_atoms_residue_names.push_back(r.str());
                break;
            }
        }

        TopologyAngle* angle = new TopologyAngle(angle_atoms, angle_atoms_residue_names);

        angle->SetIncludingHydrogen(true);
        angle->SetAnlgeType(angle_types_[angles_inc_hydrogens.at(i*4+3) - 1]);

        stringstream key;
        key << angle_atoms_residue_names.at(0) << ":" << angle_atoms.at(0) << "-" << angle_atoms_residue_names.at(1) << ":" << angle_atoms.at(1) << "-"
            << angle_atoms_residue_names.at(2) << ":" << angle_atoms.at(2) << "_" << angle->GetAngleType();
        angles_[key.str()] = angle;
    }
    // Angles excluding hydrogen
    for(int i = 0; i < number_of_angles_excluding_hydrogen_; i++)
    {
        vector<string> angle_atoms = vector<string>();
        vector<string> angle_atoms_residue_names = vector<string>();

        int atom_index_1 = angles_without_hydrogens.at(i*4) / 3;
        int atom_index_2 = angles_without_hydrogens.at(i*4+1) / 3;
        int atom_index_3 = angles_without_hydrogens.at(i*4+2) / 3;

        angle_atoms.push_back(atom_names.at(atom_index_1));
        angle_atoms.push_back(atom_names.at(atom_index_2));
        angle_atoms.push_back(atom_names.at(atom_index_3));

        for(int j = 0; j < number_of_residues_; j++)
        {
            int start_index = residue_pointers.at(j) - 1;
            int end_index;
            if(j < number_of_residues_ - 1)
            {
                end_index = residue_pointers.at(j+1) - 1;
            }
            else
            {
                end_index = number_of_atoms_;
            }
            if(atom_index_1 >= start_index && atom_index_1 < end_index)
            {
                stringstream r;
                r << residue_labels.at(j) << "(" << (j+1) << ")";
                angle_atoms_residue_names.push_back(r.str());
                break;
            }
        }
        for(int j = 0; j < number_of_residues_; j++)
        {
            int start_index = residue_pointers.at(j) - 1;
            int end_index;
            if(j < number_of_residues_ - 1)
            {
                end_index = residue_pointers.at(j+1) - 1;
            }
            else
            {
                end_index = number_of_atoms_;
            }
            if(atom_index_2 >= start_index && atom_index_2 < end_index)
            {
                stringstream r;
                r << residue_labels.at(j) << "(" << (j+1) << ")";
                angle_atoms_residue_names.push_back(r.str());
                break;
            }
        }
        for(int j = 0; j < number_of_residues_; j++)
        {
            int end_index;
            int start_index = residue_pointers.at(j) - 1;
            if(j < number_of_residues_ - 1)
            {
                end_index = residue_pointers.at(j+1) - 1;
            }
            else
            {
                end_index = number_of_atoms_;
            }
            if(atom_index_3 >= start_index && atom_index_3 < end_index)
            {
                stringstream r;
                r << residue_labels.at(j) << "(" << (j+1) << ")";
                angle_atoms_residue_names.push_back(r.str());
                break;
            }
        }

        TopologyAngle* angle = new TopologyAngle(angle_atoms, angle_atoms_residue_names);

        angle->SetIncludingHydrogen(false);
        angle->SetAnlgeType(angle_types_[angles_without_hydrogens.at(i*4+3) - 1]);

        stringstream key;
        key << angle_atoms_residue_names.at(0) << ":" << angle_atoms.at(0) << "-" << angle_atoms_residue_names.at(1) << ":" << angle_atoms.at(1) << "-"
            << angle_atoms_residue_names.at(2) << ":" << angle_atoms.at(2) << "_" << angle->GetAngleType();
        angles_[key.str()] = angle;
    }
    // Dihedrals in topology file
    // Dihedrals including hydrogen
    for(int i = 0; i < number_of_dihedrals_including_hydrogen_; i++)
    {
        vector<string> dihedral_atoms = vector<string>();
        vector<string> dihedral_atoms_residue_names = vector<string>();

        int atom_index_1 = abs(dihedrals_inc_hydrogens.at(i*5)) / 3;
        int atom_index_2 = abs(dihedrals_inc_hydrogens.at(i*5+1)) / 3;
        int atom_index_3 = abs(dihedrals_inc_hydrogens.at(i*5+2)) / 3;
        int atom_index_4 = abs(dihedrals_inc_hydrogens.at(i*5+3)) / 3;

        dihedral_atoms.push_back(atom_names.at(atom_index_1));
        dihedral_atoms.push_back(atom_names.at(atom_index_2));
        dihedral_atoms.push_back(atom_names.at(atom_index_3));
        dihedral_atoms.push_back(atom_names.at(atom_index_4));

        for(int j = 0; j < number_of_residues_; j++)
        {
            int start_index = residue_pointers.at(j) - 1;
            int end_index;
            if(j < number_of_residues_ - 1)
            {
                end_index = residue_pointers.at(j+1) - 1;
            }
            else
            {
                end_index = number_of_atoms_;
            }
            if(atom_index_1 >= start_index && atom_index_1 < end_index)
            {
                stringstream r;
                r << residue_labels.at(j) << "(" << (j+1) << ")";
                dihedral_atoms_residue_names.push_back(r.str());
                break;
            }
        }
        for(int j = 0; j < number_of_residues_; j++)
        {
            int start_index = residue_pointers.at(j) - 1;
            int end_index;
            if(j < number_of_residues_ - 1)
            {
                end_index = residue_pointers.at(j+1) - 1;
            }
            else
            {
                end_index = number_of_atoms_;
            }
            if(atom_index_2 >= start_index && atom_index_2 < end_index)
            {
                stringstream r;
                r << residue_labels.at(j) << "(" << (j+1) << ")";
                dihedral_atoms_residue_names.push_back(r.str());
                break;
            }
        }
        for(int j = 0; j < number_of_residues_; j++)
        {
            int start_index = residue_pointers.at(j) - 1;
            int end_index;
            if(j < number_of_residues_ - 1)
            {
                end_index = residue_pointers.at(j+1) - 1;
            }
            else
            {
                end_index = number_of_atoms_;
            }
            if(atom_index_3 >= start_index && atom_index_3 < end_index)
            {
                stringstream r;
                r << residue_labels.at(j) << "(" << (j+1) << ")";
                dihedral_atoms_residue_names.push_back(r.str());
                break;
            }
        }
        for(int j = 0; j < number_of_residues_; j++)
        {
            int start_index = residue_pointers.at(j) - 1;
            int end_index;
            if(j < number_of_residues_ - 1)
            {
                end_index = residue_pointers.at(j+1) - 1;
            }
            else
            {
                end_index = number_of_atoms_;
            }
            if(atom_index_4 >= start_index && atom_index_4 < end_index)
            {
                stringstream r;
                r << residue_labels.at(j) << "(" << (j+1) << ")";
                dihedral_atoms_residue_names.push_back(r.str());
                break;
            }
        }

        TopologyDihedral* dihedral = new TopologyDihedral(dihedral_atoms, dihedral_atoms_residue_names);

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
        stringstream key;
        key << dihedral_atoms_residue_names.at(0) << ":" << dihedral_atoms.at(0) << "-" << dihedral_atoms_residue_names.at(1) << ":" << dihedral_atoms.at(1) << "-"
            << dihedral_atoms_residue_names.at(2) << ":" << dihedral_atoms.at(2) << "-" << dihedral_atoms_residue_names.at(3) << ":" << dihedral_atoms.at(3) << "_"
            << dihedral->GetDihedralType() << "_"
            << (dihedral->GetIsImproper() ? "IY" : "IN") << "_" << (dihedral->GetIgnoredGroupInteraction() ? "GY" : "GN");
        dihedrals_[key.str()] = dihedral;

    }
    // Dihedrals excluding hydrogen
    for(int i = 0; i < number_of_dihedrals_excluding_hydrogen_; i++)
    {
        vector<string> dihedral_atoms = vector<string>();
        vector<string> dihedral_atoms_residue_names = vector<string>();

        int atom_index_1 = abs(dihedrals_without_hydrogens.at(i*5)) / 3;
        int atom_index_2 = abs(dihedrals_without_hydrogens.at(i*5+1)) / 3;
        int atom_index_3 = abs(dihedrals_without_hydrogens.at(i*5+2)) / 3;
        int atom_index_4 = abs(dihedrals_without_hydrogens.at(i*5+3)) / 3;

        dihedral_atoms.push_back(atom_names.at(atom_index_1));
        dihedral_atoms.push_back(atom_names.at(atom_index_2));
        dihedral_atoms.push_back(atom_names.at(atom_index_3));
        dihedral_atoms.push_back(atom_names.at(atom_index_4));

        for(int j = 0; j < number_of_residues_; j++)
        {
            int start_index = residue_pointers.at(j) - 1;
            int end_index;
            if(j < number_of_residues_ - 1)
            {
                end_index = residue_pointers.at(j+1) - 1;
            }
            else
            {
                end_index = number_of_atoms_;
            }
            if(atom_index_1 >= start_index && atom_index_1 < end_index)
            {
                stringstream r;
                r << residue_labels.at(j) << "(" << (j+1) << ")";
                dihedral_atoms_residue_names.push_back(r.str());
                break;
            }
        }
        for(int j = 0; j < number_of_residues_; j++)
        {
            int start_index = residue_pointers.at(j) - 1;
            int end_index;
            if(j < number_of_residues_ - 1)
            {
                end_index = residue_pointers.at(j+1) - 1;
            }
            else
            {
                end_index = number_of_atoms_;
            }
            if(atom_index_2 >= start_index && atom_index_2 < end_index)
            {
                stringstream r;
                r << residue_labels.at(j) << "(" << (j+1) << ")";
                dihedral_atoms_residue_names.push_back(r.str());
                break;
            }
        }
        for(int j = 0; j < number_of_residues_; j++)
        {
            int start_index = residue_pointers.at(j) - 1;
            int end_index;
            if(j < number_of_residues_ - 1)
            {
                end_index = residue_pointers.at(j+1) - 1;
            }
            else
            {
                end_index = number_of_atoms_;
            }
            if(atom_index_3 >= start_index && atom_index_3 < end_index)
            {
                stringstream r;
                r << residue_labels.at(j) << "(" << (j+1) << ")";
                dihedral_atoms_residue_names.push_back(r.str());
                break;
            }
        }
        for(int j = 0; j < number_of_residues_; j++)
        {
            int start_index = residue_pointers.at(j) - 1;
            int end_index;
            if(j < number_of_residues_ - 1)
            {
                end_index = residue_pointers.at(j+1) - 1;
            }
            else
            {
                end_index = number_of_atoms_;
            }
            if(atom_index_4 >= start_index && atom_index_4 < end_index)
            {
                stringstream r;
                r << residue_labels.at(j) << "(" << (j+1) << ")";
                dihedral_atoms_residue_names.push_back(r.str());
                break;
            }
        }

        TopologyDihedral* dihedral = new TopologyDihedral(dihedral_atoms, dihedral_atoms_residue_names);

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

        stringstream key;
        key << dihedral_atoms_residue_names.at(0) << ":" << dihedral_atoms.at(0) << "-" << dihedral_atoms_residue_names.at(1) << ":" << dihedral_atoms.at(1) << "-"
            << dihedral_atoms_residue_names.at(2) << ":" << dihedral_atoms.at(2) << "-" << dihedral_atoms_residue_names.at(3) << ":" << dihedral_atoms.at(3) << "_"
            << dihedral->GetDihedralType() << "_"
            << (dihedral->GetIsImproper() ? "IY" : "IN") << "_" << (dihedral->GetIgnoredGroupInteraction() ? "GY" : "GN");
        dihedrals_[key.str()] = dihedral;
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
            ending_atom_index = number_of_atoms_ + 1;
        }
        for(int i = starting_atom_index - 1; i < ending_atom_index - 1; i++)
        {
            TopologyAtom::ExcludedAtomNames excluded_atoms = TopologyAtom::ExcludedAtomNames();
            if(i > 0)
                start_index += number_excluded_atoms.at(i-1);
            for(int j = start_index; j < start_index + number_excluded_atoms.at(i); j++)
            {
                string excluded_atom_residue_name;
                int index = (excluded_atoms_lists.at(j) - 1 == -1) ? 0 : excluded_atoms_lists.at(j) - 1;
                for(int k = 0; k < number_of_residues_; k++)
                {
                    int start_atom_index = residue_pointers.at(k) - 1;
                    int end_atom_index;
                    if(k < number_of_residues_ - 1)
                    {
                        end_atom_index = residue_pointers.at(k+1) - 1;
                    }
                    else
                    {
                        end_atom_index = number_of_atoms_;
                    }
                    if(index >= start_atom_index && index < end_atom_index)
                    {
                        excluded_atom_residue_name = residue_labels.at(k);
                        break;
                    }
                }
                excluded_atoms.push_back(excluded_atom_residue_name + ":" + atom_names.at(index));
            }
            atoms[atom_names.at(i)] = new TopologyAtom(i + 1, atom_names.at(i), amber_atom_types.at(i), charges.at(i), atomic_numbers.at(i), masses.at(i), excluded_atoms,
                                                       number_excluded_atoms.at(i), radiis.at(i), screens.at(i), tree_chain_classifications.at(i), residue_name);
        }
        stringstream residue_key;
        residue_key << residue_name << "_" << residue_index;
        residues[residue_key.str()] = new TopologyResidue(residue_name, atoms, residue_index, starting_atom_index);
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
        for(int i = 0; i < number_of_items && item_length * (i+1) <= (int)line.length(); i++)
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
        for(int i = 0; i < number_of_items && item_length * (i+1) <= (int)line.length(); i++)
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
        for(int i = 0; i < number_of_items && item_length * (i+1) <= (int)line.length(); i++)
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
        for(int i = 0; i < number_of_items && item_length * (i+1) <= (int)line.length(); i++)
        {
            string token = line.substr(i*item_length, item_length);
            std::replace_if(token.begin(), token.end(), ::isspace, '#');
            token = Trim(token);
            items.push_back(ConvertString<T>(token));
        }

        return items;
    }
    if(format.compare("1I8") == 0)
    {
        int number_of_items = 1;
        int item_length = 8;
        for(int i = 0; i < number_of_items && item_length * (i+1) <= (int)line.length(); i++)
        {
            string token = line.substr(i*item_length, item_length);
            token = Trim(token);
            items.push_back(ConvertString<T>(token));
        }
        return items;
    }
    return items;
}

void TopologyFile::Write(const string &top_file)
{
    std::ofstream out_file;
    try
    {
        out_file.open(top_file.c_str());
    }
    catch(...)
    {
        throw TopologyFileProcessingException(__LINE__,"File could not be created");
    }
    try
    {
        this->ResolveSections(out_file);
    }
    catch(...)
    {
        out_file.close();            /// Close the parm7 files
    }
}

void TopologyFile::ResolveSections(ofstream &out_stream)
{
    this->ResolveTitleSection(out_stream);
    this->ResolvePointersSection(out_stream);
    this->ResolveAtomNameSection(out_stream);
    this->ResolveChargeSection(out_stream);
    this->ResolveAtomicNumberSection(out_stream);
    this->ResolveMassSection(out_stream);
    this->ResolveAtomTypeIndexSection(out_stream);
    this->ResolveNumberExcludedAtomsSection(out_stream);    
    this->ResolveNonbondedParmIndexSection(out_stream);
    this->ResolveResidueLabelSection(out_stream);
    this->ResolveResiduePointersSection(out_stream);
    this->ResolveBondForceConstantSection(out_stream);
    this->ResolveBondEquilValueSection(out_stream);
    this->ResolveAngleForceConstantSection(out_stream);
    this->ResolveAngleEquilValueSection(out_stream);
    this->ResolveDihedralForceConstantSection(out_stream);
    this->ResolveDihedralPeriodicitySection(out_stream);
    this->ResolveDihedralPhaseSection(out_stream);
    this->ResolveSceeScaleFactorSection(out_stream);
    this->ResolveScnbScaleFactorSection(out_stream);
    this->ResolveSoltySection(out_stream);
    this->ResolveLennardJonesACoefSection(out_stream);
    this->ResolveLennardJonesBCoefSection(out_stream);
    this->ResolveBondsIncHydrogenSection(out_stream);
    this->ResolveBondsWithoutHydrogenSection(out_stream);
    this->ResolveAnglesIncHydrogenSection(out_stream);
    this->ResolveAnglesWithoutHydrogenSection(out_stream);
    this->ResolveDihedralsIncHydrogenSection(out_stream);
    this->ResolveDihedralsWithoutHydrogenSection(out_stream);
    this->ResolveExcludedAtomsListSection(out_stream);
    this->ResolveHydrogenBondACoefSection(out_stream);
    this->ResolveHydrogenBondBCoefSection(out_stream);
    this->ResolveHBCutSection(out_stream);
    this->ResolveAmberAtomTypeSection(out_stream);
    this->ResolveTreeChainClassificationSection(out_stream);
    this->ResolveJoinArraySection(out_stream);
    this->ResolveIRotatSection(out_stream);
    this->ResolveRadiusSetSection(out_stream);
    this->ResolveRadiiSection(out_stream);
    this->ResolveScreenSection(out_stream);

}

void TopologyFile::ResolveTitleSection(ofstream &out)
{
    out << "%FLAG TITLE" << endl
        << "%FORMAT(20a4)" << endl
        << title_ << endl;
}

void TopologyFile::ResolvePointersSection(ofstream &out)
{
    out << "%FLAG POINTERS" << endl
        << "%FORMAT(10I8)" << endl;
    if(number_of_atoms_ != iNotSet)
        out << setw(8) << right << number_of_atoms_;
    else
        out << setw(8) << right << 0;
    if(number_of_types_ != iNotSet)
        out << setw(8) << right << number_of_types_;
    else
        out << setw(8) << right << 0;
    if(number_of_bonds_including_hydrogen_ != iNotSet)
        out << setw(8) << right << number_of_bonds_including_hydrogen_;
    else
        out << setw(8) << right << 0;
    if(number_of_bonds_excluding_hydrogen_ != iNotSet)
        out << setw(8) << right << number_of_bonds_excluding_hydrogen_;
    else
        out << setw(8) << right << 0;
    if(number_of_angles_including_hydrogen_ != iNotSet)
        out << setw(8) << right << number_of_angles_including_hydrogen_;
    else
        out << setw(8) << right << 0;
    if(number_of_angles_excluding_hydrogen_ != iNotSet)
        out << setw(8) << right << number_of_angles_excluding_hydrogen_;
    else
        out << setw(8) << right << 0;
    if(number_of_dihedrals_including_hydrogen_ != iNotSet)
        out << setw(8) << right << number_of_dihedrals_including_hydrogen_;
    else
        out << setw(8) << right << 0;
    if(number_of_dihedrals_excluding_hydrogen_ != iNotSet)
        out << setw(8) << right << number_of_dihedrals_excluding_hydrogen_;
    else
        out << setw(8) << right << 0;
    if(number_of_hydrogen_parameters_ != iNotSet)
        out << setw(8) << right << number_of_hydrogen_parameters_;
    else
        out << setw(8) << right << 0;
    if(number_of_parameters_ != iNotSet)
        out << setw(8) << right << number_of_parameters_;
    else
        out << setw(8) << right << 0;
    out << endl;
    if(number_of_excluded_atoms_ != iNotSet)
        out << setw(8) << right << number_of_excluded_atoms_;
    else
        out << setw(8) << right << 0;
    if(number_of_residues_ != iNotSet)
        out << setw(8) << right << number_of_residues_;
    else
        out << setw(8) << right << 0;
    if(total_number_of_bonds_ != iNotSet)
        out << setw(8) << right << total_number_of_bonds_;
    else
        out << setw(8) << right << 0;
    if(total_number_of_angles_ != iNotSet)
        out << setw(8) << right << total_number_of_angles_;
    else
        out << setw(8) << right << 0;
    if(total_number_of_dihedrals_ != iNotSet)
        out << setw(8) << right << total_number_of_dihedrals_;
    else
        out << setw(8) << right << 0;
    if(number_of_bond_types_ != iNotSet)
        out << setw(8) << right << number_of_bond_types_;
    else
        out << setw(8) << right << 0;
    if(number_of_angle_types_ != iNotSet)
        out << setw(8) << right << number_of_angle_types_;
    else
        out << setw(8) << right << 0;
    if(number_of_dihedral_types_ != iNotSet)
        out << setw(8) << right << number_of_dihedral_types_;
    else
        out << setw(8) << right << 0;
    if(number_of_atom_types_in_parameter_file_ != iNotSet)
        out << setw(8) << right << number_of_atom_types_in_parameter_file_;
    else
        out << setw(8) << right << 0;
    if(number_of_distinct_hydrogen_bonds_ != iNotSet)
        out << setw(8) << right << number_of_distinct_hydrogen_bonds_;
    else
        out << setw(8) << right << 0;
    out << endl;
    if(perturbation_option_ != iNotSet)
        out << setw(8) << right << perturbation_option_;
    else
        out << setw(8) << right << 0;
    if(number_of_bonds_perturbed_ != iNotSet)
        out << setw(8) << right << number_of_bonds_perturbed_;
    else
        out << setw(8) << right << 0;
    if(number_of_angles_perturbed_ != iNotSet)
        out << setw(8) << right << number_of_angles_perturbed_;
    else
        out << setw(8) << right << 0;
    if(number_of_dihedrals_perturbed_ != iNotSet)
        out << setw(8) << right << number_of_dihedrals_perturbed_;
    else
        out << setw(8) << right << 0;
    if(number_of_bonds_group_perturbed_ != iNotSet)
        out << setw(8) << right << number_of_bonds_group_perturbed_;
    else
        out << setw(8) << right << 0;
    if(number_of_angles_group_perturbed_ != iNotSet)
        out << setw(8) << right << number_of_angles_group_perturbed_;
    else
        out << setw(8) << right << 0;
    if(number_of_dihedrals_group_perturbed_ != iNotSet)
        out << setw(8) << right << number_of_dihedrals_group_perturbed_;
    else
        out << setw(8) << right << 0;
    if(standard_periodic_box_option_ != iNotSet)
        out << setw(8) << right << standard_periodic_box_option_;
    else
        out << setw(8) << right << 0;
    if(number_of_atoms_in_largest_residue_ != iNotSet)
        out << setw(8) << right << number_of_atoms_in_largest_residue_;
    else
        out << setw(8) << right << 0;
    if(cap_option_ != iNotSet)
        out << setw(8) << right << cap_option_;
    else
        out << setw(8) << right << 0;
    out << endl;
    if(number_of_extra_points_ != iNotSet)
        out << setw(8) << number_of_extra_points_;
    else
        out << setw(8) << right << 0;
    if(number_of_beads_ != iNotSet)
        out << setw(8) << number_of_beads_ << endl;
    else
        out << setw(8) << right << 0 << endl;
}

void TopologyFile::ResolveAtomNameSection(ofstream& out)
{
    out << "%FLAG ATOM_NAME" << endl
        << "%FORMAT(20a4)" << endl;
    int count = 0;
    int index = 0;
    const int MAX_IN_LINE = 20;
    const int ITEM_LENGTH = 4;
    for(int i = 0; i < number_of_residues_; i++)
    {
        TopologyResidue* residue = this->assembly_->GetResidueByIndex(i+1);
        for(unsigned int j = 0; j < residue->GetAtoms().size(); j++)
        {
            TopologyAtom* atom = residue->GetAtomByIndex(index+1);
            index++;
            out << setw(ITEM_LENGTH) << left << atom->GetAtomName();
            count++;
            if(count == MAX_IN_LINE)
            {
                count = 0;
                out << endl;
            }
        }
    }
    if(count < MAX_IN_LINE)
        out << endl;
}

void TopologyFile::ResolveChargeSection(ofstream& out)
{
    out << "%FLAG CHARGE" << endl
        << "%FORMAT(5E16.8)" << endl;
    int count = 0;
    int index = 0;
    const int MAX_IN_LINE = 5;
    const int ITEM_LENGTH = 16;
    for(int i = 0; i < number_of_residues_; i++)
    {
        TopologyResidue* residue = this->assembly_->GetResidueByIndex(i+1);
        for(unsigned int j = 0; j < residue->GetAtoms().size(); j++)
        {
            TopologyAtom* atom = residue->GetAtomByIndex(index+1);
            index++;
            out << setw(ITEM_LENGTH) << right << scientific << setprecision(8) << atom->GetAtomCharge();
            count++;
            if(count == MAX_IN_LINE)
            {
                count = 0;
                out << endl;
            }
        }
    }
    if(count < MAX_IN_LINE)
        out << endl;
}

void TopologyFile::ResolveAtomicNumberSection(ofstream& out)
{
    out << "%FLAG ATOMIC_NUMBER" << endl
        << "%FORMAT(10I8)" << endl;
    int count = 0;
    int index = 0;
    const int MAX_IN_LINE = 10;
    const int ITEM_LENGTH = 8;
    for(int i = 0; i < number_of_residues_; i++)
    {
        TopologyResidue* residue = this->assembly_->GetResidueByIndex(i+1);
        for(unsigned int j = 0; j < residue->GetAtoms().size(); j++)
        {
            TopologyAtom* atom = residue->GetAtomByIndex(index+1);
            index++;
            if(atom->GetAtomicNumber() != iNotSet)
                out << setw(ITEM_LENGTH) << right << atom->GetAtomicNumber();
            else
                out << setw(ITEM_LENGTH) << right << 0;
            count++;
            if(count == MAX_IN_LINE)
            {
                count = 0;
                out << endl;
            }
        }
    }
    if(count < MAX_IN_LINE)
        out << endl;
}

void TopologyFile::ResolveMassSection(ofstream& out)
{
    out << "%FLAG MASS" << endl
        << "%FORMAT(5E16.8)" << endl;
    int count = 0;
    int index = 0;
    const int MAX_IN_LINE = 5;
    const int ITEM_LENGTH = 16;
    for(int i = 0; i < number_of_residues_; i++)
    {
        TopologyResidue* residue = this->assembly_->GetResidueByIndex(i+1);
        for(unsigned int j = 0; j < residue->GetAtoms().size(); j++)
        {
            TopologyAtom* atom = residue->GetAtomByIndex(index+1);
            index++;
            out << setw(ITEM_LENGTH) << right << scientific << setprecision(8) << atom->GetAtomMass();
            count++;
            if(count == MAX_IN_LINE)
            {
                count = 0;
                out << endl;
            }
        }
    }
    if(count < MAX_IN_LINE)
        out << endl;
}

void TopologyFile::ResolveAtomTypeIndexSection(ofstream& out)
{
    out << "%FLAG ATOM_TYPE_INDEX" << endl
        << "%FORMAT(10I8)" << endl;
    int count = 0;
    int index = 0;
    const int MAX_IN_LINE = 10;
    const int ITEM_LENGTH = 8;
    TopologyAtomTypeIndexMap atom_type_index_map = this->GetAtomTypeIndexMap();
    for(int i = 0; i < number_of_residues_; i++)
    {
        TopologyResidue* residue = this->assembly_->GetResidueByIndex(i+1);
        for(unsigned int j = 0; j < residue->GetAtoms().size(); j++)
        {
            TopologyAtom* atom = residue->GetAtomByIndex(index+1);
            index++;
            out << setw(ITEM_LENGTH) << right << atom_type_index_map[atom->GetType()];
            count++;
            if(count == MAX_IN_LINE)
            {
                count = 0;
                out << endl;
            }
        }
    }
    if(count < MAX_IN_LINE)
        out << endl;
}

void TopologyFile::ResolveNumberExcludedAtomsSection(ofstream& out)
{
    out << "%FLAG NUMBER_EXCLUDED_ATOMS" << endl
        << "%FORMAT(10I8)" << endl;
    int count = 0;
    int index = 0;
    const int MAX_IN_LINE = 10;
    const int ITEM_LENGTH = 8;
    for(int i = 0; i < number_of_residues_; i++)
    {
        TopologyResidue* residue = this->assembly_->GetResidueByIndex(i+1);
        for(unsigned int j = 0; j < residue->GetAtoms().size(); j++)
        {
            TopologyAtom* atom = residue->GetAtomByIndex(index+1);
            index++;
            out << setw(ITEM_LENGTH) << right << atom->GetExcludedAtoms().size();
            count++;
            if(count == MAX_IN_LINE)
            {
                count = 0;
                out << endl;
            }
        }
    }
    if(count < MAX_IN_LINE)
        out << endl;
}

void TopologyFile::ResolveNonbondedParmIndexSection(ofstream& out)
{
    out << "%FLAG NONBONDED_PARM_INDEX" << endl
        << "%FORMAT(10I8)" << endl;
    int index = 0;
    int count = 0;
    const int MAX_IN_LINE = 10;
    const int ITEM_LENGTH = 8;

    TopologyAtomTypeIndexMap atom_type_index_map = this->GetAtomTypeIndexMap();
    vector<string> atom_types = vector<string>();
    for(int i = 0; i < number_of_residues_; i++)
    {
        TopologyResidue* residue = this->assembly_->GetResidueByIndex(i+1);
        for(unsigned int j = 0; j < residue->GetAtoms().size(); j++)
        {
            TopologyAtom* atom = residue->GetAtomByIndex(index+1);
            index++;
            atom_types.push_back(atom->GetType());
        }
    }
    map<string, bool> pair_set = map<string, bool>();
    map<int, int> nonbond_index_map = map<int, int>();
    for(int i = 0; i < number_of_atoms_; i++)
    {
        string atom_type_1 = atom_types.at(i);
        for(int j = 0; j < number_of_atoms_; j++)
        {
            string atom_type_2 = atom_types.at(j);
            if(!pair_set[atom_type_1 + "-" + atom_type_2])
            {
                int nonbond_index = number_of_types_ * (atom_type_index_map[atom_type_1] - 1) + atom_type_index_map[atom_type_2];
                int parameter_index = 0;
                if(pairs_.find(atom_type_1 + "-" + atom_type_2) != pairs_.end())
                {
                    parameter_index = pairs_[atom_type_1 + "-" + atom_type_2]->GetIndex();
                }
                else if(pairs_.find(atom_type_2 + "-" + atom_type_1) != pairs_.end())
                {
                    parameter_index = pairs_[atom_type_2 + "-" + atom_type_1]->GetIndex();
                }
                nonbond_index_map[nonbond_index] = parameter_index;
                pair_set[atom_type_1 + "-" + atom_type_2] = true;
            }
        }
    }
    int size = nonbond_index_map.size();
    for(unsigned int i = 0; i < size; i++)
    {
        out << setw(ITEM_LENGTH) << right << nonbond_index_map[i+1];
        count++;
        if(count == MAX_IN_LINE)
        {
            count = 0;
            out << endl;
        }
    }
    if(count < MAX_IN_LINE)
        out << endl;
}

void TopologyFile::ResolveResidueLabelSection(ofstream& out)
{
    out << "%FLAG RESIDUE_LABEL" << endl
        << "%FORMAT(20a4)" << endl;
    int count = 0;
    const int MAX_IN_LINE = 20;
    const int ITEM_LENGTH = 4;
    for(int i = 0; i < number_of_residues_; i++)
    {
        TopologyResidue* residue = this->assembly_->GetResidueByIndex(i+1);

        out << setw(ITEM_LENGTH) << left << residue->GetResidueName();
        count++;
        if(count == MAX_IN_LINE)
        {
            count = 0;
            out << endl;
        }
    }
    if(count < MAX_IN_LINE)
        out << endl;
}

void TopologyFile::ResolveResiduePointersSection(ofstream& out)
{
    out << "%FLAG RESIDUE_POINTER" << endl
        << "%FORMAT(10I8)" << endl;
    int count = 0;
    const int MAX_IN_LINE = 10;
    const int ITEM_LENGTH = 8;
    for(int i = 0; i < number_of_residues_; i++)
    {
        TopologyResidue* residue = this->assembly_->GetResidueByIndex(i+1);

        out << setw(ITEM_LENGTH) << right << residue->GetStartingAtomIndex();
        count++;
        if(count == MAX_IN_LINE)
        {
            count = 0;
            out << endl;
        }
    }
    if(count < MAX_IN_LINE)
        out << endl;
}

void TopologyFile::ResolveBondForceConstantSection(ofstream& out)
{
    out << "%FLAG BOND_FORCE_CONSTANT" << endl
        << "%FORMAT(5E16.8)" << endl;
    int count = 0;
    const int MAX_IN_LINE = 5;
    const int ITEM_LENGTH = 16;
    for(int i = 0; i < number_of_bond_types_; i++)
    {
        TopologyBondType* bond_type = this->GetBondTypeByIndex(i);
        if(bond_type != NULL)
        {
            out << setw(ITEM_LENGTH) << right << scientific << setprecision(8) << bond_type->GetForceConstant();
            count++;
            if(count == MAX_IN_LINE)
            {
                count = 0;
                out << endl;
            }
        }
    }
    if(count < MAX_IN_LINE)
        out << endl;
}

void TopologyFile::ResolveBondEquilValueSection(ofstream& out)
{
    out << "%FLAG BOND_EQUIL_VALUE" << endl
        << "%FORMAT(5E16.8)" << endl;
    int count = 0;
    const int MAX_IN_LINE = 5;
    const int ITEM_LENGTH = 16;
    for(int i = 0; i < number_of_bond_types_; i++)
    {
        TopologyBondType* bond_type = this->GetBondTypeByIndex(i);
        if(bond_type != NULL)
        {
            out << setw(ITEM_LENGTH) << right << scientific << setprecision(8) << bond_type->GetEquilibriumValue();
            count++;
            if(count == MAX_IN_LINE)
            {
                count = 0;
                out << endl;
            }
        }
    }
    if(count < MAX_IN_LINE)
        out << endl;
}

void TopologyFile::ResolveAngleForceConstantSection(ofstream& out)
{
    out << "%FLAG ANGLE_FORCE_CONSTANT" << endl
        << "%FORMAT(5E16.8)" << endl;
    int count = 0;
    const int MAX_IN_LINE = 5;
    const int ITEM_LENGTH = 16;
    for(int i = 0; i < number_of_angle_types_; i++)
    {
        TopologyAngleType* angle_type = this->GetAngleTypeByIndex(i);
        if(angle_type != NULL)
        {
            out << setw(ITEM_LENGTH) << right << scientific << setprecision(8) << angle_type->GetForceConstant();
            count++;
            if(count == MAX_IN_LINE)
            {
                count = 0;
                out << endl;
            }
        }
    }
    if(count < MAX_IN_LINE)
        out << endl;
}

void TopologyFile::ResolveAngleEquilValueSection(ofstream& out)
{
    out << "%FLAG ANGLE_EQUIL_VALUE" << endl
        << "%FORMAT(5E16.8)" << endl;
    int count = 0;
    const int MAX_IN_LINE = 5;
    const int ITEM_LENGTH = 16;
    for(int i = 0; i < number_of_angle_types_; i++)
    {
        TopologyAngleType* angle_type = this->GetAngleTypeByIndex(i);
        if(angle_type != NULL)
        {
            out << setw(ITEM_LENGTH) << right << scientific << setprecision(8) << angle_type->GetEquilibriumValue();
            count++;
            if(count == MAX_IN_LINE)
            {
                count = 0;
                out << endl;
            }
        }
    }
    if(count < MAX_IN_LINE)
        out << endl;
}

void TopologyFile::ResolveDihedralForceConstantSection(ofstream& out)
{
    out << "%FLAG DIHEDRAL_FORCE_CONSTANT" << endl
        << "%FORMAT(5E16.8)" << endl;
    int count = 0;
    const int MAX_IN_LINE = 5;
    const int ITEM_LENGTH = 16;
    for(int i = 0; i < number_of_dihedral_types_; i++)
    {
        TopologyDihedralType* dihedral_type = this->GetDihedralTypeByIndex(i);
        if(dihedral_type != NULL)
        {
            out << setw(ITEM_LENGTH) << right << scientific << setprecision(8) << dihedral_type->GetForceConstant();
            count++;
            if(count == MAX_IN_LINE)
            {
                count = 0;
                out << endl;
            }
        }
    }
    if(count < MAX_IN_LINE)
        out << endl;
}

void TopologyFile::ResolveDihedralPeriodicitySection(ofstream& out)
{
    out << "%FLAG DIHEDRAL_PERIODICITY" << endl
        << "%FORMAT(5E16.8)" << endl;
    int count = 0;
    const int MAX_IN_LINE = 5;
    const int ITEM_LENGTH = 16;
    for(int i = 0; i < number_of_dihedral_types_; i++)
    {
        TopologyDihedralType* dihedral_type = this->GetDihedralTypeByIndex(i);
        if(dihedral_type != NULL)
        {
            out << setw(ITEM_LENGTH) << right << scientific << setprecision(8) << dihedral_type->GetPeriodicity();
            count++;
            if(count == MAX_IN_LINE)
            {
                count = 0;
                out << endl;
            }
        }
    }
    if(count < MAX_IN_LINE)
        out << endl;
}

void TopologyFile::ResolveDihedralPhaseSection(ofstream& out)
{
    out << "%FLAG DIHEDRAL_PHASE" << endl
        << "%FORMAT(5E16.8)" << endl;
    int count = 0;
    const int MAX_IN_LINE = 5;
    const int ITEM_LENGTH = 16;
    for(int i = 0; i < number_of_dihedral_types_; i++)
    {
        TopologyDihedralType* dihedral_type = this->GetDihedralTypeByIndex(i);
        if(dihedral_type != NULL)
        {
            out << setw(ITEM_LENGTH) << right << scientific << setprecision(8) << dihedral_type->GetPhase();
            count++;
            if(count == MAX_IN_LINE)
            {
                count = 0;
                out << endl;
            }
        }
    }
    if(count < MAX_IN_LINE)
        out << endl;
}

void TopologyFile::ResolveSceeScaleFactorSection(ofstream& out)
{
    out << "%FLAG SCEE_SCALE_FACTOR" << endl
        << "%FORMAT(5E16.8)" << endl;
    int count = 0;
    const int MAX_IN_LINE = 5;
    const int ITEM_LENGTH = 16;
    for(int i = 0; i < number_of_dihedral_types_; i++)
    {
        TopologyDihedralType* dihedral_type = this->GetDihedralTypeByIndex(i);
        if(dihedral_type != NULL)
        {
            out << setw(ITEM_LENGTH) << right << scientific << setprecision(8) << dihedral_type->GetScee();
            count++;
            if(count == MAX_IN_LINE)
            {
                count = 0;
                out << endl;
            }
        }
    }
    if(count < MAX_IN_LINE)
        out << endl;
}

void TopologyFile::ResolveScnbScaleFactorSection(ofstream& out)
{
    out << "%FLAG SCNB_SCALE_FACTOR" << endl
        << "%FORMAT(5E16.8)" << endl;
    int count = 0;
    const int MAX_IN_LINE = 5;
    const int ITEM_LENGTH = 16;
    for(int i = 0; i < number_of_dihedral_types_; i++)
    {
        TopologyDihedralType* dihedral_type = this->GetDihedralTypeByIndex(i);
        if(dihedral_type != NULL)
        {
            out << setw(ITEM_LENGTH) << right << scientific << setprecision(8) << dihedral_type->GetScnb();
            count++;
            if(count == MAX_IN_LINE)
            {
                count = 0;
                out << endl;
            }
        }
    }
    if(count < MAX_IN_LINE)
        out << endl;
}

void TopologyFile::ResolveSoltySection(ofstream& out)
{
    out << "%FLAG SOLTY" << endl
        << "%FORMAT(5E16.8)" << endl;
    //
}

void TopologyFile::ResolveLennardJonesACoefSection(ofstream& out)
{
    out << "%FLAG LENNARD_JONES_ACOEF" << endl
        << "%FORMAT(5E16.8)" << endl;
    int count = 0;
    const int MAX_IN_LINE = 5;
    const int ITEM_LENGTH = 16;
    for(unsigned int i = 0; i < pairs_.size(); i++)
    {
        TopologyAtomPair* atom_pair = this->GetAtomPairByIndex(i+1);
        if(atom_pair != NULL)
        {
            out << setw(ITEM_LENGTH) << right << scientific << setprecision(8) << atom_pair->GetCoefficientA();
            count++;
            if(count == MAX_IN_LINE)
            {
                count = 0;
                out << endl;
            }
        }
    }

    if(count < MAX_IN_LINE)
        out << endl;
}

void TopologyFile::ResolveLennardJonesBCoefSection(ofstream& out)
{
    out << "%FLAG LENNARD_JONES_BCOEF" << endl
        << "%FORMAT(5E16.8)" << endl;
    int count = 0;
    const int MAX_IN_LINE = 5;
    const int ITEM_LENGTH = 16;
    for(unsigned int i = 0; i < pairs_.size(); i++)
    {
        TopologyAtomPair* atom_pair = this->GetAtomPairByIndex(i+1);
        if(atom_pair != NULL)
        {
            out << setw(ITEM_LENGTH) << right << scientific << setprecision(8) << atom_pair->GetCoefficientB();
            count++;
            if(count == MAX_IN_LINE)
            {
                count = 0;
                out << endl;
            }
        }
    }

    if(count < MAX_IN_LINE)
        out << endl;
}

void TopologyFile::ResolveBondsIncHydrogenSection(ofstream& out)
{
    out << "%FLAG BONDS_INC_HYDROGEN" << endl
        << "%FORMAT(10I8)" << endl;
    int count = 0;
    const int MAX_IN_LINE = 10;
    const int ITEM_LENGTH = 8;
    for(TopologyBondMap::iterator it = bonds_.begin(); it != bonds_.end(); it++)
    {
        TopologyBond* bond = (*it).second;
        vector<string> atom_names = bond->GetBonds();
        if(bond->GetIncludingHydrogen())
        {
            int atom_index_1 = this->assembly_->GetAtomIndexByName(atom_names.at(0));
            int atom_index_2 = this->assembly_->GetAtomIndexByName(atom_names.at(1));
            int bond_type_index = bond->GetBondType()->GetIndex();
            out << setw(ITEM_LENGTH) << right << (atom_index_1-1)*3;
            count++;
            if(count == MAX_IN_LINE)
            {
                count = 0;
                out << endl;
            }
            out << setw(ITEM_LENGTH) << right << (atom_index_2-1)*3;
            count++;
            if(count == MAX_IN_LINE)
            {
                count = 0;
                out << endl;
            }
            out << setw(ITEM_LENGTH) << right << bond_type_index+1;
            count++;
            if(count == MAX_IN_LINE)
            {
                count = 0;
                out << endl;
            }
        }
    }
    if(count < MAX_IN_LINE)
        out << endl;
}

void TopologyFile::ResolveBondsWithoutHydrogenSection(ofstream& out)
{
    out << "%FLAG BONDS_WITHOUT_HYDROGEN" << endl
        << "%FORMAT(10I8)" << endl;
    int count = 0;
    const int MAX_IN_LINE = 10;
    const int ITEM_LENGTH = 8;
    for(TopologyBondMap::iterator it = bonds_.begin(); it != bonds_.end(); it++)
    {
        TopologyBond* bond = (*it).second;
        vector<string> atom_names = bond->GetBonds();
        if(!bond->GetIncludingHydrogen())
        {
            int atom_index_1 = this->assembly_->GetAtomIndexByName(atom_names.at(0));
            int atom_index_2 = this->assembly_->GetAtomIndexByName(atom_names.at(1));
            int bond_type_index = bond->GetBondType()->GetIndex();
            out << setw(ITEM_LENGTH) << right << (atom_index_1-1)*3;
            count++;
            if(count == MAX_IN_LINE)
            {
                count = 0;
                out << endl;
            }
            out << setw(ITEM_LENGTH) << right << (atom_index_2-1)*3;
            count++;
            if(count == MAX_IN_LINE)
            {
                count = 0;
                out << endl;
            }
            out << setw(ITEM_LENGTH) << right << bond_type_index+1;
            count++;
            if(count == MAX_IN_LINE)
            {
                count = 0;
                out << endl;
            }
        }
    }
    if(count < MAX_IN_LINE)
        out << endl;
}

void TopologyFile::ResolveAnglesIncHydrogenSection(ofstream& out)
{
    out << "%FLAG ANGLES_INC_HYDROGEN" << endl
        << "%FORMAT(10I8)" << endl;
    int count = 0;
    const int MAX_IN_LINE = 10;
    const int ITEM_LENGTH = 8;
    for(TopologyAngleMap::iterator it = angles_.begin(); it != angles_.end(); it++)
    {
        TopologyAngle* angle = (*it).second;
        vector<string> atom_names = angle->GetAngles();
        if(angle->GetIncludingHydrogen())
        {
            int atom_index_1 = this->assembly_->GetAtomIndexByName(atom_names.at(0));
            int atom_index_2 = this->assembly_->GetAtomIndexByName(atom_names.at(1));
            int atom_index_3 = this->assembly_->GetAtomIndexByName(atom_names.at(2));
            int angle_type_index = angle->GetAngleType()->GetIndex();
            out << setw(ITEM_LENGTH) << right << (atom_index_1-1)*3;
            count++;
            if(count == MAX_IN_LINE)
            {
                count = 0;
                out << endl;
            }
            out << setw(ITEM_LENGTH) << right << (atom_index_2-1)*3;
            count++;
            if(count == MAX_IN_LINE)
            {
                count = 0;
                out << endl;
            }
            out << setw(ITEM_LENGTH) << right << (atom_index_3-1)*3;
            count++;
            if(count == MAX_IN_LINE)
            {
                count = 0;
                out << endl;
            }
            out << setw(ITEM_LENGTH) << right << angle_type_index+1;
            count++;
            if(count == MAX_IN_LINE)
            {
                count = 0;
                out << endl;
            }
        }
    }
    if(count < MAX_IN_LINE)
        out << endl;
}

void TopologyFile::ResolveAnglesWithoutHydrogenSection(ofstream& out)
{
    out << "%FLAG ANGLES_WITHOUT_HYDROGEN" << endl
        << "%FORMAT(10I8)" << endl;
    int count = 0;
    const int MAX_IN_LINE = 10;
    const int ITEM_LENGTH = 8;
    for(TopologyAngleMap::iterator it = angles_.begin(); it != angles_.end(); it++)
    {
        TopologyAngle* angle = (*it).second;
        vector<string> atom_names = angle->GetAngles();
        if(!angle->GetIncludingHydrogen())
        {
            int atom_index_1 = this->assembly_->GetAtomIndexByName(atom_names.at(0));
            int atom_index_2 = this->assembly_->GetAtomIndexByName(atom_names.at(1));
            int atom_index_3 = this->assembly_->GetAtomIndexByName(atom_names.at(2));
            int angle_type_index = angle->GetAngleType()->GetIndex();
            out << setw(ITEM_LENGTH) << right << (atom_index_1-1)*3;
            count++;
            if(count == MAX_IN_LINE)
            {
                count = 0;
                out << endl;
            }
            out << setw(ITEM_LENGTH) << right << (atom_index_2-1)*3;
            count++;
            if(count == MAX_IN_LINE)
            {
                count = 0;
                out << endl;
            }
            out << setw(ITEM_LENGTH) << right << (atom_index_3-1)*3;
            count++;
            if(count == MAX_IN_LINE)
            {
                count = 0;
                out << endl;
            }
            out << setw(ITEM_LENGTH) << right << angle_type_index+1;
            count++;
            if(count == MAX_IN_LINE)
            {
                count = 0;
                out << endl;
            }
        }
    }
    if(count < MAX_IN_LINE)
        out << endl;
}

void TopologyFile::ResolveDihedralsIncHydrogenSection(ofstream& out)
{
    out << "%FLAG DIHEDRALS_INC_HYDROGEN" << endl
        << "%FORMAT(10I8)" << endl;
    int count = 0;
    const int MAX_IN_LINE = 10;
    const int ITEM_LENGTH = 8;
    for(TopologyDihedralMap::iterator it = dihedrals_.begin(); it != dihedrals_.end(); it++)
    {
        TopologyDihedral* dihedral = (*it).second;
        vector<string> atom_names = dihedral->GetDihedrals();
        if(dihedral->GetIncludingHydrogen())
        {
            int atom_index_1 = this->assembly_->GetAtomIndexByName(atom_names.at(0));
            int atom_index_2 = this->assembly_->GetAtomIndexByName(atom_names.at(1));
            int atom_index_3 = this->assembly_->GetAtomIndexByName(atom_names.at(2));
            int atom_index_4 = this->assembly_->GetAtomIndexByName(atom_names.at(3));
            int dihedral_type_index = dihedral->GetDihedralType()->GetIndex();
            out << setw(ITEM_LENGTH) << right << (atom_index_1-1)*3;
            count++;
            if(count == MAX_IN_LINE)
            {
                count = 0;
                out << endl;
            }
            out << setw(ITEM_LENGTH) << right << (atom_index_2-1)*3;
            count++;
            if(count == MAX_IN_LINE)
            {
                count = 0;
                out << endl;
            }
            if(dihedral->GetIgnoredGroupInteraction())
                out << setw(ITEM_LENGTH) << right << -(atom_index_3-1)*3;
            else
                out << setw(ITEM_LENGTH) << right << (atom_index_3-1)*3;
            count++;
            if(count == MAX_IN_LINE)
            {
                count = 0;
                out << endl;
            }
            if(dihedral->GetIsImproper())
                out << setw(ITEM_LENGTH) << right << -(atom_index_4-1)*3;
            else
                out << setw(ITEM_LENGTH) << right << (atom_index_4-1)*3;
            count++;
            if(count == MAX_IN_LINE)
            {
                count = 0;
                out << endl;
            }
            out << setw(ITEM_LENGTH) << right << dihedral_type_index+1;
            count++;
            if(count == MAX_IN_LINE)
            {
                count = 0;
                out << endl;
            }
        }
    }
    if(count < MAX_IN_LINE)
        out << endl;
}

void TopologyFile::ResolveDihedralsWithoutHydrogenSection(ofstream& out)
{
    out << "%FLAG DIHEDRALS_WITHOUT_HYDROGEN" << endl
        << "%FORMAT(10I8)" << endl;
    int count = 0;
    const int MAX_IN_LINE = 10;
    const int ITEM_LENGTH = 8;
    for(TopologyDihedralMap::iterator it = dihedrals_.begin(); it != dihedrals_.end(); it++)
    {
        TopologyDihedral* dihedral = (*it).second;
        vector<string> atom_names = dihedral->GetDihedrals();
        if(!dihedral->GetIncludingHydrogen())
        {
            int atom_index_1 = this->assembly_->GetAtomIndexByName(atom_names.at(0));
            int atom_index_2 = this->assembly_->GetAtomIndexByName(atom_names.at(1));
            int atom_index_3 = this->assembly_->GetAtomIndexByName(atom_names.at(2));
            int atom_index_4 = this->assembly_->GetAtomIndexByName(atom_names.at(3));
            int dihedral_type_index = dihedral->GetDihedralType()->GetIndex();
            out << setw(ITEM_LENGTH) << right << (atom_index_1-1)*3;
            count++;
            if(count == MAX_IN_LINE)
            {
                count = 0;
                out << endl;
            }
            out << setw(ITEM_LENGTH) << right << (atom_index_2-1)*3;
            count++;
            if(count == MAX_IN_LINE)
            {
                count = 0;
                out << endl;
            }
            if(dihedral->GetIgnoredGroupInteraction())
                out << setw(ITEM_LENGTH) << right << -(atom_index_3-1)*3;
            else
                out << setw(ITEM_LENGTH) << right << (atom_index_3-1)*3;
            count++;
            if(count == MAX_IN_LINE)
            {
                count = 0;
                out << endl;
            }
            if(dihedral->GetIsImproper())
                out << setw(ITEM_LENGTH) << right << -(atom_index_4-1)*3;
            else
                out << setw(ITEM_LENGTH) << right << (atom_index_4-1)*3;
            count++;
            if(count == MAX_IN_LINE)
            {
                count = 0;
                out << endl;
            }
            out << setw(ITEM_LENGTH) << right << dihedral_type_index+1;
            count++;
            if(count == MAX_IN_LINE)
            {
                count = 0;
                out << endl;
            }
        }
    }
    if(count < MAX_IN_LINE)
        out << endl;
}

void TopologyFile::ResolveExcludedAtomsListSection(ofstream& out)
{
    out << "%FLAG EXCLUDED_ATOM_LIST" << endl
        << "%FORMAT(10I8)" << endl;
    int count = 0;
    int index = 0;
    const int MAX_IN_LINE = 10;
    const int ITEM_LENGTH = 8;
    for(int i = 0; i < number_of_residues_; i++)
    {
        TopologyResidue* residue = this->assembly_->GetResidueByIndex(i+1);
        for(unsigned int j = 0; j < residue->GetAtoms().size(); j++)
        {
            TopologyAtom* atom = residue->GetAtomByIndex(index+1);
            index++;
            vector<string> excluded_atoms = atom->GetExcludedAtoms();
            for(vector<string>::iterator it = excluded_atoms.begin(); it != excluded_atoms.end(); it++)
            {
                string atom_id = (*it);
                int index = this->assembly_->GetAtomIndexByName(Split(atom_id, ":").at(1));
                out << setw(ITEM_LENGTH) << right << index;
                count++;
                if(count == MAX_IN_LINE)
                {
                    count = 0;
                    out << endl;
                }
            }
        }
    }
    if(count < MAX_IN_LINE)
        out << endl;
}

void TopologyFile::ResolveHydrogenBondACoefSection(ofstream& out)
{
    out << "%FLAG HBOND_ACOEF" << endl
        << "%FORMAT(5E16.8)" << endl;
    //
}

void TopologyFile::ResolveHydrogenBondBCoefSection(ofstream& out)
{
    out << "%FLAG HBOND_BCOEF" << endl
        << "%FORMAT(5E16.8)" << endl;
    //
}

void TopologyFile::ResolveHBCutSection(ofstream& out)
{
    out << "%FLAG HBCUT" << endl
        << "%FORMAT(5E16.8)" << endl;
    //
}

void TopologyFile::ResolveAmberAtomTypeSection(ofstream& out)
{
    out << "%FLAG AMBER_ATOM_TYPE" << endl
        << "%FORMAT(20a4)" << endl;
    int count = 0;
    int index = 0;
    const int MAX_IN_LINE = 20;
    const int ITEM_LENGTH = 4;
    for(int i = 0; i < number_of_residues_; i++)
    {
        TopologyResidue* residue = this->assembly_->GetResidueByIndex(i+1);
        for(unsigned int j = 0; j < residue->GetAtoms().size(); j++)
        {
            TopologyAtom* atom = residue->GetAtomByIndex(index+1);
            index++;
            out << setw(ITEM_LENGTH) << left << atom->GetType();
            count++;
            if(count == MAX_IN_LINE)
            {
                count = 0;
                out << endl;
            }
        }
    }
    if(count < MAX_IN_LINE)
        out << endl;
}

void TopologyFile::ResolveTreeChainClassificationSection(ofstream& out)
{
    out << "%FLAG TREE_CHAIN_CLASSIFICATION" << endl
        << "%FORMAT(20a4)" << endl;
    int count = 0;
    int index = 0;
    const int MAX_IN_LINE = 20;
    const int ITEM_LENGTH = 4;
    for(int i = 0; i < number_of_residues_; i++)
    {
        TopologyResidue* residue = this->assembly_->GetResidueByIndex(i+1);
        for(unsigned int j = 0; j < residue->GetAtoms().size(); j++)
        {
            TopologyAtom* atom = residue->GetAtomByIndex(index+1);
            index++;
            if(!atom->GetTreeChainClassification().empty())
                out << setw(ITEM_LENGTH) << left << atom->GetTreeChainClassification();
            else
                out << setw(ITEM_LENGTH) << left << "0";
            count++;
            if(count == MAX_IN_LINE)
            {
                count = 0;
                out << endl;
            }
        }
    }
    if(count < MAX_IN_LINE)
        out << endl;
}

void TopologyFile::ResolveJoinArraySection(ofstream& out)
{
    out << "%FLAG JOIN_ARRAY" << endl
        << "%FORMAT(10I8)" << endl;
    //
    out << endl;
}

void TopologyFile::ResolveIRotatSection(ofstream& out)
{
    out << "%FLAG IROTAT" << endl
        << "%FORMAT(10I8)" << endl;
    //
}

void TopologyFile::ResolveRadiusSetSection(ofstream& out)
{
    if(radius_set_.size() != 0)
    {
        string radius_set = radius_set_.at(0);
        std::replace(radius_set.begin(), radius_set.end(), '#', ' ');
        out << "%FLAG RADIUS_SET" << endl
            << "%FORMAT(1a80)" << endl
            << radius_set << endl;
    }
    else
    {
        out << "%FLAG RADIUS_SET" << endl
            << "%FORMAT(1a80)" << endl;
    }

}

void TopologyFile::ResolveRadiiSection(ofstream& out)
{
    out << "%FLAG RADII" << endl
        << "%FORMAT(5E16.8)" << endl;
    int count = 0;
    int index = 0;
    const int MAX_IN_LINE = 5;
    const int ITEM_LENGTH = 16;
    for(int i = 0; i < number_of_residues_; i++)
    {
        TopologyResidue* residue = this->assembly_->GetResidueByIndex(i+1);
        for(unsigned int j = 0; j < residue->GetAtoms().size(); j++)
        {
            TopologyAtom* atom = residue->GetAtomByIndex(index+1);
            index++;
            if(atom->GetRadii() != dNotSet)
                out << setw(ITEM_LENGTH) << right << scientific << setprecision(8) << atom->GetRadii();
            else
                out << setw(ITEM_LENGTH) << right << scientific << setprecision(8) << 0.0;
            count++;
            if(count == MAX_IN_LINE)
            {
                count = 0;
                out << endl;
            }
        }
    }
    if(count < MAX_IN_LINE)
        out << endl;
}

void TopologyFile::ResolveScreenSection(ofstream& out)
{
    out << "%FLAG SCREEN" << endl
        << "%FORMAT(5E16.8)" << endl;
    int count = 0;
    int index = 0;
    const int MAX_IN_LINE = 5;
    const int ITEM_LENGTH = 16;
    for(int i = 0; i < number_of_residues_; i++)
    {
        TopologyResidue* residue = this->assembly_->GetResidueByIndex(i+1);
        for(unsigned int j = 0; j < residue->GetAtoms().size(); j++)
        {
            TopologyAtom* atom = residue->GetAtomByIndex(index+1);
            index++;
            if(atom->GetScreen() != dNotSet)
                out << setw(ITEM_LENGTH) << right << scientific << setprecision(8) << atom->GetScreen();
            else
                out << setw(ITEM_LENGTH) << right << scientific << setprecision(8) << 0.0;
            count++;
            if(count == MAX_IN_LINE)
            {
                count = 0;
                out << endl;
            }
        }
    }
    if(count < MAX_IN_LINE)
        out << endl;
}

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void TopologyFile::Print(ostream &out)
{
    out << "============================ " << title_ << " ===========================" << endl;
    assembly_->Print(out);
    out << "------------------------------ Atom Pairs ------------------------------" << endl;
    for(TopologyAtomPairMap::iterator it = pairs_.begin(); it != pairs_.end(); it++)
    {
        TopologyAtomPair* atom_pair = (*it).second;
        atom_pair->Print(out);
    }
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
