#ifndef TOPOLOGYFILE_HPP
#define TOPOLOGYFILE_HPP

#include <string>
#include <iostream>
#include <map>
#include <vector>

namespace TopologyFileSpace
{
    class TopologyAtomPair;
    class TopologyBondType;
    class TopologyAngleType;
    class TopologyDihedralType;
    class TopologyAssembly;
    class TopologyBond;
    class TopologyAngle;
    class TopologyDihedral;
    class TopologyFile
    {
        public:
            //////////////////////////////////////////////////////////
            //                    TYPE DEFINITION                   //
            //////////////////////////////////////////////////////////
            typedef std::map<std::string, TopologyAtomPair*> TopologyAtomPairMap;
            typedef std::map<int, TopologyBondType*> TopologyBondTypeMap;
            typedef std::map<int, TopologyAngleType*> TopologyAngleTypeMap;
            typedef std::map<int, TopologyDihedralType*> TopologyDihedralTypeMap;
            typedef std::vector<std::string> RadiusSet;
            typedef std::map<std::string, TopologyBond*> TopologyBondMap;
            typedef std::map<std::string, TopologyAngle*> TopologyAngleMap;
            typedef std::map<std::string, TopologyDihedral*> TopologyDihedralMap;
            typedef std::map<std::string, int> TopologyAtomTypeIndexMap;

            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            TopologyFile();
            /*! \fn
              * Constructor to build the current topology file from the given file name
              * @param top_file An existing topology file path to be read
              */
            TopologyFile(const std::string& top_file);

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
            /*! \fn
              * An accessor function in order to access to the path
              * @return path_ attribute of the current object of this class
              */
            std::string GetPath();
            /*! \fn
              * An accessor function in order to access to the title
              * @return title_ attribute of the current object of this class
              */
            std::string GetTitle();
            /*! \fn
              * An accessor function in order to access to the number of atoms
              * @return number_of_atoms_ attribute of the current object of this class
              */
            int GetNumberOfAtoms();
            /*! \fn
              * An accessor function in order to access to the number of types
              * @return number_of_types_ attribute of the current object of this class
              */
            int GetNumberOfTypes();
            /*! \fn
              * An accessor function in order to access to the number of bonds including hydrogen
              * @return number_of_bonds_including_hydrogen_ attribute of the current object of this class
              */
            int GetNumberOfBondsIncludingHydrogen();
            /*! \fn
              * An accessor function in order to access to the number of angles including hydrogen
              * @return number_of_angles_including_hydrogen_ attribute of the current object of this class
              */
            int GetNumberOfAnglesIncludingHydrogen();
            /*! \fn
              * An accessor function in order to access to the number of dihedrals including hydrogen
              * @return number_of_dihedrals_including_hydrogen_ attribute of the current object of this class
              */
            int GetNumberOfDihedralsIncludingHydrogen();
            /*! \fn
              * An accessor function in order to access to the number of bonds excluding hydrogen
              * @return number_of_bonds_excluding_hydrogen_ attribute of the current object of this class
              */
            int GetNumberOfBondsExcludingHydrogen();
            /*! \fn
              * An accessor function in order to access to the number of angles excluding hydrogen
              * @return number_of_angles_excluding_hydrogen_ attribute of the current object of this class
              */
            int GetNumberOfAnglesExcludingHydrogen();
            /*! \fn
              * An accessor function in order to access to the number of dihedrals excluding hydrogen
              * @return number_of_dihedrals_excluding_hydrogen_ attribute of the current object of this class
              */
            int GetNumberOfDihedralsExcludingHydrogen();
            /*! \fn
              * An accessor function in order to access to the number of hydrogen parameters
              * @return number_of_hydrogen_parameters_ attribute of the current object of this class
              */
            int GetNumberOfHydrogenParameters();
            /*! \fn
              * An accessor function in order to access to the number of parameters
              * @return number_of_parameters_ attribute of the current object of this class
              */
            int GetNumberOfParameters();
            /*! \fn
              * An accessor function in order to access to the number of excluded atoms
              * @return number_of_excluded_atoms_ attribute of the current object of this class
              */
            int GetNumberOfExcludedAtoms();
            /*! \fn
              * An accessor function in order to access to the number of residues
              * @return number_of_residues_ attribute of the current object of this class
              */
            int GetNumberOfResidues();
            /*! \fn
              * An accessor function in order to access to the total number of bonds
              * @return total_number_of_bonds_ attribute of the current object of this class
              */
            int GetTotalNumberOfBonds();
            /*! \fn
              * An accessor function in order to access to the total number of angles
              * @return total_number_of_angles_ attribute of the current object of this class
              */
            int GetTotalNumberOfAngles();
            /*! \fn
              * An accessor function in order to access to the total number of dihedrals
              * @return total_number_of_dihedrals_ attribute of the current object of this class
              */
            int GetTotalNumberOfDihedrals();
            /*! \fn
              * An accessor function in order to access to the number of bond types
              * @return number_of_bond_types_ attribute of the current object of this class
              */
            int GetNumberOfBondTypes();
            /*! \fn
              * An accessor function in order to access to the number of angle types
              * @return number_of_angles_types_ attribute of the current object of this class
              */
            int GetNumberOfAnglesTypes();
            /*! \fn
              * An accessor function in order to access to the number of dihedral types
              * @return number_of_dihedral_types_ attribute of the current object of this class
              */
            int GetNumberOfDihedralTypes();
            /*! \fn
              * An accessor function in order to access to the number of atom types in parameter file
              * @return number_of_atom_types_in_parameter_file_ attribute of the current object of this class
              */
            int GetNumberOfAtomTypesInParameterFile();
            /*! \fn
              * An accessor function in order to access to the number of distinct hydrogen bonds
              * @return number_of_distinct_hydrogen_bonds_ attribute of the current object of this class
              */
            int GetNumberOfDistinctHydrogenBonds();
            /*! \fn
              * An accessor function in order to access to the boolean attribute that indicates if topology file has perturbation
              * @return perturbation_option attribute of the current object of this class
              */
            int GetPerturbationOption();
            /*! \fn
              * An accessor function in order to access to the number of bonds perturbed
              * @return number_of_bonds_perturbed_ attribute of the current object of this class
              */
            int GetNumberOfBondsPerturbed();
            /*! \fn
              * An accessor function in order to access to the number of angles perturbed
              * @return number_of_angles_perturbed_ attribute of the current object of this class
              */
            int GetNumberOfAnglesPerturbed();
            /*! \fn
              * An accessor function in order to access to the number of dihedrals perturbed
              * @return number_of_dihedrals_perturbed_ attribute of the current object of this class
              */
            int GetNumberOfDihedralsPerturbed();
            /*! \fn
              * An accessor function in order to access to the number of bonds group perturbed
              * @return number_of_bonds_group_perturbed_ attribute of the current object of this class
              */
            int GetNumberOfBondsGroupPerturbed();
            /*! \fn
              * An accessor function in order to access to the number of angles group perturbed
              * @return number_of_angles_group_perturbed_ attribute of the current object of this class
              */
            int GetNumberOfAnglesGroupPerturbed();
            /*! \fn
              * An accessor function in order to access to the number of dihedrals group perturbed
              * @return number_of_dihedrals_group_perturbed_ attribute of the current object of this class
              */
            int GetNumberOfDihedralsGroupPerturbed();
            /*! \fn
              * An accessor function in order to access to the boolean attribute that indicates if topology file has standard periodic box
              * @return standard_periodic_box_option attribute of the current object of this class
              */
            int GetStandardPeriodicBoxOption();
            /*! \fn
              * An accessor function in order to access to the number of atoms in largest residue
              * @return number_of_atoms_in_largest_residue_ attribute of the current object of this class
              */
            int GetNumberOfAtomsInLargestResidue();
            /*! \fn
              * An accessor function in order to access to the boolean attribute that indicates if topology file has CAP option
              * @return cap_option_ attribute of the current object of this class
              */
            int GetCapOption();
            /*! \fn
              * An accessor function in order to access to the number of extra points
              * @return number_of_extra_points_ attribute of the current object of this class
              */
            int GetNumberOfExtraPoints();
            /*! \fn
              * An accessor function in order to access to the number of beads
              * @return number_of_beads_ attribute of the current object of this class
              */
            int GetNumberOfBeads();
            /*! \fn
              * An accessor function in order to access to the pairs
              * @return pairs_ attribute of the current object of this class
              */
            TopologyAtomPairMap GetPairs();
            /*! \fn
              * An accessor function in order to access to the bond types
              * @return bond_types_ attribute of the current object of this class
              */
            TopologyBondTypeMap GetBondTypes();
            /*! \fn
              * An accessor function in order to access to the angle types
              * @return angle_types_ attribute of the current object of this class
              */
            TopologyAngleTypeMap GetAngleTypes();
            /*! \fn
              * An accessor function in order to access to the atom types
              * @return atom_types_ attribute of the current object of this class
              */
            TopologyDihedralTypeMap GetDihedralTypes();
            /*! \fn
              * An accessor function in order to access to the assembly
              * @return assembly_ attribute of the current object of this class
              */
            TopologyAssembly* GetAssembly();

            RadiusSet GetRadiusSet();
            /*! \fn
              * An accessor function in order to access to the bonds
              * @return bonds_ attribute of the current object of this class
              */
            TopologyBondMap GetBonds();
            /*! \fn
              * An accessor function in order to access to the angles
              * @return angles_ attribute of the current object of this class
              */
            TopologyAngleMap GetAngles();
            /*! \fn
              * An accessor function in order to access to the dihedrals
              * @return dihedrals_ attribute of the current object of this class
              */
            TopologyDihedralMap GetDihedrals();

            /*! \fn
              * A function to return a map with atom type key and corresponding index
              * @return Mapping of all atom types in the topology file and their corresponding index
              */
            TopologyAtomTypeIndexMap GetAtomTypeIndexMap();
            /*! \fn
              * A function to retreive a bond type correspondence to the given index
              * @param index Target bond index
              * @return A pointer to a topology bond type object corresponding to the given index
              */
            TopologyBondType* GetBondTypeByIndex(int index);
            /*! \fn
              * A function to retreive an angle type correspondence to the given index
              * @param index Target angle index
              * @return A pointer to a topology angle type object corresponding to the given index
              */
            TopologyAngleType* GetAngleTypeByIndex(int index);
            /*! \fn
              * A function to retreive a dihedral type correspondence to the given index
              * @param index Target dihedral index
              * @return A pointer to a topology dihedral type corresponding to the given index
              */
            TopologyDihedralType* GetDihedralTypeByIndex(int index);
            /*! \fn
              * A function to retreive a pair of atoms correspondence to the given index
              * @param index Target topology pair
              * @return A Pointer to a topology pair corresponding to the given index
              */
            TopologyAtomPair* GetAtomPairByIndex(int index);

            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A mutator function in order to set the path of the current object
              * Set the path_ attribute of the current topology file
              * @param path The path attribute of the current object
              */
            void SetPath(std::string path);
            /*! \fn
              * A mutator function in order to set the title of the current object
              * Set the title_ attribute of the current topology file
              * @param title The title attribute of the current object
              */
            void SetTitle(std::string title);
            /*! \fn
              * A mutator function in order to set the number of atoms of the current object
              * Set the number_of_atoms_ attribute of the current topology file
              * @param number_of_atoms The number of atoms attribute of the current object
              */
            void SetNumberOfAtoms(int number_of_atoms);
            /*! \fn
              * A mutator function in order to set the number_of_types of the current object
              * Set the number_of_types_ attribute of the current topology file
              * @param number_of_types The number of types attribute of the current object
              */
            void SetNumberOfTypes(int number_of_types);
            /*! \fn
              * A mutator function in order to set the number of bonds including hydrogen of the current object
              * Set the number_of_bonds_including_hydrogen_ attribute of the current topology file
              * @param number_of_bonds_including_hydrogen The numberof bonds including hydrogen attribute of the current object
              */
            void SetNumberOfBondsIncludingHydrogen(int number_of_bonds_including_hydrogen);
            /*! \fn
              * A mutator function in order to set the number_of_angles including hydrogen of the current object
              * Set the number_of_angles_including_hydrogen_ attribute of the current topology file
              * @param number_of_angles_including_hydrogen The number of angles including hydrogen attribute of the current object
              */
            void SetNumberOfAnglesIncludingHydrogen(int number_of_angles_including_hydrogen);
            /*! \fn
              * A mutator function in order to set the number of dihedrals including hydrogen of the current object
              * Set the number_of_dihedrals_including_hydrogen_ attribute of the current topology file
              * @param number_of_dihedrals_including_hydrogen The number of dihedrals including hydrogen attribute of the current object
              */
            void SetNumberOfDihedralsIncludingHydrogen(int number_of_dihedrals_including_hydrogen);
            /*! \fn
              * A mutator function in order to set the number of bonds excluding hydrogen of the current object
              * Set the number_of_bonds_excluding_hydrogen_ attribute of the current topology file
              * @param number_of_bonds_excluding_hydrogen The number of bonds excluding hydrogen attribute of the current object
              */
            void SetNumberOfBondsExcludingHydrogen(int number_of_bonds_excluding_hydrogen);
            /*! \fn
              * A mutator function in order to set the number of angles excluding hydrogen of the current object
              * Set the number_of_angles_excluding_hydrogen_ attribute of the current topology file
              * @param number_of_angles_excluding_hydrogen The number of angles excluding_hydrogen attribute of the current object
              */
            void SetNumberOfAnglesExcludingHydrogen(int number_of_angles_excluding_hydrogen);
            /*! \fn
              * A mutator function in order to set the number of dihedrals excluding hydrogen of the current object
              * Set the number_of_dihedrals_excluding_hydrogen_ attribute of the current topology file
              * @param number_of_dihedrals_excluding_hydrogen The number of dihedrals excluding hydrogen attribute of the current object
              */
            void SetNumberOfDihedralsExcludingHydrogen(int number_of_dihedrals_excluding_hydrogen);
            /*! \fn
              * A mutator function in order to set the number of hydrogen parameters of the current object
              * Set the number_of_hydrogen_parameters_ attribute of the current topology file
              * @param number_of_hydrogen_parameters The number of hydrogen parameters attribute of the current object
              */
            void SetNumberOfHydrogenParameters(int number_of_hydrogen_parameters);
            /*! \fn
              * A mutator function in order to set the number of parameters of the current object
              * Set the number_of_parameters_ attribute of the current topology file
              * @param number_of_parameters The number of parameters attribute of the current object
              */
            void SetNumberOfParameters(int number_of_parameters);
            /*! \fn
              * A mutator function in order to set the number of excluded atoms of the current object
              * Set the number_of_excluded_atoms_ attribute of the current topology file
              * @param number_of_excluded_atoms The number of excluded atoms attribute of the current object
              */
            void SetNumberOfExcludedAtoms(int number_of_excluded_atoms);
            /*! \fn
              * A mutator function in order to set the number of residues of the current object
              * Set the number_of_residues_ attribute of the current topology file
              * @param number_of_residues The number of residues attribute of the current object
              */
            void SetNumberOfResidues(int number_of_residues);
            /*! \fn
              * A mutator function in order to set the total number of bonds of the current object
              * Set the total_number_of_bonds_ attribute of the current topology file
              * @param total_number_of_bonds The total number of bonds attribute of the current object
              */
            void SetTotalNumberOfBonds(int total_number_of_bonds);
            /*! \fn
              * A mutator function in order to set the total number of angles of the current object
              * Set the total_number_of_angles_ attribute of the current topology file
              * @param total_number_of_angles The total number of angles attribute of the current object
              */
            void SetTotalNumberOfAngles(int total_number_of_angles);
            /*! \fn
              * A mutator function in order to set the total number of dihedrals of the current object
              * Set the total_number_of_dihedrals_ attribute of the current topology file
              * @param total_number_of_dihedrals The total number of dihedrals attribute of the current object
              */
            void SetTotalNumberOfDihedrals(int total_number_of_dihedrals);
            /*! \fn
              * A mutator function in order to set the number_of_bond_types of the current object
              * Set the number_of_bond_types_ attribute of the current topology file
              * @param number_of_bond_types The number of bond types attribute of the current object
              */
            void SetNumberOfBondTypes(int number_of_bond_types);
            /*! \fn
              * A mutator function in order to set the number of angle types of the current object
              * Set the number_of_angle_types_ attribute of the current topology file
              * @param number_of_angle_types The number of angle types attribute of the current object
              */
            void SetNumberOfAngleTypes(int number_of_angle_types);
            /*! \fn
              * A mutator function in order to set the number of dihedrals of the current object
              * Set the number_of_dihedrals_ attribute of the current topology file
              * @param number_of_dihedrals The number of dihedrals attribute of the current object
              */
            void SetNumberOfDihedralTypes(int number_of_dihedral_types);
            /*! \fn
              * A mutator function in order to set the number of atom types in parameter file of the current object
              * Set the number_of_atom_types_in_parameter_file_ attribute of the current topology file
              * @param number_of_atom_types_in_parameter_file The number of atom types in parameter file attribute of the current object
              */
            void SetNumberOfAtomTypesInParameterFile(int number_of_atom_types_in_parameter_file);
            /*! \fn
              * A mutator function in order to set the number of distinct hydrogen bonds of the current object
              * Set the number_of_distinct_hydrogen_bonds_ attribute of the current topology file
              * @param number_of_distinct_hydrogen_bonds The number of distinct hydrogen bonds attribute of the current object
              */
            void SetNumberOfDistinctHydrogenBonds(int number_of_distinct_hydrogen_bonds);
            /*! \fn
              * A mutator function in order to set the has perturbation attribute of the current object
              * Set the perturbation_option_ attribute of the current topology file
              * @param perturbation_option The has perturbation attribute of the current object
              */
            void SetPerturbationOption(int perturbation_option);
            /*! \fn
              * A mutator function in order to set the number of bonds perturbed of the current object
              * Set the number_of_bonds_perturbed_ attribute of the current topology file
              * @param number_of_bonds_perturbed The number of bonds perturbed attribute of the current object
              */
            void SetNumberOfBondsPerturbed(int number_of_bonds_perturbed);
            /*! \fn
              * A mutator function in order to set the number of angles perturbed of the current object
              * Set the number_of_angles_perturbed_ attribute of the current topology file
              * @param number_of_angles_perturbed The number of angles perturbed attribute of the current object
              */
            void SetNumberOfAnglesPerturbed(int number_of_angles_perturbed);
            /*! \fn
              * A mutator function in order to set the number of dihedrals perturbed of the current object
              * Set the number_of_dihedrals_perturbed_ attribute of the current topology file
              * @param number_of_dihedrals_perturbed The number of dihedrals perturbed attribute of the current object
              */
            void SetNumberOfDihedralsPerturbed(int number_of_dihedrals_perturbed);
            /*! \fn
              * A mutator function in order to set the number of bonds group perturbed of the current object
              * Set the number_of_bonds_group_perturbed_ attribute of the current topology file
              * @param number_of_bonds_group_perturbed The number of bonds group perturbed attribute of the current object
              */
            void SetNumberOfBondsGroupPerturbed(int number_of_bonds_group_perturbed);
            /*! \fn
              * A mutator function in order to set the number of angles group perturbed of the current object
              * Set the number_of_angles_group_perturbed_ attribute of the current topology file
              * @param number_of_angles_group_perturbed The number of angles group perturbed attribute of the current object
              */
            void SetNumberOfAnglesGroupPerturbed(int number_of_angles_group_perturbed);
            /*! \fn
              * A mutator function in order to set the number of dihedrals group perturbed of the current object
              * Set the number_of_dihedrals_group_perturbed_ attribute of the current topology file
              * @param number_of_dihedrals_group_perturbed The number of dihedrals group perturbed attribute of the current object
              */
            void SetNumberOfDihedralsGroupPerturbed(int number_of_dihedrals_group_perturbed);
            /*! \fn
              * A mutator function in order to set the has standard periodic box attribute of the current object
              * Set the standard_periodic_box_option_ attribute of the current topology file
              * @param standard_periodic_box_option The has standard periodic box attribute of the current object
              */
            void SetStandardPeriodicBoxOption(int standard_periodic_box_option);
            /*! \fn
              * A mutator function in order to set the number of atoms in largest residue attribute of the current object
              * Set the number_of_atoms_in_largest_residue_ attribute of the current topology file
              * @param number_of_atoms_in_largest_residue The number of atoms in largest residue attribute of the current object
              */
            void SetNumberOfAtomsInLargestResidue(int number_of_atoms_in_largest_residue);
            /*! \fn
              * A mutator function in order to set the has cap option attribute of the current object
              * Set the cap_option_ attribute of the current topology file
              * @param cap_option The has cap option attribute of the current object
              */
            void SetCapOption(int cap_option);
            /*! \fn
              * A mutator function in order to set the number of extra points of the current object
              * Set the number_of_extra_points_ attribute of the current topology file
              * @param number_of_extra_points The number_of_extra_points attribute of the current object
              */
            void SetNumberOfExtraPoints(int number_of_extra_points);
            /*! \fn
              * A mutator function in order to set the number of beads of the current object
              * Set the number_of_beads_ attribute of the current topology file
              * @param number_of_beads The number of beads attribute of the current object
              */
            void SetNumberOfBeads(int number_of_beads);
            /*! \fn
              * A mutator function in order to set the assembly for the current object
              * Set the assembly_ attribute of the current topology file
              * @param assembly The number of beads attribute of the current object
              */
            void SetAssembly(TopologyAssembly* assembly);
            /*! \fn
              * A mutator function in order to set the residue set for the current object
              * Set the residue_set_ attribute of the current topology file
              * @param residue_set Residue set
              */
            void SetRadiusSet(RadiusSet radius_set);
            /*! \fn
              * A mutator function in order to set the bonds of the current object
              * Set the bonds_ attribute of the current topology file
              * @param bonds The bonds attribute of the current object
              */
            void SetBonds(TopologyBondMap bonds);
            /*! \fn
              * A function in order to add the bond to the current object
              * Set the bonds_ attribute of the current topology file
              * @param bond The bond that will be added to the current object
              */
            void AddBond(TopologyBond* bond);
            /*! \fn
              * A mutator function in order to set the bond types of the current object
              * Set the bond_types_ attribute of the current topology file
              * @param bond_types The bond_types attribute of the current object
              */
            void SetBondTypes(TopologyBondTypeMap bond_types);
            /*! \fn
              * A function in order to add the bond type to the current object
              * Set the bond_Types_ attribute of the current topology file
              * @param bond_type The bond_type that will be added to the current object
              */
            void AddBondType(TopologyBondType* bond_type);
            /*! \fn
              * A mutator function in order to set the angles of the current object
              * Set the angles_ attribute of the current topology file
              * @param angles The angles attribute of the current object
              */
            void SetAngles(TopologyAngleMap angles);
            /*! \fn
              * A function in order to add the angle to the current object
              * Set the angles_ attribute of the current topology file
              * @param angle The angle that will be added to the current object
              */
            void AddAngle(TopologyAngle* angle);
            /*! \fn
              * A mutator function in order to set the angle types of the current object
              * Set the angle_types_ attribute of the current topology file
              * @param angle_types The angle_types attribute of the current object
              */
            void SetAngleTypes(TopologyAngleTypeMap angle_types);
            /*! \fn
              * A function in order to add the angle type to the current object
              * Set the angle_Types_ attribute of the current topology file
              * @param angle_type The angle_type that will be added to the current object
              */
            void AddAngleType(TopologyAngleType* angle_type);
            /*! \fn
              * A mutator function in order to set the dihedrals of the current object
              * Set the dihedrals_ attribute of the current topology file
              * @param dihedrals The dihedrals attribute of the current object
              */
            void SetDihedrals(TopologyDihedralMap dihedrals);
            /*! \fn
              * A function in order to add the dihedral to the current object
              * Set the dihedrals_ attribute of the current topology file
              * @param dihedral The dihedral that will be added to the current object
              */
            void AddDihedral(TopologyDihedral* dihedral);
            /*! \fn
              * A mutator function in order to set the dihedral types of the current object
              * Set the dihedral_types_ attribute of the current topology file
              * @param dihedral_types The dihedral_types attribute of the current object
              */
            void SetDihedralTypes(TopologyDihedralTypeMap dihedral_types);
            /*! \fn
              * A function in order to add the dihedral type to the current object
              * Set the dihedral_Types_ attribute of the current topology file
              * @param dihedral_type The dihedral_type that will be added to the current object
              */
            void AddDihedralType(TopologyDihedralType* dihedral_type);
            /*! \fn
              * A function in order to set the atom pairs of the current object
              * Set the pairs_ attribute of the current topology file
              * @param pairs The pairs attribute of the current object
              */
            void SetAtomPairs(TopologyAtomPairMap pairs);
            /*! \fn
              * A function in order to add the atom pairs to the current object
              * Set the pairs_ attribute of the current topology file
              * @param pair The atom pair that will be added to the current object
              */
            void AddAtomPair(TopologyAtomPair* pair);

            //////////////////////////////////////////////////////////
            //                        FUNCTIONS                     //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to parse the contents of a given stream of a file
              * Parse the given stream and set the attributes of the current object accordingly
              * @param in_file A stream contains whole contents of a topology file
              */
            void Read(std::ifstream& in_file);
            /*! \fn
              * A function to parse the contents of a given stream of a file
              * @param in_stream A stream contains whole contents of a topology file
              */
            void ParseSections(std::ifstream& in_stream);
            /*! \fn
              * A function to partition a topology file to its sections
              * @param stream Input stream built from topology input file
              * @param line The line that the reader is currently at
              * @param section Stream that the output section will be written on
              */
            void PartitionSection(std::ifstream& stream, std::string& line, std::stringstream& section);
            /*! \fn
              * A function to parse title section of a pdb file
              * @param stream Title section of a topology file as stream
              */
            void ParseTitlePartition(std::stringstream& stream);
            /*! \fn
              * A function to parse the pointer section of a topology file
              * @param stream Pointer section of a topology file as stream
              */
            void ParsePointersPartition(std::stringstream& stream);
            /*! \fn
              * A function to parse the other partitions of a topology file which have same template
              * @param stream Stream of a section that has to be parsed
              * @return Vector of types according to the section
              */
            template<typename T>
            std::vector<T> ParsePartition(std::stringstream& stream);
            /*! \fn
              * A function to parse a line of a section to its values
              * @param line Line of a section
              * @param format Format of the line which indicates the number and type of the values
              * @return A vector of values in the given line parsed based on the given format
              */
            template<typename T>
            std::vector<T> PartitionLine(std::string line, std::string format);
            /*! \fn
              * A function to create an output topology file with the given name
              * @param top_file Output topology file name
              */
            void Write(const std::string& top_file);
            /*! \fn
              * A function that resolve sections of a topology file in order to write into an output stream
              * @param out_stream Stream the the result has to be written into
              */
            void ResolveSections(std::ofstream& out_stream);
            /*! \fn
              * A function to write back title section of the topology file into an output stream
              * @param out Intermediate output stream in order to write title section
              */
            void ResolveTitleSection(std::ofstream& out);
            /*! \fn
              * A function to write back pointer section of the topology file into an output stream
              * @param out Intermediate output stream in order to write pointer section
              */
            void ResolveGitVersion(std::ofstream& out);
            /*! \fn
              * A function to write back version section of the topology file into an output stream
              * @param out Intermediate output stream in order to write pointer section
              */
            void ResolvePointersSection(std::ofstream& out);
            /*! \fn
              * A function to write back atom name section of the topology file into an output stream
              * @param out Intermediate output stream in order to write atom name section
              */
            void ResolveAtomNameSection(std::ofstream& out);
            /*! \fn
              * A function to write back charge section of the topology file into an output stream
              * @param out Intermediate output stream in order to write charge section
              */
            void ResolveChargeSection(std::ofstream& out);
            /*! \fn
              * A function to write back atomic number section of the topology file into an output stream
              * @param out Intermediate output stream in order to write atomic number section
              */
            void ResolveAtomicNumberSection(std::ofstream& out);
            /*! \fn
              * A function to write back mass section of the topology file into an output stream
              * @param out Intermediate output stream in order to write mass section
              */
            void ResolveMassSection(std::ofstream& out);
            /*! \fn
              * A function to write back atom type index section of the topology file into an output stream
              * @param out Intermediate output stream in order to write atom type index section
              */
            void ResolveAtomTypeIndexSection(std::ofstream& out);
            /*! \fn
              * A function to write back number excluded atoms section of the topology file into an output stream
              * @param out Intermediate output stream in order to write number excluded atoms section
              */
            void ResolveNumberExcludedAtomsSection(std::ofstream& out);
            /*! \fn
              * A function to write back non bonded parm index section of the topology file into an output stream
              * @param out Intermediate output stream in order to write non bonded parm index section
              */
            void ResolveNonbondedParmIndexSection(std::ofstream& out);
            /*! \fn
              * A function to write back residue lab section of the topology file into an output stream
              * @param out Intermediate output stream in order to write residue lab section
              */
            void ResolveResidueLabelSection(std::ofstream& out);
            /*! \fn
              * A function to write back residue pointers section of the topology file into an output stream
              * @param out Intermediate output stream in order to write residue pointers section
              */
            void ResolveResiduePointersSection(std::ofstream& out);
            /*! \fn
              * A function to write back bond force constant section of the topology file into an output stream
              * @param out Intermediate output stream in order to write bond force constant section
              */
            void ResolveBondForceConstantSection(std::ofstream& out);
            /*! \fn
              * A function to write back bond equil value section of the topology file into an output stream
              * @param out Intermediate output stream in order to write bond equil value section
              */
            void ResolveBondEquilValueSection(std::ofstream& out);
            /*! \fn
              * A function to write back angle force constant section of the topology file into an output stream
              * @param out Intermediate output stream in order to write angle force constant section
              */
            void ResolveAngleForceConstantSection(std::ofstream& out);
            /*! \fn
              * A function to write back angle equil value section of the topology file into an output stream
              * @param out Intermediate output stream in order to write angle equil value section
              */
            void ResolveAngleEquilValueSection(std::ofstream& out);
            /*! \fn
              * A function to write back dihedral force constant section of the topology file into an output stream
              * @param out Intermediate output stream in order to write dihedral force constant section
              */
            void ResolveDihedralForceConstantSection(std::ofstream& out);
            /*! \fn
              * A function to write back dihedral periodicity section of the topology file into an output stream
              * @param out Intermediate output stream in order to write dihedral periodicity section
              */
            void ResolveDihedralPeriodicitySection(std::ofstream& out);
            /*! \fn
              * A function to write back dihedral phase section of the topology file into an output stream
              * @param out Intermediate output stream in order to write dihedral phase section
              */
            void ResolveDihedralPhaseSection(std::ofstream& out);
            /*! \fn
              * A function to write back scee scale factor section of the topology file into an output stream
              * @param out Intermediate output stream in order to write scee scale factor section
              */
            void ResolveSceeScaleFactorSection(std::ofstream& out);
            /*! \fn
              * A function to write back scnb scale factor section of the topology file into an output stream
              * @param out Intermediate output stream in order to write scnb scale factor section
              */
            void ResolveScnbScaleFactorSection(std::ofstream& out);
            /*! \fn
              * A function to write back solty section of the topology file into an output stream
              * @param out Intermediate output stream in order to write solty section
              */
            void ResolveSoltySection(std::ofstream& out);
            /*! \fn
              * A function to write back lennard jones a coefficient section of the topology file into an output stream
              * @param out Intermediate output stream in order to write lennard jones a coefficient section
              */
            void ResolveLennardJonesACoefSection(std::ofstream& out);
            /*! \fn
              * A function to write back lennard jones b coefficient section of the topology file into an output stream
              * @param out Intermediate output stream in order to write lennard jones b coefficient section
              */
            void ResolveLennardJonesBCoefSection(std::ofstream& out);
            /*! \fn
              * A function to write back bonds inc hydrogen section of the topology file into an output stream
              * @param out Intermediate output stream in order to write bonds inc hydrogen section
              */
            void ResolveBondsIncHydrogenSection(std::ofstream& out);
            /*! \fn
              * A function to write back bonds without hydrogen section of the topology file into an output stream
              * @param out Intermediate output stream in order to write bonds without hydrogen section
              */
            void ResolveBondsWithoutHydrogenSection(std::ofstream& out);
            /*! \fn
              * A function to write back angles inc hydrogen section of the topology file into an output stream
              * @param out Intermediate output stream in order to write angles including hydrogen section
              */
            void ResolveAnglesIncHydrogenSection(std::ofstream& out);
            /*! \fn
              * A function to write back angles without hydrogen section of the topology file into an output stream
              * @param out Intermediate output stream in order to write angles without hydrogen section
              */
            void ResolveAnglesWithoutHydrogenSection(std::ofstream& out);
            /*! \fn
              * A function to write back dihedrals inc hydrogen section of the topology file into an output stream
              * @param out Intermediate output stream in order to write dihedrals inc hydrogen section
              */
            void ResolveDihedralsIncHydrogenSection(std::ofstream& out);
            /*! \fn
              * A function to write back dihedrals without hydrogen section of the topology file into an output stream
              * @param out Intermediate output stream in order to write dihedral without hydrogen section
              */
            void ResolveDihedralsWithoutHydrogenSection(std::ofstream& out);
            /*! \fn
              * A function to write back excluded atoms list section of the topology file into an output stream
              * @param out Intermediate output stream in order to write excluded atoms list section
              */
            void ResolveExcludedAtomsListSection(std::ofstream& out);
            /*! \fn
              * A function to write back hydrogen bond a coefficient section of the topology file into an output stream
              * @param out Intermediate output stream in order to write hydrogen bond a coefficient section
              */
            void ResolveHydrogenBondACoefSection(std::ofstream& out);
            /*! \fn
              * A function to write back hydrogen bond b coefficient section of the topology file into an output stream
              * @param out Intermediate output stream in order to write hydrogen bond coefficient section
              */
            void ResolveHydrogenBondBCoefSection(std::ofstream& out);
            /*! \fn
              * A function to write back hb cut section of the topology file into an output stream
              * @param out Intermediate output stream in order to write hb cut section
              */
            void ResolveHBCutSection(std::ofstream& out);
            /*! \fn
              * A function to write back amber atom type section of the topology file into an output stream
              * @param out Intermediate output stream in order to write amber atom type section
              */
            void ResolveAmberAtomTypeSection(std::ofstream& out);
            /*! \fn
              * A function to write back tree chain classification section of the topology file into an output stream
              * @param out Intermediate output stream in order to write tree chain classification section
              */
            void ResolveTreeChainClassificationSection(std::ofstream& out);
            /*! \fn
              * A function to write back join array section of the topology file into an output stream
              * @param out Intermediate output stream in order to write join array section
              */
            void ResolveJoinArraySection(std::ofstream& out);
            /*! \fn
              * A function to write back irotat section of the topology file into an output stream
              * @param out Intermediate output stream in order to write title section
              */
            void ResolveIRotatSection(std::ofstream& out);
            /*! \fn
              * A function to write back residue set section of the topology file into an output stream
              * @param out Intermediate output stream in order to write residue set section
              */
            void ResolveRadiusSetSection(std::ofstream& out);
            /*! \fn
              * A function to write back radii section of the topology file into an output stream
              * @param out Intermediate output stream in order to write radii section
              */
            void ResolveRadiiSection(std::ofstream& out);
            /*! \fn
              * A function to write solvent pointer section of the topology file into an output stream
              * @param out Intermediate output stream in order to write screen section
              */
           void ResolveSolventPointersSection(std::ofstream& out);
            /*! \fn
              * A function to write back screen section of the topology file into an output stream
              * @param out Intermediate output stream in order to write screen section
              */
            void ResolveScreenSection(std::ofstream& out);

            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the topology file contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cerr);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            std::string path_;                              /*!< Path of the file >*/
            std::string title_;                             /*!< Title of the file >*/
            int number_of_atoms_;                           /*!< Number of atoms >*/
            int number_of_types_;                           /*!< Number of types in a topology file >*/
            int number_of_bonds_including_hydrogen_;        /*!< Number of bonds including hydrogen in a topology file >*/
            int number_of_angles_including_hydrogen_;       /*!< Number of angles including hydrogen in a topology file >*/
            int number_of_dihedrals_including_hydrogen_;    /*!< Number of dihedrals including hydrogen in a topology file >*/
            int number_of_bonds_excluding_hydrogen_;        /*!< Number of bonds excluding hydrogen in a topology file >*/
            int number_of_angles_excluding_hydrogen_;       /*!< Number of angles excluding hydrogen in a topology file >*/
            int number_of_dihedrals_excluding_hydrogen_;    /*!< Number of dihedrals excluding hydrogen in a topology file >*/
            int number_of_hydrogen_parameters_;             /*!< Number of hydrogen parameters in a topology file >*/
            int number_of_parameters_;                      /*!< Number of parameters in a topology file >*/
            int number_of_excluded_atoms_;                  /*!< Number of excluded atoms in a topology file >*/
            int number_of_residues_;                        /*!< Number of residues in a topology file >*/
            int total_number_of_bonds_;                     /*!< Total number of bonds in a topology file >*/
            int total_number_of_angles_;                    /*!< Total number of angles in a topology file >*/
            int total_number_of_dihedrals_;                 /*!< Total number of dihedrals in a topology file >*/
            int number_of_bond_types_;                      /*!< Number of bond types in a topology file >*/
            int number_of_angle_types_;                     /*!< Number of angle types in a topology file >*/
            int number_of_dihedral_types_;                  /*!< Number of dihedral types in a topology file >*/
            int number_of_atom_types_in_parameter_file_;    /*!< Number of atom types in parameter file in a topology file >*/
            int number_of_distinct_hydrogen_bonds_;         /*!< Number of distinct hydrogen bonds in a topology file >*/
            int perturbation_option_;                       /*!< Perturbation option in a topology file >*/
            int number_of_bonds_perturbed_;                 /*!< Number of bonds perturbed in a topology file >*/
            int number_of_angles_perturbed_;                /*!< Number of angles perturbed in a topology file >*/
            int number_of_dihedrals_perturbed_;             /*!< Number of dihedrals perturbed in a topology file >*/
            int number_of_bonds_group_perturbed_;           /*!< Number of bonds group perturbed in a topology file >*/
            int number_of_angles_group_perturbed_;          /*!< Number of angles group perturbed in a topology file >*/
            int number_of_dihedrals_group_perturbed_;       /*!< Number of dihedrals group perturbed in a topology file >*/
            int standard_periodic_box_option_;              /*!< Standard periodic box option in a topology file >*/
            int number_of_atoms_in_largest_residue_;        /*!< Number of atoms in largest residue in a topology file >*/
            int cap_option_;                                /*!< Cap option in a topology file >*/
            int number_of_extra_points_;                    /*!< Number of extra points in a topology file >*/
            int number_of_beads_;                           /*!< Number of beads in a topology file >*/
            TopologyAtomPairMap pairs_;
            TopologyBondTypeMap bond_types_;
            TopologyAngleTypeMap angle_types_;
            TopologyDihedralTypeMap dihedral_types_;
            TopologyAssembly* assembly_;
            RadiusSet radius_set_;
            TopologyBondMap bonds_;
            TopologyAngleMap angles_;
            TopologyDihedralMap dihedrals_;
    };
}

#endif // TOPOLOGYFILE_HPP
