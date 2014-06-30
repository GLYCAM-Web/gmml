#ifndef TOPOLOGYFILE_HPP
#define TOPOLOGYFILE_HPP

#include <string>
#include <iostream>
#include <map>
#include <vector>

namespace TopologyFileSpace
{
    class TopologyAtomType;
    class TopologyBondType;
    class TopologyAngleType;
    class TopologyDihedralType;
    class TopologyAssembly;

    class TopologyFile
    {
        public:
            //////////////////////////////////////////////////////////
            //                    TYPE DEFINITION                   //
            //////////////////////////////////////////////////////////
            typedef std::map<int, TopologyAtomType*> TopologyAtomTypeMap;
            typedef std::map<int, TopologyBondType*> TopologyBondTypeMap;
            typedef std::map<int, TopologyAngleType*> TopologyAngleTypeMap;
            typedef std::map<int, TopologyDihedralType*> TopologyDihedralTypeMap;
            typedef std::vector<std::string> RadiusSet;

            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            TopologyFile();

            TopologyFile(const std::string& top_file);

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
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
              * An accessor function in order to access to the atom types
              * @return atom_types_ attribute of the current object of this class
              */
            TopologyAtomTypeMap GetAtomTypes();
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

            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
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

            void SetRadiusSet(RadiusSet radius_set);

            //////////////////////////////////////////////////////////
            //                        FUNCTIONS                     //
            //////////////////////////////////////////////////////////
            void Read(std::ifstream& in_file);
            void ParseSections(std::ifstream& in_stream);
            void PartitionSection(std::ifstream& stream, std::string& line, std::stringstream& section);

            void ParseTitlePartition(std::stringstream& stream);
            void ParsePointersPartition(std::stringstream& stream);
            template<typename T>
            std::vector<T> ParsePartition(std::stringstream& stream);

            template<typename T>
            std::vector<T> PartitionLine(std::string line, std::string format);

            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the topology file contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cout);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            std::string title_;
            int number_of_atoms_;
            int number_of_types_;
            int number_of_bonds_including_hydrogen_;
            int number_of_angles_including_hydrogen_;
            int number_of_dihedrals_including_hydrogen_;
            int number_of_bonds_excluding_hydrogen_;
            int number_of_angles_excluding_hydrogen_;
            int number_of_dihedrals_excluding_hydrogen_;
            int number_of_hydrogen_parameters_;
            int number_of_parameters_;
            int number_of_excluded_atoms_;
            int number_of_residues_;
            int total_number_of_bonds_;
            int total_number_of_angles_;
            int total_number_of_dihedrals_;
            int number_of_bond_types_;
            int number_of_angle_types_;
            int number_of_dihedral_types_;
            int number_of_atom_types_in_parameter_file_;
            int number_of_distinct_hydrogen_bonds_;
            int perturbation_option_;
            int number_of_bonds_perturbed_;
            int number_of_angles_perturbed_;
            int number_of_dihedrals_perturbed_;
            int number_of_bonds_group_perturbed_;
            int number_of_angles_group_perturbed_;
            int number_of_dihedrals_group_perturbed_;
            int standard_periodic_box_option_;
            int number_of_atoms_in_largest_residue_;
            int cap_option_;
            int number_of_extra_points_;
            int number_of_beads_;
            TopologyAtomTypeMap atom_types_;
            TopologyBondTypeMap bond_types_;
            TopologyAngleTypeMap angle_types_;
            TopologyDihedralTypeMap dihedral_types_;
            TopologyAssembly* assembly_;
            RadiusSet radius_set_;
    };
}

#endif // TOPOLOGYFILE_HPP
