#ifndef ASSEMBLY_HPP
#define ASSEMBLY_HPP

#include <string>
#include <iostream>
#include <vector>

#include "../Geometry/coordinate.hpp"
#include "../common.hpp"
#include "../FileSet/PdbFileSpace/pdbfile.hpp"
#include "../FileSet/TopologyFileSpace/topologyfile.hpp"
#include "../FileSet/CoordinateFileSpace/coordinatefile.hpp"
#include "../ParameterSet/PrepFileSpace/prepfile.hpp"
#include "../ParameterSet/LibraryFileSpace/libraryfile.hpp"
#include "../ParameterSet/ParameterFileSpace/parameterfile.hpp"
#include "../FileSet/PdbFileSpace/pdbmodelresidueset.hpp"
#include "../FileSet/PdbFileSpace/pdbmodelcard.hpp"
#include "../FileSet/PdbFileSpace/pdbmodel.hpp"

namespace MolecularModeling
{
    class Residue;
    class Atom;
    class Assembly
    {
        public:
            //////////////////////////////////////////////////////////
            //                    TYPE DEFINITION                   //
            //////////////////////////////////////////////////////////
            typedef std::vector<Assembly*> AssemblyVector;
            typedef std::vector<Residue*> ResidueVector;
            typedef std::vector<Atom*> AtomVector;
            typedef std::vector<Geometry::Coordinate*> CoordinateVector;

            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            Assembly();
            /*! \fn
              * Constructor to build a structure from the given set of input files
              * @param file_paths Set of file paths that are required to build a structure
              *         Most of the time it has just one element (topology, prep, lib, pdb) but at some point it needs more than one file (topology+coordinate)
              * @param type Type of the input which is selected from InputFileType enumerator
              */
            Assembly(std::vector<std::string> file_paths, gmml::InputFileType type);
            /*! \fn
              * Constructor to build a structure from multiple file types, a general version of the previous one
              * @param file_paths Set of set of file paths that are required to build a structure
              * @param types Set of inoput file types of the inputs which are selected from InputFileType enumerator
              */
            Assembly(std::vector<std::vector<std::string> > file_paths, std::vector<gmml::InputFileType> types);

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
            /*! \fn
              * An accessor function in order to access to the name
              * @return name_ attribute of the current object of this class
              */
            std::string GetName();
            /*! \fn
              * An accessor function in order to access to the assemblies
              * @return assemblies_ attribute of the current object of this class
              */
            AssemblyVector GetAssemblies();
            /*! \fn
              * An accessor function in order to access to the residues
              * @return residues_ attribute of the current object of this class
              */
            ResidueVector GetResidues();
            /*! \fn
              * An accessor function in order to access to the chemical type
              * @return chemical_type_ attribute of the current object of this class
              */
            std::string GetChemicalType();
            /*! \fn
              * An accessor function in order to access to the sequence number
              * @return sequence_number_ attribute of the current object of this class
              */
            int GetSequenceNumber();
            /*! \fn
              * An accessor function in order to access to the total mass
              * @return total_mass_ attribute of the current object of this class
              */
            double GetTotalMass();
            /*! \fn
              * An accessor function in order to access to the center of mass
              * @return center_of_mass_ attribute of the current object of this class
              */
            Geometry::Coordinate GetCenterOfMass();
            /*! \fn
              * An accessor function in order to access to the center of geometry
              * @return center_of_geometry_ attribute of the current object of this class
              */
            Geometry::Coordinate GetCenterOfGeometry();
            /*! \fn
              * An accessor function in order to access to the center of description
              * @return description_ attribute of the current object of this class
              */
            std::string GetDescription();
            /*! \fn
              * An accessor function in order to access to the source file
              * @return source_file_ attribute of the current object of this class
              */
            std::string GetSourceFile();
            /*! \fn
              * An accessor function in order to access to the source file type
              * @return source_file_type_ attribute of the current object of this class
              */
            gmml::InputFileType GetSourceFileType();
            /*! \fn
              * An accessor function in order to access to the model index
              * @return model_index_ attribute of the current object of this class
              */
            int GetModelIndex();
            /*! \fn
              * A functions that extracts all atoms of an assembly
              * @return Vector of all atoms in the current object of assembly
              */
            AtomVector GetAllAtomsOfAssembly();
            /*! \fn
              * A functions that extracts all residues of an assembly
              * @return Vector of all residues in the current object of assembly
              */
            ResidueVector GetAllResiduesOfAssembly();
            /*! \fn
              * A function to return all coordinates of all atoms in all residues and assemblies of an assembly
              * @return List of all coordinates of all atoms in all residues and assemblies of an assembly
              */
            CoordinateVector GetAllCoordinates();
            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A mutator function in order to set the name of the current object
              * Set the name_ attribute of the current assembly
              * @param name The name attribute of the current object
              */
            void SetName(std::string name);
            /*! \fn
              * A mutator function in order to set the assemblies of the current object
              * Set the assemblies_ attribute of the current assembly
              * @param assemblies The assemblies attribute of the current object
              */
            void SetAssemblies(AssemblyVector assemblies);
            /*! \fn
              * A function in order to add the assembly to the current object
              * Set the assemblies_ attribute of the current assembly
              * @param assembly The assembly of the current object
              */
            void AddAssembly(Assembly* assembly);
            /*! \fn
              * A mutator function in order to set the residues of the current object
              * Set the residues_ attribute of the current assembly
              * @param residues The residues attribute of the current object
              */
            void SetResidues(ResidueVector residues);
            /*! \fn
              * A function in order to add the residue to the current object
              * Set the residues_ attribute of the current assembly
              * @param residue The residue of the current object
              */
            void AddResidue(Residue* residue);
            /*! \fn
              * A mutator function in order to set the chemical type of the current object
              * Set the chemical_type_ attribute of the current assembly
              * @param chemical_type The chemical_type attribute of the current object
              */
            void SetChemicalType(std::string chemical_type);
            /*! \fn
              * A mutator function in order to set the sequence number of the current object
              * Set the sequence_number_ attribute of the current assembly
              * @param sequence_number The sequence number attribute of the current object
              */
            void SetSequenceNumber(int sequence_number);
            /*! \fn
              * A mutator function in order to set the total mass of the current object
              * Set the total_mass_ attribute of the current assembly
              * @param total_mass The total mass attribute of the current object
              */
            void SetTtoalMass(double total_mass);
            /*! \fn
              * A mutator function in order to set the center of mass of the current object
              * Set the center_of_mass_ attribute of the current assembly
              * @param center_of_mass The center of mass attribute of the current object
              */
            void SetCenterOfMass(Geometry::Coordinate center_of_mass);
            /*! \fn
              * A mutator function in order to set the center of geometry of the current object
              * Set the residues_ attribute of the current assembly
              * @param center_of_geometry The center of geometry attribute of the current object
              */
            void SetCenterOfGeometry(Geometry::Coordinate center_of_geometry);
            /*! \fn
              * A mutator function in order to set the description of the current object
              * Set the description_ attribute of the current assembly
              * @param description The description attribute of the current object
              */
            void SetDescription(std::string description);
            /*! \fn
              * A mutator function in order to set the source file of the current object
              * Set the source_file_ attribute of the current assembly
              * @param source_file The source file attribute of the current object
              */
            void SetSourceFile(std::string source_file);
            /*! \fn
              * A mutator function in order to set the source file type of the current object
              * Set the source_file_type_ attribute of the current assembly
              * @param source_file_type The source file type attribute of the current object
              */
            void SetSourceFileType(gmml::InputFileType source_file_type);
            /*! \fn
              * A mutator function in order to set the selected model index of the current object
              * Set the model_index_ attribute of the current assembly
              * @param model_index The target model index attribute of the current object
              */
            void SetModelIndex(int model_index);
            //////////////////////////////////////////////////////////
            //                       FUNCTIONS                      //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to build a structure from a single pdb file
              * Imports data from pdb file data structure into central data structure
              * @param pdb_file_path Path to a pdb file
              */
            void BuildAssemblyFromPdbFile(std::string pdb_file_path);
            /*! \fn
              * A function to build a structure from a single topology file
              * Imports data from topology file data structure into central data structure
              * @param topology_file_path Path to a topology file
              */
            void BuildAssemblyFromTopologyFile(std::string topology_file_path);
            /*! \fn
              * A function to build a structure from a single library file
              * Imports data from library file data structure into central data structure
              * @param library_file_path Path to a library file
              */
            void BuildAssemblyFromLibraryFile(std::string library_file_path);
            /*! \fn
              * A function to build a structure from a combination of a topology file and its corresponding coordinate file
              * Imports data from topology file data structure into central data structure and assign the atom coordinates based on the coordinate file
              * @param topology_file_path Path to a topology file
              * @param coordinate_file_path Path to a coordinate file corresponds to the given topology file
              */
            void BuildAssemblyFromTopologyCoordinateFile(std::string topology_file_path, std::string coordinate_file_path);
            /*! \fn
              * A function to build a structure from a single prep file
              * Imports data from prep file data structure into central data structure
              * @param prep_file_path Path to a prep file
              */
            void BuildAssemblyFromPrepFile(std::string prep_file_path);
            /*! \fn
              * A function to build a pdb file structure from the current assembly object
              * Exports data from assembly data structure into pdb file structure
              */
            PdbFileSpace::PdbFile* BuildPdbFileStructureFromAssembly();

            void ExtractPdbModelCardFromAssembly(PdbFileSpace::PdbModelResidueSet* residue_set, int &serial_number, int &sequence_number, int model_number = 0);
            /*! \fn
              * A function to extract bonds from the current assembly object
              * @param inserted_bond_types Bond types that have been already detected in the assembly structure
              * @param assembly_atom A source atom in a bond in the assembly structure
              * @param neighbor Second atom in a bond which is a neighbor of assembly_atom
              * @param bonds All known bonds from parameter file
              * @param bond_type_counter A counter that indicates the number of bond types that have been already detected and also determines the index associated with it
              * @param topology_file Output topology file structure that the detected bond types belong to
              */
            void ExtractTopologyBondTypesFromAssembly(std::vector<std::vector<std::string> > inserted_bond_types, Atom* assembly_atom, Atom* neighbor, ParameterFileSpace::ParameterFile::BondMap bonds, int bond_type_counter, TopologyFileSpace::TopologyFile* topology_file);
            /*! \fn
              * A function to extract bond types from the current assembly object
              * @param inserted_bonds Bonds that have been already detected in the assembly structure
              * @param inserted_bond_types Bond types that have been already detected in the assembly structure
              * @param assembly_atom A source atom in a bond in the assembly structure
              * @param neighbor Second atom in a bond which is a neighbor of assembly_atom
              * @param topology_file Output topology file structure that the detected bond types belong to
              */
            void ExtractTopologyBondsFromAssembly(std::vector<std::vector<std::string> > inserted_bonds, std::vector<std::vector<std::string> > inserted_bond_tyoes, Atom* assembly_atom, Atom* neighbor, TopologyFileSpace::TopologyFile* topology_file);
            /*! \fn
              * A function to build a topology file structure from the current assembly object
              * Exports data from assembly data structure into topology file structure
              */
            TopologyFileSpace::TopologyFile* BuildTopologyFileStructureFromAssembly(std::string parameter_file_path);
            /*! \fn
              * A function to build agnle types of topology file structure from the current assembly object
              * Exports data from assembly data structure to generate topology atom types
              * @param assembly_atom A source atom in an angle in the assembly structure
              * @param neighbor Second atom in an angle which is a neighbor of assembly_atom
              * @param neighbor_of_neighbor Third atom in an angle which is a neighbor of neighbor atom and it is not identical to assmebly_atom
              * @param inserted_angle_types Angle types that have been already detected in an assembly structure
              * @param angle_type_counter A counter that indicates the number of angle types that have been already detected and also determines the index associated with it
              * @param topology_file Output topology file structure that the detected angle types belong to
              * @param angles All known angles from parameter file
              */
            void ExtractTopologyAngleTypesFromAssembly(Atom* assembly_atom, Atom* neighbor, Atom* neighbor_of_neighbor,
                                                       std::vector<std::vector<std::string> > inserted_angle_types,
                                                       int angle_type_counter, TopologyFileSpace::TopologyFile* topology_file,
                                                       ParameterFileSpace::ParameterFile::AngleMap angles);
            /*! \fn
              * A function to build agnle types of topology file structure from the current assembly object
              * Exports data from assembly data structure to generate topology atom types
              * @param assembly_atom A source atom in an angle in the assembly structure
              * @param neighbor Second atom in an angle which is a neighbor of assembly_atom
              * @param neighbor_of_neighbor Third atom in an angle which is a neighbor of neighbor atom and it is not identical to assmebly_atom
              * @param inserted_angles Angles that have been already detected in an assembly structure
              * @param inserted_angle_types Angle types that have been already detected in an assembly structure
              * @param angle_type_counter A counter that indicates the number of angle types that have been already detected and also determines the index associated with it
              * @param topology_file Output topology file structure that the detected angle types belong to
              */
            void ExtractTopologyAnglesFromAssembly(Atom* assembly_atom, Atom* neighbor, Atom* neighbor_of_neighbor,
                                                             std::vector<std::vector<std::string> > inserted_angles,
                                                             std::vector<std::vector<std::string> > inserted_angle_types, TopologyFileSpace::TopologyFile* topology_file);
            /*! \fn
              * A function to build a coordinate file structure from the current assembly object
              * Exports data from assembly data structure into coordinate file structure
              */
            CoordinateFileSpace::CoordinateFile* BuildCoordinateFileStructureFromAssembly();
            /*! \fn
              * A function to build a library file structure from the current assembly object
              * Exports data from assembly data structure into library file structure
              */
            LibraryFileSpace::LibraryFile* BuildLibraryFileStructureFromAssembly();

            /*! \fn
              * A function to build a graph structure (bonding information) for the current object of central data structure
              * @param building_option A building option that can be selected from BuildingStructureOption enumerator
              * @param options List of additional options that can be defined by user
              * @param file_paths In the case that building structure based on the information of database files is selected this argument is a list of database file names
              */
            void BuildStructure(gmml::BuildingStructureOption building_option, std::vector<std::string> options, std::vector<std::string> file_paths);
            /*! \fn
              * A function to build a graph structure for the current object of central data structure based on the distance between the atoms of the structure
              * @param cutoff Threshold of closeness of the atoms to be considered as bonded
              * @param model_index In the case that the structure has multiple model (multiple coordinates for atoms, such as pdb) this arguments indicates the desired model index
              */
            void BuildStructureByDistance(double cutoff = gmml::dCutOff, int model_index = 0);
            /*! \fn
              * A function to build a graph structure for the current object of central data structure based on the bonding information provided in the original file
              */
            void BuildStructureByOriginalFileBondingInformation();
            /*! \fn
              * A function to build a graph structure for the current object of central data structure based on the bonding information provided in
              * the original file in the case that the original file is a pdb file
              */
            void BuildStructureByPDBFileInformation();
            /*! \fn
              * A function to build a graph structure for the current object of central data structure based on the bonding information provided in
              * the original file in the case that the original file is a topology file
              */
            void BuildStructureByTOPFileInformation();
            /*! \fn
              * A function to build a graph structure for the current object of central data structure based on the bonding information provided in
              * the original file in the case that the original file is a lib file
              */
            void BuildStructureByLIBFileInformation();
            /*! \fn
              * A function to build a graph structure for the current object of central data structure based on the bonding information provided in
              * the original file in the case that the original file is a prep file
              */
            void BuildStructureByPrepFileInformation();
            /*! \fn
              * A function to build a graph structure for the current object of central data structure based on the bonding information provided in
              * the given database files
              * @param types List of types of the database files
              * @param file_paths List of the database file paths
              */
            void BuildStructureByDatabaseFilesBondingInformation(std::vector<gmml::InputFileType> types, std::vector<std::string> file_paths);            
            /*! \fn
              * A function to calculate the center of geometry of the assembly
              */
            void CalculateCenterOfGeometry();
            /*! \fn
              * A function that counts the number of atoms in all assemblies and residues of the assembly
              * @return counter Number of atoms in all assemblies and residues in the current object of assembly
              */
            int CountNumberOfAtoms();            
            /*! \fn
              * A function that counts the number of atoms types in all assemblies and residues of the assembly
              * @return counter Number of atoms in all assemblies and residues in the current object of assembly
              */
            int CountNumberOfAtomTypes();
            /*! \fn
              * A function that counts the number of residues in all assemblies and residues of the assembly
              * @return type_list Number of residues in all assemblies and residues in the current object of assembly
              */
            int CountNumberOfResidues();
            /*! \fn
              * A function that counts the number of bonds including hydrogen in all assemblies and residues of the assembly
              * @return counter Number of bonds including hydrogen in all assemblies and residues in the current object of assembly
              */
            int CountNumberOfBondsIncludingHydrogen();
            /*! \fn
              * A function that counts the number of bonds excluding hydrogen in all assemblies and residues of the assembly
              * @return counter Number of bonds excluding hydrogen in all assemblies and residues in the current object of assembly
              */
            int CountNumberOfBondsExcludingHydrogen();
            /*! \fn
              * A function that counts the number of bonds in all assemblies and residues of the assembly
              * @return counter Number of bonds in all assemblies and residues in the current object of assembly
              */
            int CountNumberOfBonds();
            /*! \fn
              * A function that counts the number of bond types in all assemblies and residues of the assembly
              * @return type_list Number of bond types in all assemblies and residues in the current object of assembly
              */
            int CountNumberOfBondTypes();
            /*! \fn
              * A function that counts the number of angles including hydrogen in all assemblies and residues of the assembly
              * @return counter Number of angles including hydrogen in all assemblies and residues in the current object of assembly
              */
            int CountNumberOfAnglesIncludingHydrogen();
            /*! \fn
              * A function that counts the number of angles excluding hydrogen in all assemblies and residues of the assembly
              * @return counter Number of angles excluding hydrogen in all assemblies and residues in the current object of assembly
              */
            int CountNumberOfAnglesExcludingHydrogen();
            /*! \fn
              * A function that counts the number of angles in all assemblies and residues of the assembly
              * @return counter Number of angles in all assemblies and residues in the current object of assembly
              */
            int CountNumberOfAngles();
            /*! \fn
              * A function that counts the number of angle types in all assemblies and residues of the assembly
              * @return type_list Number of angle types in all assemblies and residues in the current object of assembly
              */
            int CountNumberOfAngleTypes();
            /*! \fn
              * A function that counts the number of dihedrals includig hydrogen in all assemblies and residues of the assembly
              * @param parameter_file_path The path of the parameter file
              * @return type_list Number of dihedrals including hydrogen in all assemblies and residues in the current object of assembly
              */
            int CountNumberOfDihedralsIncludingHydrogen(std::string parameter_file_path);
            /*! \fn
              * A function that counts the number of dihedrals excluding hydrogen in all assemblies and residues of the assembly
              * @param parameter_file_path The path of the parameter file
              * @return counter Number of dihedrals excluding hydrogen in all assemblies and residues in the current object of assembly
              */
            int CountNumberOfDihedralsExcludingHydrogen(std::string parameter_file_path);
            /*! \fn
              * A function that counts the number of dihedrals in all assemblies and residues of the assembly
              * @param parameter_file_path The path of the parameter file
              * @return counter Number of dihedrals in all assemblies and residues in the current object of assembly
              */
            int CountNumberOfDihedrals(std::string parameter_file_path);
            std::vector<std::vector<std::string> > CreateAllAtomTypePermutationsforDihedralType(std::string atom_type1, std::string atom_type2, std::string atom_type3, std::string atom_type4);
            std::vector<std::vector<std::string> > CreateAllAtomTypePermutationsforImproperDihedralType(std::string atom_type1, std::string atom_type2, std::string atom_type3, std::string atom_type4);
            /*! \fn
              * A function that counts the number of dihedral types in all assemblies and residues of the assembly
              * @param parameter_file_path The path of the parameter file
              * @return type_list Number of dihedral types in all assemblies and residues in the current object of assembly
              */
            int CountNumberOfDihedralTypes(std::string parameter_file_path);
            /*! \fn
              * A function that returns all atoms with at least three neighbors in all assemblies and residues of the assembly
              * @return atoms_with_at_least_three_neighbors A list of atoms with at least three neighbors in all assemblies and residues in the current object of assembly
              */
            AtomVector GetAllAtomsOfAssemblyWithAtLeastThreeNeighbors();
            /*! \fn
              * A function that counts the number of excluded atoms in all assemblies and residues of the assembly
              * @return counter Number of excluded atoms in all assemblies and residues in the current object of assembly
              */
            int CountNumberOfExcludedAtoms();
            /*! \fn
              * A function that counts the number of atoms in largest residue in all assemblies and residues of the assembly
              * @return max Number of atoms in largest residue in all assemblies and residues in the current object of assembly
              */
            int CountMaxNumberOfAtomsInLargestResidue();

            void ClearAssembly();
            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the assembly contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cout);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            std::string name_;                              /*!< Assembly name >*/
            AssemblyVector assemblies_;                     /*!< List of assemblies that may involved in a big structure as a new assembly >*/
            ResidueVector residues_;                        /*!< List of residues involved in the current object of assembly >*/
            std::string chemical_type_;                     /*!< A descriptor for the chemical type of the current object of assembly >*/
            int sequence_number_;                           /*!< An integer number to indicates sequence of importing the current assembly >*/
            double total_mass_;                             /*!< Total mass of the assembly >*/
            Geometry::Coordinate center_of_mass_;           /*!< Center of mass of the assembly >*/
            Geometry::Coordinate center_of_geometry_;       /*!< Center of geometry of the assembly >*/
            std::string description_;                       /*!< Short description for the current assembly >*/
            std::string source_file_;                       /*!< File name that the current assembly has been built upon >*/
            gmml::InputFileType source_file_type_;          /*!< Type of the file that the current assembly has been built upon >*/
            int model_index_;                               /*!< In case that there are more than one models for an assembly, this attribute indicated which model is the target model >*/
    };
}

#endif // ASSEMBLY_HPP
