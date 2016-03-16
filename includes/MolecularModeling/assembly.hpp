#ifndef ASSEMBLY_HPP
#define ASSEMBLY_HPP

#include <string>
#include <iostream>
#include <vector>

#include "../GeometryTopology/coordinate.hpp"
#include "../GeometryTopology/plane.hpp"
#include "../common.hpp"
#include "../Glycan/chemicalcode.hpp"
#include "../Glycan/sugarname.hpp"
#include "../Glycan/monosaccharide.hpp"
#include "../Glycan/ontologyvocabulary.hpp"
#include "../InputSet/PdbFileSpace/pdbfile.hpp"
#include "../InputSet/PdbqtFileSpace/pdbqtfile.hpp"
#include "../InputSet/TopologyFileSpace/topologyfile.hpp"
#include "../InputSet/CoordinateFileSpace/coordinatefile.hpp"
#include "../ParameterSet/PrepFileSpace/prepfile.hpp"
#include "../ParameterSet/PrepFileSpace/prepfileresidue.hpp"
#include "../ParameterSet/PrepFileSpace/prepfileatom.hpp"
#include "../ParameterSet/LibraryFileSpace/libraryfile.hpp"
#include "../ParameterSet/ParameterFileSpace/parameterfile.hpp"
#include "../InputSet/PdbqtFileSpace/pdbqtmodelresidueset.hpp"
#include "../InputSet/PdbFileSpace/pdbmodelresidueset.hpp"
#include "../InputSet/PdbFileSpace/pdbmodelcard.hpp"
#include "../InputSet/PdbFileSpace/pdbmodel.hpp"
#include "../Glycan/oligosaccharide.hpp"
#include "../Glycan/note.hpp"
#include "../InputSet/CondensedSequenceSpace/condensedsequence.hpp"

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
            typedef std::vector<GeometryTopology::Coordinate*> CoordinateVector;
            typedef std::map<std::string, gmml::GraphSearchNodeStatus> AtomStatusMap;
            typedef std::map<std::string, Atom*> AtomIdAtomMap;
//            typedef std::vector<AtomVector> AtomVectorVector;
            typedef std::map<std::string, AtomVector> CycleMap;
            typedef std::map<std::string, std::map<std::string, std::vector<std::string> > > SelectPatternMap;
            typedef std::map<std::string, ResidueVector> HierarchicalContainmentMap;
//            typedef std::map<Residue*, ResidueVector> ResidueAttachmentMap;
            typedef std::vector<Glycan::Oligosaccharide*> OligosaccharideVector;
            typedef std::map<int, int> AssemblytoPdbSequenceNumberMap;
            typedef std::map<int, int> AssemblytoPdbSerialNumberMap;
            typedef std::map<std::string, std::string> DerivativeModificationMap;
            typedef std::vector<std::vector<std::string> > AttachedGlycanStructuresVector;
            typedef std::vector<Glycan::Note*> NoteVector;

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
            Assembly(Assembly* assembly);

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
           void ODBC();

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
            std::string GetId();
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
              * A functions that extracts all atoms of an assembly except atoms of water residues
              * @return Vector of all atoms in the current object of assembly except atoms of water residues
              */
            AtomVector GetAllAtomsOfAssemblyExceptProteinWaterResiduesAtoms();
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
            /*! \fn
              * A function to return all issues/notes within an assembly
              * @return List of all notes of an assembly
              */
            NoteVector GetNotes();

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
            void UpdateIds(std::string new_id);
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
            void SetId(std::string id);
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
            /*! \fn
              * A mutator function in order to set the list of notes of the current object
              * Set the notes_ attribute of the current assembly
              * @param notes The notes attribute of the current object
              */
            void SetNotes(NoteVector notes);
            /*! \fn
              * A function in order to add a note instance to the current object
              * Set the notes_ attribute of the current assembly
              * @param note The note instance of the current object
              */
            void AddNote(Glycan::Note* note);

            //////////////////////////////////////////////////////////
            //                       FUNCTIONS                      //
            //////////////////////////////////////////////////////////
            void BuildAssemblyFromCondensedSequence(std::string sequence, std::string prep_file, std::string parameter_file, bool structure = false);
            void AttachResidues(Residue* residue, Residue* parent_residue, int branch_index, std::string parameter_file);
            void SetAttachedResidueBond(Residue* residue, Residue* parent_residue, int branch_index, std::string parameter_file);
            void SetAttachedResidueAngle(Residue* residue, Residue* parent_residue, int branch_index, std::string parameter_file);
            void SetAttachedResidueTorsion(Residue* residue, Residue* parent_residue, int branch_index);
            void SetDihedral(Atom* atom1, Atom* atom2, Atom* atom3, Atom* atom4, double torsion);
            /*! \fn
              * A function to build a structure from a single pdb file
              * Imports data from pdb file data structure into central data structure
              * @param pdb_file_path Path to a pdb file
              */
            void BuildAssemblyFromPdbFile(std::string pdb_file_path, std::vector<std::string> amino_lib_files = std::vector<std::string>(),
                                          std::vector<std::string> glycam_lib_files = std::vector<std::string>(),
                                          std::vector<std::string> other_lib_files = std::vector<std::string>(),
                                          std::vector<std::string> prep_files = std::vector<std::string>(),
                                          std::string parameter_file = "");
            /*! \fn
              * A function to build a structure from a single pdb file
              * Imports data from pdb file data structure into central data structure
              * @param pdb_file Pdb file object
              */
            void BuildAssemblyFromPdbFile(PdbFileSpace::PdbFile* pdb_file,
                                          std::vector<std::string> amino_lib_files = std::vector<std::string>(),
                                          std::vector<std::string> glycam_lib_files = std::vector<std::string>(),
                                          std::vector<std::string> other_lib_files = std::vector<std::string>(),
                                          std::vector<std::string> prep_files = std::vector<std::string>(),
                                          std::string parameter_file = "");
            /*! \fn
              * A function to build a structure from a single pdbqt file
              * Imports data from pdbqt file data structure into central data structure
              * @param pdbqt_file_path Path to a pdbqt file
              */
            void BuildAssemblyFromPdbqtFile(std::string pdbqt_file_path, std::string parameter_file = "");
            /*! \fn
              * A function to build a structure from a single pdbqt file
              * Imports data from pdbqt file data structure into central data structure
              * @param pdbqt_file Pdbqt file object
              */
            void BuildAssemblyFromPdbqtFile(PdbqtFileSpace::PdbqtFile* pdbqt_file, std::string parameter_file = "");
            /*! \fn
              * A function to build a structure from a single topology file
              * Imports data from topology file data structure into central data structure
              * @param topology_file_path Path to a topology file
              */
            void BuildAssemblyFromTopologyFile(std::string topology_file_path, std::string parameter_file = "");
            /*! \fn
              * A function to build a structure from a single topology file
              * Imports data from topology file data structure into central data structure
              * @param topology_file Topology file object
              */
            void BuildAssemblyFromTopologyFile(TopologyFileSpace::TopologyFile* topology_file, std::string parameter_file = "");
            /*! \fn
              * A function to build a structure from a single library file
              * Imports data from library file data structure into central data structure
              * @param library_file_path Path to a library file
              */
            void BuildAssemblyFromLibraryFile(std::string library_file_path, std::string parameter_file = "");
            /*! \fn
              * A function to build a structure from a single library file
              * Imports data from library file data structure into central data structure
              * @param library_file Library file object
              */
            void BuildAssemblyFromLibraryFile(LibraryFileSpace::LibraryFile* library_file, std::string parameter_file = "");
            /*! \fn
              * A function to build a structure from a combination of a topology file and its corresponding coordinate file
              * Imports data from topology file data structure into central data structure and assign the atom coordinates based on the coordinate file
              * @param topology_file_path Path to a topology file
              * @param coordinate_file_path Path to a coordinate file corresponding to the given topology file
              */
            void BuildAssemblyFromTopologyCoordinateFile(std::string topology_file_path, std::string coordinate_file_path, std::string parameter_file = "");
            /*! \fn
              * A function to build a structure from a combination of a topology file and its corresponding coordinate file
              * Imports data from topology file data structure into central data structure and assign the atom coordinates based on the coordinate file
              * @param topology_file Topology file object
              * @param coordinate_file Coordinate file object corresponding to the given topology file
              */
            void BuildAssemblyFromTopologyCoordinateFile(TopologyFileSpace::TopologyFile* topology_file, CoordinateFileSpace::CoordinateFile* coordinate_file,
                                                         std::string parameter_file = "");
            /*! \fn
              * A function to build a structure from a single prep file
              * Imports data from prep file data structure into central data structure
              * @param prep_file_path Path to a prep file
              */
            void BuildAssemblyFromPrepFile(std::string prep_file_path, std::string parameter_file = "");
            /*! \fn
              * A function to build a structure from a single prep file
              * Imports data from prep file data structure into central data structure
              * @param prep_file Prep file object
              */
            void BuildAssemblyFromPrepFile(PrepFileSpace::PrepFile* prep_file, std::string parameter_file = "");
//            PdbFileSpace::PdbFile* BuildPdbFileStructureFromCondensedSequence
            /*! \fn
              * A function to build a pdb file structure from the current assembly object
              * Exports data from assembly data structure into pdb file structure
              * @param link_card_direction An integer to define the direction in the link cards (-1: O -> C and 1: C -> O)
              */
            PdbFileSpace::PdbFile* BuildPdbFileStructureFromAssembly(int link_card_direction = -1, int connect_card_existance = 1);
            /*! \fn
              * A function to build a pdbqt file structure from the current assembly object
              * Exports data from assembly data structure into pdbqt file structure
              */
            PdbqtFileSpace::PdbqtFile* BuildPdbqtFileStructureFromAssembly();

            void ExtractPdbModelCardFromAssembly(PdbFileSpace::PdbModelResidueSet* residue_set, int &serial_number, int &sequence_number, int model_number,
                                                 AssemblytoPdbSequenceNumberMap& assembly_to_pdb_sequence_number_map, AssemblytoPdbSerialNumberMap& assembly_to_pdb_serial_number_map);
            void ExtractPdbqtModelCardFromAssembly(PdbqtFileSpace::PdbqtModelResidueSet* residue_set, int &serial_number, int &sequence_number, int model_number);
            void ExtractPdbLinkCardFromAssembly(PdbFileSpace::PdbLinkCard* link_card, int model_index, AssemblytoPdbSequenceNumberMap assembly_to_pdb_sequence_number,
                                                int link_card_direction = -1);
            void ExtractPdbConnectCardFromAssembly(PdbFileSpace::PdbConnectCard* connect_card,
                                                   AssemblytoPdbSerialNumberMap assembly_to_pdb_serial_number);
            /*! \fn
              * A function to extract bonds from the current assembly object
              * @param inserted_bond_types Bond types that have been already detected in the assembly structure
              * @param assembly_atom A source atom in a bond in the assembly structure
              * @param neighbor Second atom in a bond which is a neighbor of assembly_atom
              * @param bonds All known bonds from parameter file
              * @param bond_type_counter A counter that indicates the number of bond types that have been already detected and also determines the index associated with it
              * @param topology_file Output topology file structure that the detected bond types belong to
              */
            void ExtractTopologyBondTypesFromAssembly(std::vector<std::vector<std::string> > &inserted_bond_types, Atom* assembly_atom, Atom* neighbor,
                                                      ParameterFileSpace::ParameterFile::BondMap &bonds, int &bond_type_counter, TopologyFileSpace::TopologyFile* topology_file);
            /*! \fn
              * A function to extract bond types from the current assembly object
              * @param inserted_bonds Bonds that have been already detected in the assembly structure
              * @param inserted_bond_types Bond types that have been already detected in the assembly structure
              * @param assembly_atom A source atom in a bond in the assembly structure
              * @param neighbor Second atom in a bond which is a neighbor of assembly_atom
              * @param topology_file Output topology file structure that the detected bond types belong to
              */
            void ExtractTopologyBondsFromAssembly(std::vector<std::vector<std::string> > &inserted_bonds, std::vector<std::vector<std::string> > &inserted_bond_tyoes,
                                                  Atom* assembly_atom, Atom* neighbor, TopologyFileSpace::TopologyFile* topology_file);

            /*! \fn
              * A function to build a prep file structure from the current assembly object
              * Exports data from assembly data structure into prep file structure
              * @param parameter_file_path A path to a parameter file
              */
            PrepFileSpace::PrepFile* BuildPrepFileStructureFromAssembly(std::string parameter_file_path);
            /*! \fn
              * A function to build improper dihedral types of prep file structure from the current assembly object
              * Exports data from assembly data structure to generate prep improper dihedral types
              * @param assembly_atom A source atom in a dihedral in the assembly structure
              * @param inserted_improper_dihedrals_types Improper dihedral types that have been already detected in an assembly structure
              * @param dihedrals All known dihedrals from parameter file
              */
            void ExtractPrepImproperDihedralTypesFromAssembly(Atom* assembly_atom, std::vector<std::string> &inserted_improper_dihedral_types, ParameterFileSpace::ParameterFile::DihedralMap &dihedrals);
            /*! \fn
              * A function to build dihedrals of prep file structure from the current assembly object
              * Exports data from assembly data structure to generate prep dihedrals
              * @param assembly_atom A source atom in a dihedral in the assembly structure
              * @param inserted_improper_dihedrals Improper dihedrals that have been already detected in an assembly structure
              * @param inserted_improper_dihedral_types Improper dihedral types that have been already detected in an assembly structure
              * @param dihedrals All known dihedrals from parameter file
              */
            void ExtractPrepImproperDihedralsFromAssembly(Atom* assembly_atom, std::vector<std::vector<std::string> > &inserted_improper_dihedrals, std::vector<std::string> inserted_improper_dihedral_types,
                                                          ParameterFileSpace::ParameterFile::DihedralMap &dihedrals);
            /*! \fn
              * A function to extract topological types of atoms of a residue in assembly object
              * @param assembly_atoms Atoms of a residue in the assembly structure
              * @return List of all topological types of atoms in a residue of an assembly
              */
            std::vector<gmml::TopologicalType> GetAllTopologicalTypesOfAtomsOfResidue(AtomVector assembly_atoms,
                                                                                      PrepFileSpace::PrepFileResidue::Loop& loops,
                                                                                      std::vector<int>& bond_index,
                                                                                      int dummy_atoms = gmml::DEFAULT_DUMMY_ATOMS);
            /*! \fn
              * A function to build a topology file structure from the current assembly object
              * Exports data from assembly data structure into topology file structure
              */
            TopologyFileSpace::TopologyFile* BuildTopologyFileStructureFromAssembly(std::string parameter_file_path, std::string ion_parameter_file_path = "");
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
                                                       std::vector<std::vector<std::string> > &inserted_angle_types,
                                                       int &angle_type_counter, TopologyFileSpace::TopologyFile* topology_file,
                                                       ParameterFileSpace::ParameterFile::AngleMap &angles);
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
                                                             std::vector<std::vector<std::string> > &inserted_angles,
                                                             std::vector<std::vector<std::string> > &inserted_angle_types, TopologyFileSpace::TopologyFile* topology_file);
            /*! \fn
              * A function to build dihedral types of topology file structure from the current assembly object
              * Exports data from assembly data structure to generate topology dihedral types
              * @param assembly_atom A source atom in a dihedral in the assembly structure
              * @param neighbor Second atom in a dihedral which is a neighbor of assembly_atom
              * @param neighbor_of_neighbor Third atom in a dihedral which is a neighbor of neighbor atom and it is not identical to assmebly_atom
              * @param neighbor_of_neighbor_of_neighbor Fourth atom in a dihedral which is a neighbor of neighbor of neighbor atom and it is not identical to neighbor atom
              * @param inserted_dihedrals_types Dihedral types that have been already detected in an assembly structure
              * @param dihedral_type_counter A counter that indicates the number of dihedral types that have been already detected and also determines the index associated with it
              * @param topology_file Output topology file structure that the detected angle types belong to
              * @param dihedrals All known dihedrals from parameter file
              */
            void ExtractTopologyDihedralTypesFromAssembly(Atom* assembly_atom, Atom* neighbor, Atom* neighbor_of_neighbor, Atom* neighbor_of_neighbor_of_neighbor,
                                                       std::vector<std::string> &inserted_dihedral_types,
                                                       int &dihedral_type_counter, TopologyFileSpace::TopologyFile* topology_file,
                                                       ParameterFileSpace::ParameterFile::DihedralMap &dihedrals);
            /*! \fn
              * A function to build dihedrals of topology file structure from the current assembly object
              * Exports data from assembly data structure to generate topology dihedrals
              * @param assembly_atom A source atom in a dihedral in the assembly structure
              * @param neighbor Second atom in a dihedral which is a neighbor of assembly_atom
              * @param neighbor_of_neighbor Third atom in a dihedral which is a neighbor of neighbor atom and it is not identical to assmebly_atom
              * @param neighbor_of_neighbor_of_neighbor Fourth atom in a dihedral which is a neighbor of neighbor of neighbor atom and it is not identical to neighbor atom
              * @param inserted_dihedrals Dihedrals that have been already detected in an assembly structure
              * @param inserted_dihedral_types Dihedral types that have been already detected in an assembly structure
              * @param dihedrals All known dihedrals from parameter file
              * @param topology_file Output topology file structure that the detected angle types belong to
              */
            void ExtractTopologyDihedralsFromAssembly(Atom* assembly_atom, Atom* neighbor, Atom* neighbor_of_neighbor, Atom* neighbor_of_neighbor_of_neighbor,
                                                             std::vector<std::vector<std::string> > &inserted_dihedrals, std::vector<std::string> &inserted_dihedral_types,
                                                             ParameterFileSpace::ParameterFile::DihedralMap &dihedrals,
                                                             TopologyFileSpace::TopologyFile* topology_file);

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
            void BuildStructureByDistance(int number_of_threads = 1, double cutoff = gmml::dCutOff, int model_index = 0);
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
              * @param parameter_file_path The path of the parameter file
              * @return counter Number of bonds including hydrogen in all assemblies and residues in the current object of assembly
              */
            int CountNumberOfBondsIncludingHydrogen(std::string parameter_file_path);
            /*! \fn
              * A function that counts the number of bonds excluding hydrogen in all assemblies and residues of the assembly
              * @param parameter_file_path The path of the parameter file
              * @return counter Number of bonds excluding hydrogen in all assemblies and residues in the current object of assembly
              */
            int CountNumberOfBondsExcludingHydrogen(std::string parameter_file_path);
            /*! \fn
              * A function that counts the number of bonds in all assemblies and residues of the assembly
              * @return counter Number of bonds in all assemblies and residues in the current object of assembly
              */
            int CountNumberOfBonds();
            /*! \fn
              * A function that counts the number of bond types in all assemblies and residues of the assembly
              * @param parameter_file_path The path of the parameter file
              * @return type_list Number of bond types in all assemblies and residues in the current object of assembly
              */
            int CountNumberOfBondTypes(std::string parameter_file_path);
            /*! \fn
              * A function that counts the number of angles including hydrogen in all assemblies and residues of the assembly
              * @param parameter_file_path The path of the parameter file
              * @return counter Number of angles including hydrogen in all assemblies and residues in the current object of assembly
              */
            int CountNumberOfAnglesIncludingHydrogen(std::string parameter_file_path);
            /*! \fn
              * A function that counts the number of angles excluding hydrogen in all assemblies and residues of the assembly
              * @param parameter_file_path The path of the parameter file
              * @return counter Number of angles excluding hydrogen in all assemblies and residues in the current object of assembly
              */
            int CountNumberOfAnglesExcludingHydrogen(std::string parameter_file_path);
            /*! \fn
              * A function that counts the number of angles in all assemblies and residues of the assembly
              * @return counter Number of angles in all assemblies and residues in the current object of assembly
              */
            int CountNumberOfAngles();
            /*! \fn
              * A function that counts the number of angle types in all assemblies and residues of the assembly
              * @param parameter_file_path The path of the parameter file
              * @return type_list Number of angle types in all assemblies and residues in the current object of assembly
              */
            int CountNumberOfAngleTypes(std::string parameter_file_path);
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

            /*! \fn
              * A function to select a set of atoms in an assembly by using a string pattern
              * @param pattern Pattern of selection inclusing assembly number, residue name and/or atom name
              * Sample pattern: "1.*:#520,MAN,GAL@#3740-3750,^C:ALA@O1;1.1.2,1.1.3:NAG@O$"
              * @return A list of atoms which satisfy the condition specified in the pattern
              */
            AtomVector Select(std::string pattern);
            SelectPatternMap ParsePatternString(std::string pattern);
            void GetHierarchicalMapOfAssembly(HierarchicalContainmentMap& hierarchical_map, std::stringstream& index);

            void ClearAssembly();

//            void CycleDetection();
//            std::vector<std::vector<std::string> > CreateAllCyclePermutations(std::string id1, std::string id2, std::string id3, std::string id4, std::string id5, std::string id6);

//            CycleMap DetectCyclesByDFS(std::string cycle_size = "5|6");

            /*! \fn
              * A function in order to access to the list of all residue from lib files
              * @param lib_files The list of paths to library files
              * @return all_residues_
              */
            LibraryFileSpace::LibraryFile::ResidueMap GetAllResiduesFromMultipleLibFilesMap(std::vector<std::string> lib_files);
            /*! \fn
              * A function in order to access to the list of all residue from prep files
              * @param prep_files The list of paths to prep files
              * @return all_residues_
              */
            PrepFileSpace::PrepFile::ResidueMap GetAllResiduesFromMultiplePrepFilesMap(std::vector<std::string> prep_files);

            /*! \fn
              * A function in order to extract all residue names existing in the given lib files
              * @param lib_files The list of paths to library files
              * @return all_residue_names
              */
            gmml::ResidueNameMap GetAllResidueNamesFromMultipleLibFilesMap(std::vector<std::string> lib_files);

            gmml::GlycamResidueNamingMap ExtractResidueGlycamNamingMap(OligosaccharideVector oligosaccharides);
            void ExtractOligosaccharideNamingMap(gmml::GlycamResidueNamingMap& pdb_glycam_map, Glycan::Oligosaccharide* oligosaccharide,
                                                 CondensedSequenceSpace::CondensedSequence::CondensedSequenceAmberPrepResidueTree condensed_sequence_amber_residue_tree,
                                                int& index);            
            void UpdateResidueName2GlycamName(gmml::GlycamResidueNamingMap residue_glycam_map, std::string prep_file);
            bool PatternMatching(Residue* residue, Residue* query_residue, gmml::GlycamAtomNameMap& pdb_glycam_map, gmml::GlycamAtomNameMap& glycam_atom_map);
            bool PatternMatching(Atom* atom, Atom* query_atom, gmml::GlycamAtomNameMap& pdb_glycam_map, gmml::GlycamAtomNameMap& glycam_atom_map);
            bool HasAllNeighborsOf(Atom* atom, Atom* query_atom);

            /*! \fn
              * A function in order to extract all the saccharide structures
              * @param amino_lib_files The list of paths to amino library files
              */
            OligosaccharideVector ExtractSugars(std::vector<std::string> amino_lib_files);            
            /*! \fn
            * A function in order to detec the shape of the ring using the external BFMP program
            * This function creates a pdb file and a configuration file for input arguments of the external detect_shape program.
            * the function updates the bfmp_ring_confomration attribute of the monosaccharide
            * @param cycle The list of ring atoms
            * @param mono The monosaccharide object
            */
            void DetectShape(AtomVector cycle, Glycan::Monosaccharide* mono);

            void PopulateOntology(std::ofstream& out_file, OligosaccharideVector oligosaccharides);
            void PopulateNotes(std::stringstream& pdb_stream, std::stringstream& note_stream, std::string pdb_uri, NoteVector notes, std::string id_prefix, int note_id);
            void PopulateOligosaccharide(std::stringstream& pdb_stream, std::stringstream& oligo_stream, std::stringstream& mono_stream, std::stringstream& linkage_stream, std::string pdb_uri,
                                         std::string id_prefix, int& link_id, OligosaccharideVector oligos, std::vector<std::string>& side_or_ring_atoms,
                                         std::vector<int>& visited_oligos);
            void PopulateLinkage(std::stringstream& linkage_stream, Glycan::Oligosaccharide* oligo, std::string oligo_uri, std::string id_prefix, int& link_id,
                                 std::vector<int>& visited_oligos);
            void PopulateMonosaccharide(std::stringstream& pdb_stream, std::stringstream& oligo_stream, std::string oligo_uri, std::string id_prefix,
                                        Glycan::Monosaccharide* mono, std::vector<std::string>& side_or_ring_atoms);
            void PopulateRingAtom(std::stringstream& ring_atom_stream, std::string id_prefix, std::string ring_uri, std::string ring_resource, int ring_index, Atom* ring_atom,
                                  Glycan::Monosaccharide* mono, std::vector<std::string>& side_or_ring_atoms);
            void PopulateSideAtom(std::stringstream& side_atom_stream, std::string id_prefix, std::string side_uri, std::string side_resource, int ring_index, int side_index, Atom* side_atom,
                                  Glycan::Monosaccharide* mono, std::vector<std::string>& side_or_ring_atoms);
            void PopulateSugarName(std::stringstream& mono_stream, std::string mono_uri, std::string id_prefix, int mono_id, Glycan::SugarName sugar_name);
            void PopulateResidue(std::stringstream& pdb_stream, std::stringstream& residue_stream, std::string pdb_uri, std::string id_prefix, ResidueVector residues,
                                 std::vector<std::string> side_or_ring_atoms);
            void PopulateAtom(std::stringstream& atom_stream, std::string atom_uri, std::string atom_resource, std::string id_prefix, Atom* atom);
            void CreateTitle(std::string pdb_resource, std::stringstream& pdb_stream);
            void AddTriple(std::string s, std::string p, std::string o, std::stringstream& stream);
            void AddLiteral(std::string s, std::string p, std::string o, std::stringstream& stream);
            std::string CreateURIResource(gmml::URIType resource, int number, std::string id_prefix, std::string id);
            void FormulateCURL(std::string query);
            std::string CreateURI(std::string uri_resource);
            std::string ExtractOntologyInfoByNameOfGlycan(std::string stereo_name, std::string stereo_condensed_name, std::string name, std::string condensed_name);
            std::string ExtractOntologyInfoByNamePartsOfGlycan(std::string isomer, std::string ring_type, std::string configuration);
            std::string ExtractOntologyInfoByPDBID(std::string pdb_id);
            std::string ExtractOntologyInfoByStringChemicalCode(std::string chemical_code);
            std::string ExtractOntologyInfoByOligosaccharideNameSequence(std::string oligo_name);
            std::string ExtractOntologyInfoByOligosaccharideNameSequenceByRegex(std::string oligo_name_pattern);
            std::string ExtractOntologyInfoByGlycanStructure(std::string ring_type, std::string anomeric_orientation, std::string minus_one_orientation, std::string index_two_orientation,
                                                               std::string index_three_orientation, std::string index_four_orientation = "", std::string plus_one_orientation = "");
            std::string ExtractOntologyInfoByDerivativeModificationMap(std::string ring_type, DerivativeModificationMap derivative_modification_map);
            std::string ExtractOntologyInfoByAttachedGlycanStructures(AttachedGlycanStructuresVector attached_structures);

            void TestQueries();

            /*! \fn
              * A function in order to extract and print out all saccharides ring atoms information
              */
            void ExtractRingAtomsInformation();
            /*! \fn
              * A function in order to detect cycles in the molecular graph using the exhaustive ring perception algorithm
              * @return cycles A map between the string version of atoms of cycles and the list of cycle atom objects
              */
            CycleMap DetectCyclesByExhaustiveRingPerception();
            /*! \fn
              * A function in order to prune the graph (recursively removing nodes with zero or 1 neighbors)
              * @param all_atoms The list of atoms of the graph
              */
            void PruneGraph(AtomVector& all_atoms);
            /*! \fn
              * A function in order to convert the graph into a path graph (creating list of edges between the nodes and a list of labels for those edges )
              * @param path_graph_edges The list of edges between the nodes to be filled by the function
              * @param path_graph_labels The list of edge labels to be filled by the function
              * @param atoms The list of atoms of the graph
              */
            void ConvertIntoPathGraph(std::vector<std::string>& path_graph_edges, std::vector<std::string>& path_graph_labels, AtomVector atoms);
            /*! \fn
              * A function in order to reduce the path graph such that if there is a path/walk a-b-c in the graph converting it to a-c and creating a new label
                    for the new edge and checking if the new edge makes a cycle
              * @param path_graph_edges The list of edges between the nodes
              * @param path_graph_labels The list of edge labels
              * @param reduced_path_graph_edges The list of edges of the (reduced) path graph
              * @param reduced_path_graph_labels The list of edge labels of the (reduced) path graph
              * @param common_atom The atom that needs to be checked if it is involved in a walk (a path like a-b-c)
              * @param cycles The list of cycles to be filled by the function
              */
            void ReducePathGraph(std::vector<std::string> path_graph_edges, std::vector<std::string> path_graph_labels,
                                 std::vector<std::string>& reduced_path_graph_edges, std::vector<std::string>& reduced_path_graph_labels, std::string common_atom, std::vector<std::string>& cycles);

            /*! \fn
              * A function in order to detect cycles in the molecular graph using depth first search algorithm
              * @return cycles A map between the string version of atoms of cycles and the list of cycle atom objects
              */
            CycleMap DetectCyclesByDFS();
            /*! \fn
              * A function of depth first search algorithm in order to traverse the graph
              * @param path_graph_edges The list of edges between the nodes to be filled by the function
              * @param path_graph_labels The list of edge labels to be filled by the function
              * @param atoms The list of atoms of the graph
              */
            void DFSVisit(AtomVector atoms, AtomStatusMap& atom_status_map, AtomIdAtomMap& atom_parent_map, Atom* atom, int& counter, AtomIdAtomMap& dest_srd_map);
            /*! \fn
              * A recursive function in order to back track a path from current atom to a source atom to return the atom objects using the information extracted by the DFS algorithm
              * @param src_id The identifier of the source atom
              * @param current_atom The current atom that the function has been called on
              * @param atom_parent_map A map between atom identifier and its parent atom object
              * @param cycle The list of atom objects involved in a cycle
              * @param cycle_stream The cycle path that has been back traversed so far
              */
            void ReturnCycleAtoms(std::string src_id, Atom* current_atom, AtomIdAtomMap& atom_parent_map, AtomVector& cycle, std::stringstream& cycle_stream);
            /*! \fn
              * A function in order to discard rings/cycles that are only made from carbons atoms
              * @param cycles A map between the string version of atoms of cycles and the list of cycle atom objects
              */
            void FilterAllCarbonCycles(CycleMap& cycles);
            /*! \fn
              * A function in order to discard cycles/rings that are sharing an edge (fused cycles)
              * @param cycles A map between the string version of atoms of cycles and the list of cycle atom objects
              */
            void RemoveFusedCycles(CycleMap& cycles);
            /*! \fn
              * A function in order to detect the anomeric carbon of the ring (the carbon which has two oxygon neighbors)
              * @param anomeric_carbons_status The detection status of the anomeric carbon to be filled by the function
              * @param cycle The list of cycle atoms
              * @param cycle_atom_str The string version of atom identifiers of the cycle
              * @return o_neighbor The anomeric carbon which is the neighbor of oxygen
              */
            Atom* FindAnomericCarbon(Glycan::Note* anomeric_carbons_note, std::vector<std::string>& anomeric_carbons_status, AtomVector cycle, std::string cycle_atoms_str);
            /*! \fn
              * A function in order to sort atom objects of the cycle starting from the anomeric carbon of the ring (ring oxygen will be last atom)
              * @param cycle The list of cycle atoms
              * @param anomeric_atom The anomeric carbon of the ring
              * @param sorted_cycle_stream The sorted atom of the cycle so far (to be filled with the fuction)
              * @return sorted_cycle The sorted list of cycle atom objects
              */
            AtomVector SortCycle(AtomVector cycle, Atom* anomeric_atom, std::stringstream& sorted_cycle_stream);
            /*! \fn
              * A function in order to extract the relative oriantation of the side atoms that are attached to ring atoms against the ring
              * @param mono The monosaccharide object
              * @param sorted_cycle_stream The sorted atom of the cycle so far (to be filled with the fuction)
              * @param cycle_atom_str The string version of atom identifiers of the cycle
              * @return orientations The list of side atoms orinetations
              */
            std::vector<std::string> GetSideGroupOrientations(Glycan::Monosaccharide* mono, std::string cycle_atoms_str);
            /*! \fn
              * A function in order to build a chemical code structure for the monosaccharide based on the side atom orientations
              * @param orientations The list of side atoms orinetations
              * @return code The chemical code of the monosaccharide structure
              */
            Glycan::ChemicalCode* BuildChemicalCode(std::vector<std::string> orientations);
            /*! \fn
              * A function in order to extract additional side atoms of the saccharide (+2 and +3 positions)
              * @param mono The monosaccharide object
              * @return plus_sides The list of additional side atoms of the monosaccharide
              */
            AtomVector ExtractAdditionalSideAtoms(Glycan::Monosaccharide* mono);
            /*! \fn
              * A function in order to extract the probable derivatives that are attached to the side atoms of the ring and sets the derivative_map attribute of the monosaccharide
              * Entry examples for the derivative map: [-1, derivative pattern] [a, derivative pattern] [2, derivative pattern] ... [+1, derivative pattern] [+2, derivative pattern]
              * @param mono The monosaccharide object
              */
            void ExtractDerivatives(Glycan::Monosaccharide* mono);
            /*! \fn
              * A function in order to generate a complete name for the monosaccharide structure based on its derivatives
              * @param mono The monosaccharide object
              */
            void GenerateCompleteSugarName(Glycan::Monosaccharide* mono);
            /*! \fn
              * A function in order to add the modification info to the name of the sugar based on the first group of rules (No Bracket -> 2(r:6&!-1), Warning position -> a, Error Position -> 5(r6),4(r5) )
              * @param key The index of the atom which is a key in the derivative map
              * @param pattern The pattern of the modification which is based on the value in the derivative map
              * @param mono The monosaccharide object
              * @param long_name_pattern The long name pattern of the modification pattern which should be added to the sugar's long name
              * @param cond_name_pattern The short name pattern of the modification pattern which should be added to the sugar's short name
              * @param tail The stream that will be added to the end of monosaccharide name and will be updated by this function
              * @param head The stream that will be added to the beginning of monosaccharide name and will be updated by this function
              * @param minus_one A boolean value which indicates whether anomeric carbon is attached to a non-ring carbon
              * @param in_bracket The boolean value which indicates whether anomeric carbon is attached to a non-ring carbon
              */
            void AddModificationRuleOneInfo(std::string key, std::string pattern, Glycan::Monosaccharide* mono, std::string long_name_pattern, std::string cond_name_pattern,
                                            std::stringstream& head, std::stringstream& tail, bool minus_one, std::stringstream& in_bracket);
            /*! \fn
              * A function in order to add the modification info to the name of the sugar based on the first group of rules (No Bracket -> 2(r:6&!-1), Warning position -> a, Error Position -> 5(r6),4(r5) )
              * @param key The index of the atom which is a key in the derivative map
              * @param pattern The pattern of the modification which is based on the value in the derivative map
              * @param mono The monosaccharide object
              * @param long_name_pattern_at_minus_one The long name modification pattern for non-ring carbon attached to anomeric carbon which should be added to the sugar's long name
              * @param long_name_pattern_at_minus_one The long name modification pattern for non-ring carbon attached to last ring carbon which should be added to the sugar's long name
              * @param cond_name_pattern The short name pattern of the modification pattern which should be added to the sugar's short name
              * @param tail The stream that will be added to the end of monosaccharide name and will be updated by this function
              * @param minus_one A boolean value which indicates whether anomeric carbon is attached to a non-ring carbon
              * @param in_bracket The boolean value which indicates whether anomeric carbon is attached to a non-ring carbon
              */
            void AddModificationRuleTwoInfo(std::string key, std::string pattern, Glycan::Monosaccharide* mono, std::string long_name_pattern_at_minus_one, std::string long_name_pattern_at_plus_one,
                                            std::string cond_name_pattern, std::stringstream& tail, bool minus_one, std::stringstream& in_bracket);
            /*! \fn
              * A function in order to add the derivative info to the name of the sugar based on the first group of rules (Error Position -> 5(r6),4(r5) )
              * @param key The index of the atom which is a key in the derivative map
              * @param pattern The pattern of the modification which is based on the value in the derivative map
              * @param mono The monosaccharide object
              * @param long_name_pattern The long name pattern of the modification pattern which should be added to the sugar's long name
              * @param cond_name_pattern The short name pattern of the modification pattern which should be added to the sugar's short name
              * @param head The stream that will be added to the beginning of monosaccharide name and will be updated by this function
              * @param minus_one A boolean value which indicates whether anomeric carbon is attached to a non-ring carbon
              * @param in_bracket The boolean value which indicates whether anomeric carbon is attached to a non-ring carbon
              */
            void AddDerivativeRuleInfo(std::string key, std::string pattern, Glycan::Monosaccharide* mono, std::string long_name_pattern, std::string cond_name_pattern,
                                            std::stringstream& head, bool minus_one, std::stringstream& in_bracket);
            /*! \fn
              * A function in order to update the chemical code structure of a complex monosaccharide (monosaccharide with side atoms at position +2 and +3)
              * @param mono The monosaccharide object
              */
            void UpdateComplexSugarChemicalCode(Glycan::Monosaccharide* mono);
            /*! \fn
              * A function in order to extract oligosacchride structure based on the linkages between monosacchrides
              * @param monos The list of extracted monosaccharide object
              * @param dataset_residue_names List of known terminal residue names
              * @return oligosacchrides The list of extracted oligosacchrides
              */
            OligosaccharideVector ExtractOligosaccharides(std::vector<Glycan::Monosaccharide*> monos, gmml::ResidueNameMap dataset_residue_names);
            /*! \fn
              * A function in order to check if the target atom is attached to OME terminal
              * @param target_atom The atom which will be checked for a terminal
              * @return pattern The discovered pattern attached to the target atom
              */
            std::string CheckOMETerminal(Atom* target, AtomVector& terminal_atoms);
            /*! \fn
              * A function in order to check if the target atom is attached to ROH terminal
              * @param target_atom The atom which will be checked for a terminal
              * @return pattern The discovered pattern attached to the target atom
              */
            std::string CheckROHTerminal(Atom* target, AtomVector& terminal_atoms);
            /*! \fn
              * A function in order to check if the target atom is attached to TBT terminal
              * @param target_atom The atom which will be checked for a terminal
              * @return pattern The discovered pattern attached to the target atom
              */
            std::string CheckTBTTerminal(Atom* target, AtomVector& terminal_atoms);
            /*! \fn
              * A function in order to check if the target atom is attached to TBT, ROH or OME patterns when the residue name is different from the terminal residue names
              * @param target_atom The atom which will be checked for a terminal
              * @return pattern The discovered pattern attached to the target atom
              */
            std::string CheckTerminals(Atom* target, AtomVector& terminal_atoms);
            void CheckLinkageNote(Glycan::Monosaccharide* mono1, Glycan::Monosaccharide* mono2, std::string linkage, std::vector<std::string>& checked_linkages);

            /*! \fn
              * A function in order to calculate the orientation(S/R) of additional side atoms
              * @param prev_atom The previous neighbor atom of the current atom
              * @param target_atom The atom which the orientation of its oxygen neighbor is going to be calculated
              * @param next_atom The next neighbor atom of the current atom
              * @return orientation The relative orientation of the additional side atom's oxygen neighbor
              */
            std::string CalculateRSOrientations(Atom* prev_atom, Atom* target, Atom* next_atom);
            /*! \fn
              * A recursive function in order to build the tree like structure of the oligosacchrides (a directed graph made from the molecular graph structure of monosaccharides.
                each monosaccharide is a node and their linkages are the edges)
              * @param key The root monosacchride(node) of each oligosacchride
              * @param value The monosaccharide attached to the key monosacharide
              * @param oligo The oligosaccharide object to be constructed by the function
              * @param visited_monos The monosacchride objects that have been visited by previous calls to the function
              * @param monos_table A mapping between a monosacchride object and the list of monosacchrides attached to it
              * @param monos_table_linkages A mapping between a monosacchride object and the list of atoms that are involved in the linkages with the monosacchrides attached to
                the current monosacchride object
              * @param visited_linkages The list of linkages(atoms involved in the linkages between monosacchrides) that have been visited so far by calls to the function
              */
            void BuildOligosaccharideTreeStructure(Glycan::Monosaccharide* key, std::vector<Glycan::Monosaccharide*> val, Glycan::Oligosaccharide* oligo,
                                                                  std::vector<int>& visited_monos, std::map<Glycan::Monosaccharide*, std::vector<Glycan::Monosaccharide*> > monos_table,
                                                                  std::map<Glycan::Monosaccharide*, std::vector<std::string> > monos_table_linkages, std::vector<std::string>& visited_linkages);
            /*! \fn
              * A function in order to check if the target atom is attached to a derivative with the pattern xCH-N
              * @param target_atom The atom which will be checked for a derivative
              * @param cycle_atom_str The string version of atom identifiers of the cycle
              * @param pattern_atoms The list of atoms that involved in the pattern
              * @return pattern The discovered pattern of the attached derivative
              */
            std::string CheckxC_N(Atom* target, std::string cycle_atoms_str, AtomVector& pattern_atoms);
            /*! \fn
              * A function in order to check if the target atom is attached to a derivative with the pattern xC-O-C=OCH3 or xC-N-C=OCH3
              * @param target_atom The atom which will be checked for a derivative
              * @param cycle_atom_str The string version of atom identifiers of the cycle
              * @param NxO The option for checking the pattern with oxygen or nitrogen
              * @param pattern_atoms The list of atoms that involved in the pattern
              * @return pattern The discovered pattern of the attached derivative
              */
            std::string CheckxC_NxO_CO_C(Atom* target, std::string cycle_atoms_str, char NxO, AtomVector& pattern_atoms);
            /*! \fn
              * A function in order to check if the target atom is attached to a derivative with the pattern xC-N-C=OCH2OH or xC-O-C=OCH2OH
              * @param target_atom The atom which will be checked for a derivative
              * @param cycle_atom_str The string version of atom identifiers of the cycle
              * @param NxO The option for checking the pattern with oxygen or nitrogen
              * @param pattern_atoms The list of atoms that involved in the pattern
              * @return pattern The discovered pattern of the attached derivative
              */
            std::string CheckxC_NxO_CO_CO(Atom* target, std::string cycle_atoms_str, char NxO, AtomVector& pattern_atoms);
            /*! \fn
              * A function in order to check if the target atom is attached to a derivative with the pattern xC-N-SO3 or xC-O-SO3
              * @param target_atom The atom which will be checked for a derivative
              * @param cycle_atom_str The string version of atom identifiers of the cycle
              * @param NxO The option for checking the pattern with oxygen or nitrogen
              * @param pattern_atoms The list of atoms that involved in the pattern
              * @return pattern The discovered pattern of the attached derivative
              */
            std::string CheckxC_NxO_SO3(Atom* target, std::string cycle_atoms_str, char NxO, AtomVector& pattern_atoms);
            /*! \fn
              * A function in order to check if the target atom is attached to a derivative with the pattern xC-N-PO3 or xC-PO3
              * @param target_atom The atom which will be checked for a derivative
              * @param cycle_atom_str The string version of atom identifiers of the cycle
              * @param NxO The option for checking the pattern with oxygen or nitrogen
              * @param pattern_atoms The list of atoms that involved in the pattern
              * @return pattern The discovered pattern of the attached derivative
              */
            std::string CheckxC_NxO_PO3(Atom* target, std::string cycle_atoms_str, char NxO, AtomVector& pattern_atoms);
            /*! \fn
              * A function in order to check if the target atom is attached to a derivative with the pattern xC-N-CH3 or xC-O-CH3
              * @param target_atom The atom which will be checked for a derivative
              * @param cycle_atom_str The string version of atom identifiers of the cycle
              * @param NxO The option for checking the pattern with oxygen or nitrogen
              * @param pattern_atoms The list of atoms that involved in the pattern
              * @return pattern The discovered pattern of the attached derivative
              */
            std::string CheckxC_NxO_C(Atom* target, std::string cycle_atoms_str, char NxO, AtomVector& pattern_atoms);
            /*! \fn
              * A function in order to check if the target atom is attached to a derivative with the pattern C-(O,O), C-(O.OH), rC-(O.O) and rC-(O.OH)
              * @param target_atom The atom which will be checked for a derivative
              * @param cycle_atom_str The string version of atom identifiers of the cycle
              * @param NxO The option for checking the pattern with oxygen or nitrogen
              * @param pattern_atoms The list of atoms that involved in the pattern
              * @return pattern The discovered pattern of the attached derivative
              */
            std::string CheckxCOO(Atom* target, std::string cycle_atoms_str, AtomVector& pattern_atoms);

            void Ionizing(std::string ion_name, std::string lib_file, std::string parameter_file, int ion_count = 0);
            void Solvation(double extension, double closeness, std::string lib_file);

            double GetTotalCharge();
            double GetRadius();
            double GetTotalMass();
            void GetCenterOfMass(GeometryTopology::Coordinate* center_of_mass);
            void GetCenterOfGeometry(GeometryTopology::Coordinate* center_of_geometry);
            void GetBoundary(GeometryTopology::Coordinate* lower_left_back_corner, GeometryTopology::Coordinate* upper_right_front_corner);
            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the assembly contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cout);

            void PrettyPrintHet(std::ostream& out = std::cout);
            void PrintHetResidues(std::ostream& out = std::cout);
            void PrintHetAtoms(std::ostream& out = std::cout);

            void WriteHetResidues(std::string file_name);
            void WriteHetAtoms(std::string file_name);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            std::string name_;                              /*!< Assembly name >*/
            AssemblyVector assemblies_;                     /*!< List of assemblies that may involved in a big structure as a new assembly >*/
            ResidueVector residues_;                        /*!< List of residues involved in the current object of assembly >*/
            std::string chemical_type_;                     /*!< A descriptor for the chemical type of the current object of assembly >*/
            int sequence_number_;                           /*!< An integer number to indicates sequence of importing the current assembly >*/
            std::string id_;
            std::string description_;                       /*!< Short description for the current assembly >*/
            std::string source_file_;                       /*!< File name that the current assembly has been built upon >*/
            gmml::InputFileType source_file_type_;          /*!< Type of the file that the current assembly has been built upon >*/
            int model_index_;                               /*!< In case that there are more than one models for an assembly, this attribute indicated which model is the target model >*/
            NoteVector notes_;                              /*!< A list of note instances from the Note struct in Glycan name space which is used for representing the potential issues within a structure >*/
    };

    struct DistanceCalculationThreadArgument{
            int thread_index;
            int number_of_threads;
            int model_index;
            double cutoff;
            Assembly* a;
            DistanceCalculationThreadArgument()
            {
                thread_index = 0;
                number_of_threads = 1;
                model_index = 0;
                cutoff = gmml::dCutOff;
                a = NULL;
            }

            DistanceCalculationThreadArgument(int ti, int tn, int mi, double c, Assembly* assembly)
            {
                thread_index = ti;
                number_of_threads = tn;
                model_index = mi;
                cutoff = c;
                a = assembly;
            }
    };

    struct DistanceCalculationByMatrixThreadArgument{
            int thread_index;
            int model_index;
            double cutoff;
            std::vector<Atom*>* first_chunk;
            std::vector<Atom*>* second_chunk;

            DistanceCalculationByMatrixThreadArgument()
            {
                thread_index = 0;
                model_index = 0;
                cutoff = gmml::dCutOff;
                first_chunk = new std::vector<Atom*>();
                second_chunk = new std::vector<Atom*>();
            }

            DistanceCalculationByMatrixThreadArgument(int ti, int mi, double c, std::vector<Atom*>* fc, std::vector<Atom*>* sc)
            {
                thread_index = ti;
                model_index = mi;
                cutoff = c;
                first_chunk = fc;
                second_chunk = sc;
            }
    };


}

#endif // ASSEMBLY_HPP
