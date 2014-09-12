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
#include "../PrepFileSpace/prepfile.hpp"
#include "../LibraryFileSpace/libraryfile.hpp"

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
              * @param source_file_type The residues attribute of the current object
              */
            void SetSourceFileType(gmml::InputFileType source_file_type);           

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
              * A functions that extracts all atoms of an assembly
              * @return Vector of atoms all in the current object of assembly
              */
            AtomVector GetAllAtomsOfAssembly();

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

    };
}

#endif // ASSEMBLY_HPP
