#ifndef ASSEMBLY_HPP
#define ASSEMBLY_HPP

#include <string>
#include <iostream>
#include <vector>

#include "../Geometry/coordinate.hpp"
#include "../common.hpp"

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

            Assembly(std::vector<std::string> file_paths, gmml::InputFileType type);

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
            void BuildAssemblyFromPdbFile(std::string pdb_file_path);
            void BuildAssemblyFromTopologyFile(std::string topology_file_path);
            void BuildAssemblyFromLibraryFile(std::string library_file_path);
            void BuildAssemblyFromTopologyCoordinateFile(std::string topology_file_path, std::string coordinate_file_path);
            void BuildAssemblyFromPrepFile(std::string prep_file_path);

            void BuildStructure(gmml::BuildingStructureOption building_option, std::vector<std::string> options, std::vector<std::string> file_paths);
            void BuildStructureByDistance(double cutoff = gmml::dCutOff, int model_index = 0);
            void BuildStructureByOriginalFileBondingInformation();
            void BuildStructureByPDBFileInformation();
            void BuildStructureByTOPFileInformation();
            void BuildStructureByLIBFileInformation();
            void BuildStructureByPrepFileInformation();
            void BuildStructureByDatabaseFilesBondingInformation(std::vector<gmml::InputFileType> types, std::vector<std::string> file_paths);
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
            std::string name_;
            AssemblyVector assemblies_;
            ResidueVector residues_;
            std::string chemical_type_;
            int sequence_number_;
            double total_mass_;
            Geometry::Coordinate center_of_mass_;
            Geometry::Coordinate center_of_geometry_;
            std::string description_;
            std::string source_file_;
            gmml::InputFileType source_file_type_;

    };
}

#endif // ASSEMBLY_HPP
