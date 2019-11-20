#ifndef LIBRARYFILE_HPP
#define LIBRARYFILE_HPP

#include <string>
#include <map>
#include <iostream>
#include <vector>
#include "../../common.hpp"

namespace LibraryFileSpace
{
    class LibraryFileResidue;
    class LibraryFileAtom;
    class LibraryFile
    {
        public:
            //////////////////////////////////////////////////////////
            //                     TYPE DEFINITION                  //
            //////////////////////////////////////////////////////////
            /*! \def
              * A mapping between a name of a residue and its residue object
              */
            typedef std::map<std::string, LibraryFileResidue*> ResidueMap;

            //////////////////////////////////////////////////////////
            //                       Constructor                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            LibraryFile();
            /*! \fn
              * Constructor
              * @param lib_file An existing library file path to be read
              */
            LibraryFile(const std::string& lib_file);

            ~LibraryFile();

            //////////////////////////////////////////////////////////
            //                           ACCESSOR                   //
            //////////////////////////////////////////////////////////
/** \addtogroup Molecular_Data_Structure
               * @{
               */
            /*! \fn
              * An accessor function in order to access to library file path of the current object
              * @return path_ attribute of the current object of this class
              */
            const std::string& GetFilePath() const;
            /*! \fn
              * An accessor function in order to access to residue map created based on the content of the given file
              * @return residues_ attribute of the current object of this class
              */
            const ResidueMap& GetResidues() const;
            /*! \fn
              * An accessor function in order to access to all residue names of the current object
              * @return residue_names residue names of the current object of this class
              */
            std::vector<std::string> GetAllResidueNames();

            //*****************************************************
            gmml::ResidueNameMap GetAllResidueNamesMap();

            //*****************************************************

            /*! \fn
              * An accessor function in order to access to all atom names of the current object
              * @param residue_name the residue name of the current object
              * @return atom_names_of_residue The atom names of the current object of this class
              */
            std::vector<std::string> GetAllAtomNamesOfResidue(std::string residue_name);
/**@}*/
            //////////////////////////////////////////////////////////
            //                           MUTATOR                    //
            //////////////////////////////////////////////////////////
/** \addtogroup Manipulators
               * @{
               */
            /*! \fn
              * A mutator function in order to set the path of the current library file
              * Set the path_ attribute of the current library file
              * @param path The path attribute of the current object
              */
            void SetPath(std::string path);
            /*! \fn
              * A mutator function in order to set the residues of the current library file
              * Set the residues_ attribute of the current library file
              * @param residues The residue map attribute of the current object
              */
            void SetResidues(ResidueMap residues);
/** @}*/
            //////////////////////////////////////////////////////////
            //                         FUNCTIONS                    //
            //////////////////////////////////////////////////////////
/** \addtogroup Input_File_Reader
               * @{
               */
            /*! \fn
              * A function to parse the contents of a given stream of a file
              * Parse the given stream and set the attributes of the current object accordingly
              * @param in_file A stream contains whole contents of a library file
              */
            void Read(std::ifstream& in_file);
            /*! \fn
              * Process a line of the atom section of a library file and create a new atom object
              * @param line A line from the atom section of a library file
              * @return A new LibraryFileAtom instance created by the information of the given line
              */
/** @}*/
            /** \addtogroup Molecular_Data_Structure
               * @{
               */
            LibraryFileAtom* ProcessAtom(std::string& line);
            /*! \fn
              * A function in order to access to library file residue by a residue name
              * @param residue_name The name of the residue
              * @return library_file_residue
              */
/** @}*/
            /** \addtogroup Output_File_Builder
               * @{
               */
            LibraryFileResidue* GetLibraryResidueByResidueName(std::string residue_name);
            /*! \fn
              * A function to write back a library file into an output file
              * @param library_file Library output file name
              */
/** @}*/
            void Write(const std::string& library_file);
            /*! \fn
              * A function to write a library file into an output stream
              * @param out_stream Output stream
              */
/** \addtogroup Output_File_Builder
               * @{
               */
            void BuildLibraryFile(std::ofstream& out_stream);
            /*! \fn
              * A function in order to write the atom section of a specified residue into an output stream
              * @param stream Output stream
              * @param residue Specified residue in the library file that has to be written into an output stream
              */
/** @}*/
             /** \addtogroup Verifiers_and_Issue_Resolvers
               * @{
               */
            void ResolveAtomSection(std::ofstream& stream, LibraryFileResidue* residue);
            /*! \fn
              * A function in order to write the atom pert info section of a specified residue into an output stream
              * @param stream Output stream
              * @param residue Specified residue in the library file that has to be written into an output stream
              */
            void ResolveAtomPertInfoSection(std::ofstream& stream, LibraryFileResidue* residue);
            /*! \fn
              * A function in order to write the bound box section of a specified residue into an output stream
              * @param stream Output stream
              * @param residue Specified residue in the library file that has to be written into an output stream
              */
            void ResolveBoundBoxSection(std::ofstream& stream, LibraryFileResidue* residue);
            /*! \fn
              * A function in order to write the child sequence section of a specified residue into an output stream
              * @param stream Output stream
              * @param residue Specified residue in the library file that has to be written into an output stream
              */
            void ResolveChildSequenceSection(std::ofstream& stream, LibraryFileResidue* residue);
            /*! \fn
              * A function in order to write the connect section of a specified residue into an output stream
              * @param stream Output stream
              * @param residue Specified residue in the library file that has to be written into an output stream
              */
            void ResolveConnectSection(std::ofstream& stream, LibraryFileResidue* residue);
            /*! \fn
              * A function in order to write the connectivity section of a specified residue into an output stream
              * @param stream Output stream
              * @param residue Specified residue in the library file that has to be written into an output stream
              */
            void ResolveConnectivitySection(std::ofstream& stream, LibraryFileResidue* residue);
            /*! \fn
              * A function in order to write the hierarchy section of a specified residue into an output stream
              * @param stream Output stream
              * @param residue Specified residue in the library file that has to be written into an output stream
              */
            void ResolveHierarchySection(std::ofstream& stream, LibraryFileResidue* residue);
            /*! \fn
              * A function in order to write the name section of a specified residue into an output stream
              * @param stream Output stream
              * @param residue Specified residue in the library file that has to be written into an output stream
              */
            void ResolveNameSection(std::ofstream& stream, LibraryFileResidue* residue);
            /*! \fn
              * A function in order to write the position section of a specified residue into an output stream
              * @param stream Output stream
              * @param residue Specified residue in the library file that has to be written into an output stream
              */
            void ResolvePositionSection(std::ofstream& stream, LibraryFileResidue* residue);
            /*! \fn
              * A function in order to write the residue connect section of a specified residue into an output stream
              * @param stream Output stream
              * @param residue Specified residue in the library file that has to be written into an output stream
              */
            void ResolveResidueConnectSection(std::ofstream& stream, LibraryFileResidue* residue);
            /*! \fn
              * A function in order to write the residues section of a specified residue into an output stream
              * @param stream Output stream
              * @param residue Specified residue in the library file that has to be written into an output stream
              */
            void ResolveResiduesSection(std::ofstream& stream, LibraryFileResidue* residue);
            /*! \fn
              * A function in order to write the residue pdb sequence number section of a specified residue into an output stream
              * @param stream Output stream
              * @param residue Specified residue in the library file that has to be written into an output stream
              */
            void ResolveResiduePdbSequenceNumberSection(std::ofstream& stream, LibraryFileResidue* residue);
            /*! \fn
              * A function in order to write the solvent cap section of a specified residue into an output stream
              * @param stream Output stream
              * @param residue Specified residue in the library file that has to be written into an output stream
              */
            void ResolveSolventCapSection(std::ofstream& stream, LibraryFileResidue* residue);
            /*! \fn
              * A function in order to write the velocities section of a specified residue into an output stream
              * @param stream Output stream
              * @param residue Specified residue in the library file that has to be written into an output stream
              */
            void ResolveVelocitiesSection(std::ofstream& stream, LibraryFileResidue* residue);
/**@}*/
            //////////////////////////////////////////////////////////
            //                     DISPLAY FUNCTIONS                //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the file contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cerr);

        private:
            //////////////////////////////////////////////////////////
            //                         ATTRIBUTES                   //
            //////////////////////////////////////////////////////////
            std::string path_;          /*!< Path of the given library file */
            ResidueMap residues_;       /*!< Map of residues included in the given file mapped to their names */
    };
}

#endif // LIBRARYFILE_HPP
