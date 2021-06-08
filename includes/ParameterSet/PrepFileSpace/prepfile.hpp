#ifndef PREPFILE_HPP
#define PREPFILE_HPP

#include <map>
#include <string>
#include <iostream>
#include <vector>
#include "../../common.hpp"

namespace PrepFileSpace
{
    class PrepFileResidue;
    class PrepFile
    {
        public:
            //////////////////////////////////////////////////////////
            //                     TYPE DEFINITION                  //
            //////////////////////////////////////////////////////////
            /*! \def
              * A mapping between a residue name and its residue object
              */
            typedef std::map< std::string, PrepFileResidue* > ResidueMap;

            //////////////////////////////////////////////////////////
            //                       Constructor                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Constructor
              * @param prep_file An existing prep file path to be read
              */
            PrepFile(const std::string& prep_file);
            PrepFile(const std::string& prep_file, std::vector<std::string>& query_residue_names);
            /*! \fn
              * Default constructor
              */
            PrepFile();
            ~PrepFile();
            //////////////////////////////////////////////////////////
            //                           ACCESSOR                   //
            //////////////////////////////////////////////////////////
/** \addtogroup Molecular_Data_Structure
               * @{
               */
            /*! \fn
              * An accessor function in order to access to prep file path of the current object
              * @return path_ attribute of the current object of this class
              */
            ResidueMap& GetResidues();
            /*! \fn
              * An accessor function in order to access to all residue names of the current object
              * @return residue_names residue names of the current object of this class
              */
            std::vector<std::string> GetAllResidueNames();

            //**************************************
            gmml::ResidueNameMap GetAllResidueNamesMap();

            //**************************************

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
              * A mutator function in order to set the path of the current prep file
              * Set the path_ attribute of the current library file
              * @param path The path attribute of the current object
              */
            void SetPath(std::string path);
            /*! \fn
              * A mutator function in order to set the residues of the current prep file
              * Set the residues_ attribute of the current prep file
              * @param residues The residue map attribute of the current object
              */
            void SetResidues(ResidueMap residues);
/**@}*/
            //////////////////////////////////////////////////////////
            //                         FUNCTIONS                    //
            //////////////////////////////////////////////////////////
/** \addtogroup Input_File_Reader
               * @{
               */
            /*! \fn
              * A function to parse the contents of a given stream of a file
              * Parse the given stream and set the attributes of the current object accordingly
              * @param in_file A stream contains whole contents of a prep file
              */
            void ReadAllResidues(std::ifstream& in_file);
            void ReadOnlyQueryResidues(std::ifstream &in_file, std::vector<std::string>& query_residue_names);
/** @}*/
            /** \addtogroup Output_File_Builder
               * @{
               */
            /*! \fn
              * A function to residue section of the current prep file
              * Parse the given stream and set the attributes of the current object accordingly
              * @param in_file A stream contains whole contents of a prep file
              */
            PrepFileResidue* ProcessResidueSection(std::ifstream& in_file);
            /*! \fn
              * A function in order to write back a prep file into an output file
              * @param prep_file Output prep file name
              */
            void Write(const std::string& prep_file);
            /*! \fn
              * A function in order to write sections of a prep file in an output stream in order for writing into an output file
              * @param out_stream Output stream
              */
            void BuildPrepFile(std::ofstream& out_stream);
/** @}*/
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
            std::string path_;              /*!< Actual path of the given prep file */
            ResidueMap residues_;           /*!< Fill by PrepFileResidues */
            /*!< End of a prep file gets marked by STOP */
    };
}

#endif // PREPFILE_HPP
