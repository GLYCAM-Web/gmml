#ifndef PARAMETERFILE_HPP
#define PARAMETERFILE_HPP

#include <string>
#include <vector>
#include <map>
#include <iostream>
#include "../../common.hpp"

namespace ParameterFileSpace
{
    ///////////////////// FORWARD DECLARATION //////////////////
    class ParameterFileAtom;
    class ParameterFileBond;
    class ParameterFileAngle;
    class ParameterFileDihedral;
    class ParameterFile
    {
        public:
            //////////////////////////////////////////////////////////
            //                     TYPE DEFINITION                  //
            //////////////////////////////////////////////////////////
            /*! \def
              * A mapping between an atom type and its atom object
              */
            typedef std::map<std::string, ParameterFileAtom*> AtomTypeMap;
            /*! \def
              * A mapping between a connection of two atom types and its bond object
              */
            typedef std::map<std::vector<std::string>, ParameterFileBond*> BondMap;
            /*! /def
              * A mapping between a connection of three atom types and its angle object
              */
            typedef std::map<std::vector<std::string>, ParameterFileAngle*> AngleMap;
            /*! \def
              * A mapping between a connection of four atom types and its dihedral object
              */
            typedef std::map<std::vector<std::string>, ParameterFileDihedral*> DihedralMap;

            //////////////////////////////////////////////////////////
            //                       Constructor                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Constructor
              * @param param_file An existing library file path to be read
              */
            ParameterFile(std::string param_file, int type = gmml::MAIN);

            //////////////////////////////////////////////////////////
            //                           ACCESSOR                   //
            //////////////////////////////////////////////////////////
            /*! \fn
              * An accessor function in order to access to parameter file path of the current object
              * @return path_ attribute of the current object of this class
              */
            const std::string& GetFilePath();
            /*! \fn
              * An accessor function in order to access to the title of the current object
              * @return title_ attribute of the current object of this class
              */
            const std::string& GetTitle();
            /*! \fn
              * An accessor function in order to access to the map of atom types (as string) to their atom objects of the current object
              * @return atom_types_ attribute of the current object of this class
              */
            const AtomTypeMap& GetAtomTypes();
            /*! \fn
              * An accessor function in order to access to the map of bonds (dual string of atom types) to their bond objects of the current object
              * @return bonds_ attribute of the current object of this class
              */
            const BondMap& GetBonds();
            /*! \fn
              * An accessor function in order to access to the map of angles (triple string of atom types) of the current object
              * @return angles_ attribute of the current object of this class
              */
            const AngleMap& GetAngles();
            /*! \fn
              * An accessor function in order to access to the map of dihedrals (fourple string of atom types) of the current object
              * @return path_ attribute of the current object of this class
              */
            const DihedralMap& GetDihedrals();
            /*! \fn
              * An accessor function in order to access to the parameter file type of the current object
              * @return file_type_ attribute of the current object of this class
              */
            const int GetParameterFileType();
            /*! \fn
              * An accessor function in order to access all improper dihedrals of the current object
              * @return improper_dihedrals Improper dihedrals the current object of this class
              */
            DihedralMap GetAllImproperDihedrals();
            /*! \fn
              * An accessor function in order to access all proper dihedrals of the current object
              * @return proper_dihedrals proper dihedrals the current object of this class
              */
            DihedralMap GetAllproperDihedrals();

            //////////////////////////////////////////////////////////
            //                           MUTATOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A mutator function in order to set the file type of the current object
              * Set the file_type_ attribute of the current parameter file
              * @param file_type The file type attribute of the current object
              */
            void SetParameterFileType(int file_type);
            //////////////////////////////////////////////////////////
            //                         FUNCTIONS                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to parse the contents of a given stream of a file
              * Parse the given stream and set the attributes of the current object accordingly
              * @param in_file A stream contains whole contents of a parameter file
              */
            void ReadMainParameter(std::ifstream& in_file);
            /*! \fn
              * A function to parse the contents of a given stream of a file
              * Parse the given stream and set the attributes of the current object accordingly
              * @param in_file A stream contains whole contents of a frcmod parameter file
              */
            void ReadModifiedParameter(std::ifstream& in_file);
            /*! \fn
              * A function to parse the contents of a given stream of a file
              * Parse the given stream and set the attributes of the current object accordingly
              * @param in_file A stream contains whole contents of a frcmod parameter file containing ions
              */
            void ReadIonicModifiedParameter(std::ifstream& in_file);
            /*! \fn
              * A function that parses a line of atom type section of the current object
              * Process the atom type lines of the parameter file
              * @param line A line of the atom type section of the current parameter file
              */
            void ProcessAtomType(const std::string& line);
            /*! \fn
              * A function to paras a line of hydrophicil atoms section of the current object
              * Process the hydrophilic atom type line of the parameter file
              * @param line A line of the hydrophilic atom types section of the current parameter file
              */
            void ProcessHydrophilicAtomType(const std::string& line);
            /*! \fn
              * A function that parses a line of bond section of the current object
              * Process the bond lines of the parameter file
              * @param line A line of the bond section of the current parameter file
              */
            void ProcessBond(const std::string& line);
            /*! \fn
              * A function to parse a line corresponding to the angle section of the current object
              * Process the angle lines of the parameter file
              * @param line A line of the angle section of the current parameter file
              */
            void ProcessAngle(const std::string& line);
            /*! \fn
              * A function that parses a line of dihedral section of the current object
              * Process the dihedral lines of the parameter file
              * @param line A line of dihedral section of the current parameter file
              * @param line_number The line number of the current read line from the parameter file
              * @param in_file A stream of the current parameter file
              */
            void ProcessDihedral(std::string& line, int& line_number, std::ifstream& in_file);
            /*! \fn
              * A function that parses a line of improper dihedral section of the current object
              * Process the improper dihedral lines of the parameter file
              * @param line A line of improper dihedral section of the current parameter file
              * @param line_number The line number of the current read line from the parameter file
              * @param in_file A stream of the current parameter file
              */
            void ProcessImproperDihedral(std::string& line, int& line_number, std::ifstream& in_file);
            /*! \fn
              * A function that parses a line of hydrogen bond section of the current object
              * Process the hydrogen bond lines of the parameter file
              * @param line A line of hydrogen bond section of the current parameter file
              */
            void ProcessHydrogenBond(const std::string& line);
            /*! \fn
              * A function that parses a line of equivalent symbols section of the current object
              * Process the equivalent symbols lines of the parameter file
              * @param line A line of equivalent symbols section of the current parameter file
              */
            void ProcessEquivalentSymbols(const std::string& line);
            /*! \fn
              * A function that parses a line of potential parameter section of the current object
              * Process the potential parameter lines of the parameter file located at the end of the file
              * @param line A line of potential parameter section of the current parameter file
              */
            void ProcessPotentialParameter(const std::string& line);
            /*! \fn
              * A function that parses a description of a line in dihedral section of the current object to extract double value of the given key
              * Process the description of a line in dihedral section of the parameter file
              * @param dscr The description part of a line in dihedral section of the current parameter file
              * @param key A key string which a value is assigned to it in the first parameter
              * @return A double value of the given key string
              */
            double ProcessDoubleDihedralDescription(const std::string& dscr, const std::string& key);
            /*! \fn
              * A function in order to write back a parameter file into an output file
              * @param parameter_file Parameter output file name
              */
            void Write(const std::string& parameter_file);
            /*! \fn
              * A function that writes back a main parameter file (mostly indicates as dat file) content into an output stream
              * @param out_stream Main parameter file output stream
              */
            void BuildMainParameterFile(std::ofstream& out_stream);
            /*! \fn
              * A function that writes back a modified parameter file (mostly indicates as frcmod file) content into an output stream
              * @param out_stream Modified parameter file output stream
              */
            void BuildModifiedParameterFile(std::ofstream& out_stream);
            /*! \fn
              * A function that writes back a modified parameter file containing ions (mostly indicates as frcmod file) content into an output stream
              * @param out_stream Modified parameter file output stream
              */
            void BuildIonicModifiedParameterFile(std::ofstream& out_stream);
            /*! \fn
              * A function that writes back atom type section of a parameter file into an output stream
              * @param stream Parameter file output stream
              */
            void ResolveAtomTypeSection(std::ofstream& stream);
            /*! \fn
              * A function that writes back hydrophilic atom type section of a parameter file into an output stream
              * @param stream Parameter file output stream
              */
            void ResolveHydrophilicAtomTypeSection(std::ofstream& stream);
            /*! \fn
              * A function that writes back bond section of a parameter file into an output stream
              * @param stream Parameter file output stream
              */
            void ResolveBondSection(std::ofstream& stream);
            /*! \fn
              * A function that writes back angle section of a parameter file into an output stream
              * @param stream Parameter file output stream
              */
            void ResolveAngleSection(std::ofstream& stream);
            /*! \fn
              * A function that writes back dihedral section of a parameter file into an output stream
              * @param stream Parameter file output stream
              */
            void ResolveDihedralSection(std::ofstream& stream);
            /*! \fn
              * A function that writes back improper dihedral section of a parameter file into an output stream
              * @param stream Parameter file output stream
              */
            void ResolveImproperDihedralSection(std::ofstream& stream);
            /*! \fn
              * A function that writes back hydrogen bond section of a parameter file into an output stream
              * @param stream Parameter file output stream
              */
            void ResolveHydrogenBondSection(std::ofstream& stream);
            /*! \fn
              * A function that writes back equivalent symbols section of a parameter file into an output stream
              * @param stream Parameter file output stream
              */
            void ResolveEquivalentSymbolsSection(std::ofstream& stream);
            /*! \fn
              * A function that writes back potential parameter section of a parameter file into an output stream
              * @param stream Parameter file output stream
              */
            void ResolvePotentialParameterSection(std::ofstream& stream);

            //////////////////////////////////////////////////////////
            //                     DISPLAY FUNCTIONS                //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the file contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cout);

        private:
            //////////////////////////////////////////////////////////
            //                         ATTRIBUTES                   //
            //////////////////////////////////////////////////////////
            std::string path_;           /*!< Path of the given parameter file */
            std::string title_;          /*!< Title of the parameter file */
            AtomTypeMap atom_types_;     /*!< A collection of mapping between atom type and its attributes*/
            BondMap bonds_;              /*!< A collection of mapping between bond (double atom types) and its attributes*/
            AngleMap angles_;            /*!< A collection of mapping between angle (tripple atom types) and its attributes*/
            DihedralMap dihedrals_;      /*!< A collection of mapping between dihedral (quad atom types) and its attributes*/
            int file_type_;              /*!< An integer number that indicates main parameter file and modified one from each other >*/

    };
}

#endif // PARAMETERFILE_HPP
