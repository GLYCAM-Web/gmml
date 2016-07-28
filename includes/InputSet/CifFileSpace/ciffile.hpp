#ifndef CIFFILE_HPP
#define CIFFILE_HPP

#include <string>
#include <vector>
#include <iostream>

namespace CifFileSpace{

    class CifFileAtom;
    class CifFile
    {
        public:

            //////////////////////////////////////////////////////////
            //                    TYPE DEFINITION                   //
            //////////////////////////////////////////////////////////
            /*! \typedef
              * List of cif atoms
              */
            typedef std::vector<CifFileAtom*> CifFileAtomVector;

            //////////////////////////////////////////////////////////
            //                       Constructor                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            CifFile();

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
            /*! \fn
              * An accessor function in order to access to cif file path of the current object
              * @return path_ attribute of the current object of this class
              */
            std::string GetPath();
            /*! \fn
              * An accessor function in order to access to the residue name attribute of the current object
              * @return residue_name_ attribute of the current object of this class
              */
            std::string GetResidueName();
            /*! \fn
              * An accessor function in order to access to list of atoms of the the residue in the current object
              * @return atoms The vector of atoms of the residue in the file
              */
            CifFileAtomVector GetAtoms();

            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function in order to set the path of the cif file
              * @param cif_path Cif file path
              */
            void SetPath(std::string cif_path);
            /*! \fn
              * A mutator function in order to set the residue name of the current object
              * Set the resiude_name_ attribute of the current cif file
              * @param residue_name The residue name of the current object
              */
            void SetResidueName(const std::string residue_name);
            /*! \fn
              * A mutator function in order to set atoms of the current object
              * Set the atoms_ attribute of the current cif file
              * @param atoms The list of atoms of the current object
              */
            void SetAtoms(CifFileAtomVector atoms);

            //////////////////////////////////////////////////////////
            //                        FUNCTIONS                     //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to parse the contents of a given stream of a file
              * Parse the given stream and set the attributes of the current object accordingly
              * @param in_file A stream contains whole contents of a cif file
              */
            void Read(std::ifstream& in_file);

            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the cif file contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cout);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            std::string path_;              /*!< The path of the current cif file */
            std::string residue_name_;       /*!< The name of the residue in the current cif file */
            CifFileAtomVector atoms_;        /*!< The list of atoms of the residue in the current cif file */
    };
}

#endif // CIFFILE_HPP
