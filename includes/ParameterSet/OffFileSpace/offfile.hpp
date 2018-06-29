#ifndef OFFFILE_HPP
#define OFFFILE_HPP

#include <string>
#include <map>
#include <iostream>
#include <vector>
#include "../../common.hpp"
#include "../../MolecularModeling/assembly.hpp"
#include "../../MolecularModeling/residue.hpp"
#include "../../MolecularModeling/atom.hpp"
#include "../../MolecularModeling/atomnode.hpp"

namespace OffFileSpace
{
    class OffFile
    {
        public:
            //////////////////////////////////////////////////////////
            //                     TYPE DEFINITION                  //
            //////////////////////////////////////////////////////////

            typedef std::vector<MolecularModeling::Residue*> ResidueVector;
            typedef std::vector<MolecularModeling::Atom*> AtomVector;

            //////////////////////////////////////////////////////////
            //                       Constructor                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            OffFile();

            /*! \fn
              * Default deconstructor
              */

            ~OffFile();

            //////////////////////////////////////////////////////////
            //                           ACCESSOR                   //
            //////////////////////////////////////////////////////////


            //////////////////////////////////////////////////////////
            //                           MUTATOR                    //
            //////////////////////////////////////////////////////////


            //////////////////////////////////////////////////////////
            //                         FUNCTIONS                    //
            //////////////////////////////////////////////////////////

            /*! \fn
              * A function to write an assembly in Off file format
              * @param out_stream Output stream
              */
            void Write(std::string file_name, int CoordinateIndex, MolecularModeling::Assembly* assembly);

            /*! \fn
              * A function in order to write the atom section of a specified residue into an output stream
              * @param stream Output stream
              */
            void WriteAtomSection(std::ofstream& stream, ResidueVector assembly_residues);

            /*! \fn
              * A function in order to write the atom pert info section of a specified residue into an output stream
              * @param stream Output stream
              */
            void WriteAtomPertInfoSection(std::ofstream& stream, ResidueVector assembly_residues);

            /*! \fn
              * A function in order to write the bound box section of a specified residue into an output stream
              * @param stream Output stream
              */
            void WriteBoundBoxSection(std::ofstream& stream, MolecularModeling::Assembly* assembly);

            /*! \fn
              * A function in order to write the child sequence section of a specified residue into an output stream
              * @param stream Output stream
              */
            void WriteChildSequenceSection(std::ofstream& stream, ResidueVector assembly_residues);

            /*! \fn
              * A function in order to write the connect section of a specified residue into an output stream
              * @param stream Output stream
              */
            void WriteConnectSection(std::ofstream& stream, ResidueVector assembly_residues);

            /*! \fn
              * A function in order to write the connectivity section of a specified residue into an output stream
              * @param stream Output stream
              */
            void WriteConnectivitySection(std::ofstream& stream, ResidueVector assembly_residues);

            /*! \fn
              * A function in order to write the hierarchy section of a specified residue into an output stream
              * @param stream Output stream
              */
            void WriteHierarchySection(std::ofstream& stream, ResidueVector assembly_residues);

            /*! \fn
              * A function in order to write the name section of a specified residue into an output stream
              * @param stream Output stream
              */
            void WriteNameSection(std::ofstream& stream);

            /*! \fn
              * A function in order to write the position section of a specified residue into an output stream
              * @param stream Output stream
              */
            void WritePositionSection(std::ofstream& stream, ResidueVector assembly_residues, int CoordinateIndex);

            /*! \fn
              * A function in order to write the residue connect section of a specified residue into an output stream
              * @param stream Output stream
              */
            void WriteResidueConnectSection(std::ofstream& stream, ResidueVector assembly_residues);

            /*! \fn
              * A function in order to write the residues section of a specified residue into an output stream
              * @param stream Output stream
              */
            void WriteResiduesSection(std::ofstream& stream, ResidueVector assembly_residues);

            /*! \fn
              * A function in order to write the solvent cap section of a specified residue into an output stream
              * @param stream Output stream
              */
            void WriteSolventCapSection(std::ofstream& stream);

            /*! \fn
              * A function in order to write the velocities section of a specified residue into an output stream
              * @param stream Output stream
              */
            void WriteVelocitiesSection(std::ofstream& stream, ResidueVector assembly_residues);


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

             std::string unit_name_;                  /*!< Name of OFF file unit >*/

    };
}

#endif // OFFFILE_HPP
