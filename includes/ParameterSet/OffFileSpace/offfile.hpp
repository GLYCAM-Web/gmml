#ifndef OFFFILE_HPP
#define OFFFILE_HPP

#include <string>
#include <map>
#include <iosfwd>
#include <vector>
#include "../../common.hpp"
//#include "../../gmml.hpp"
#include "../../MolecularModeling/assembly.hpp"
#include "../../MolecularModeling/residue.hpp"
#include "../../MolecularModeling/atom.hpp"
#include "../../MolecularModeling/atomnode.hpp"
#include "../../MolecularMetadata/elementattributes.hpp"

namespace OffFileSpace
{
    class OffFileResidue;
    class OffFileAtom;  
    class OffFile
    {
        public:
            //////////////////////////////////////////////////////////
            //                     TYPE DEFINITION                  //
            //////////////////////////////////////////////////////////

            //typedef std::vector<MolecularModeling::Residue*> ResidueVector;
            typedef std::vector<OffFileSpace::OffFileResidue*> OffFileResidueVector;
            typedef std::vector<OffFileSpace::OffFileAtom*> OffFileAtomVector;
            typedef std::map<int, int>AtomIndexMap;
            typedef std::map<int, int>AtomBondingMap;


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
              * A function to populate off file residues from the assembly
              * @param  assembly_residues The residues of the current assembly
              * @param CoordinateIndex The coodinate index of the coordinate vector
              * return off_file_residues_ The off file residues of the current off assembly
              */
            void PopulateOffFileResiduesFromAssembly(MolecularModeling::ResidueVector assembly_residues,int CoordinateIndex);
            /*! \fn
              * A function to write an assembly in Off file format
              * @param out_stream Output stream
              */
            void Write(std::string file_name, int CoordinateIndex, MolecularModeling::Assembly* assembly);

            /*! \fn
              * A function in order to write the atom section of a specified residue into an output stream
              * @param stream Output stream
              */
            void WriteAtomSection(std::ofstream& stream, OffFileResidueVector off_file_residues);

            /*! \fn
              * A function in order to write the atom pert info section of a specified residue into an output stream
              * @param stream Output stream
              */
            void WriteAtomPertInfoSection(std::ofstream& stream, OffFileResidueVector off_file_residues);

            /*! \fn
              * A function in order to write the bound box section of a specified residue into an output stream
              * @param stream Output stream
              */
            void WriteBoundBoxSection(std::ofstream& stream);

            /*! \fn
              * A function in order to write the child sequence section of a specified residue into an output stream
              * @param stream Output stream
              */
            void WriteChildSequenceSection(std::ofstream& stream,OffFileResidueVector off_file_residues);

            /*! \fn
              * A function in order to write the connect section of a specified residue into an output stream
              * @param stream Output stream
              */
            void WriteConnectSection(std::ofstream& stream);

            /*! \fn
              * A function in order to write the connectivity section of a specified residue into an output stream
              * @param stream Output stream
              */
            void WriteConnectivitySection(std::ofstream& stream, MolecularModeling::ResidueVector assembly_residues);

            /*! \fn
              * A function in order to write the hierarchy section of a specified residue into an output stream
              * @param stream Output stream
              */
            void WriteHierarchySection(std::ofstream& stream, OffFileResidueVector off_file_residues);

            /*! \fn
              * A function in order to write the name section of a specified residue into an output stream
              * @param stream Output stream
              */
            void WriteNameSection(std::ofstream& stream);

            /*! \fn
              * A function in order to write the position section of a specified residue into an output stream
              * @param stream Output stream
              */
            void WritePositionSection(std::ofstream& stream, MolecularModeling::ResidueVector assembly_residues, int CoordinateIndex);

            /*! \fn
              * A function in order to write the residue connect section of a specified residue into an output stream
              * @param stream Output stream
              */
            void WriteResidueConnectSection(std::ofstream& stream, MolecularModeling::ResidueVector assembly_residues);

            /*! \fn
              * A function in order to write the residues section of a specified residue into an output stream
              * @param stream Output stream
              */
            void WriteResiduesSection(std::ofstream& stream, MolecularModeling::ResidueVector assembly_residues);

            /*! \fn
              * A function in order to write the solvent cap section of a specified residue into an output stream
              * @param stream Output stream
              */
            void WriteSolventCapSection(std::ofstream& stream);

            /*! \fn
              * A function in order to write the velocities section of a specified residue into an output stream
              * @param stream Output stream
              */
            void WriteVelocitiesSection(std::ofstream& stream, OffFileResidueVector off_file_residues);


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

            std::string unit_name_;                  /*!< Name of OFF file unit >*/
            OffFileResidueVector off_file_residues_;   /*!< List of off file residues>*/
            AtomIndexMap  atom_index_map_;              /*!< A mapping of off file atom index to assembly atom index >*/
            AtomBondingMap atom_bonding_map_;            /*!< A mapping of off file atom index to bonded atom index >*/

    };
}

#endif // OFFFILE_HPP
