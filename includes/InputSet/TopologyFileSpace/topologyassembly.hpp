#ifndef TOPOLOGYASSEMBLY_HPP
#define TOPOLOGYASSEMBLY_HPP

#include <string>
#include <iostream>
#include <sstream>
#include <vector>
#include <map>

namespace TopologyFileSpace
{
    class TopologyResidue;

    class TopologyAssembly
    {
        public:
            //////////////////////////////////////////////////////////
            //                    TYPE DEFINITION                   //
            //////////////////////////////////////////////////////////
            /*! \typedef
              * A mapping between a residue name and the residue itself
              */
            typedef std::vector<TopologyResidue*> TopologyResidueVector;

            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            TopologyAssembly();

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
            /*! \fn
              * An accessor function in order to access to the assembly name
              * @return assembly_name_ Attribute of the current object of this class
              */
            std::string GetAssemblyName();
            /*! \fn
              * An accessor function in order to access to the residues
              * @return residues_ Attribute of the current object of this class
              */
            TopologyResidueVector GetResidues();
            /*! \fn
              * An accessor function in order to access to a residue of the current object using the index
              * @param index The index of the desired resiude
              * @return residue The residue with the given index
              */
            TopologyResidue* GetResidueByIndex(int index);
            /*! \fn
              * An accessor function in order to access to atom index of an atom of the current object using the name
              * @param atom_name The atom name of the desired atom index
              * @return index_ The index of the given atom name
              */
            int GetAtomIndexByName(std::string atom_name, int residue_index);

            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A mutator function in order to set the assembly name of the current object
              * Set the assembly_name_ attribute of the current topology assembly
              * @param assembly_name The assembly name attribute of the current object
              */
            void SetAssemblyName(std::string assembly_name);
            /*! \fn
              * A mutator function in order to set the residues of the current object
              * Set the residues_ attribute of the current topology assembly
              * @param residues The residues attribute of the current object
              */
            void SetResidues(TopologyResidueVector residues);
            void AddResidue(TopologyResidue* residue);

            //////////////////////////////////////////////////////////
            //                        FUNCTIONS                     //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the topology assembly contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cerr);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            std::string assembly_name_;             /*!< Name of the assembly made by the all residue involving in the assembly >*/
            TopologyResidueVector residues_;           /*!< Vector of residues that are involving in an assembly >*/

    };
}

#endif // TOPOLOGYASSEMBLY_HPP
