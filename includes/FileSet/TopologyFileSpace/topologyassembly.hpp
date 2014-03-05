#ifndef TOPOLOGYASSEMBLY_HPP
#define TOPOLOGYASSEMBLY_HPP

#include <string>
#include <iostream>
#include <vector>
#include <map>

namespace TopologyFileSpace
{
    class TopologyBond;
    class TopologyAngle;
    class TopologyDihedral;
    class TopologyResidue;

    class TopologyAssembly
    {
        public:
            //////////////////////////////////////////////////////////
            //                    TYPE DEFINITION                   //
            //////////////////////////////////////////////////////////
            typedef std::map<std::string, TopologyResidue*> TopologyResidueMap;
            typedef std::map<std::vector<std::string>, TopologyBond*> TopologyBondMap;
            typedef std::map<std::vector<std::string>, TopologyAngle*> TopologyAngleMap;
            typedef std::map<std::vector<std::string>, TopologyDihedral*> TopologyDihedralMap;

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
              * @return assembly_name_ attribute of the current object of this class
              */
            std::string GetAssemblyName();
            /*! \fn
              * An accessor function in order to access to the residues
              * @return residues_ attribute of the current object of this class
              */
            TopologyResidueMap GetResidues();
            /*! \fn
              * An accessor function in order to access to the bonds
              * @return bonds_ attribute of the current object of this class
              */
            TopologyBondMap GetBonds();
            /*! \fn
              * An accessor function in order to access to the angles
              * @return angles_ attribute of the current object of this class
              */
            TopologyAngleMap GetAngles();
            /*! \fn
              * An accessor function in order to access to the dihedrals
              * @return dihedrals_ attribute of the current object of this class
              */
            TopologyDihedralMap GetDihedrals();

            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A mutator function in order to set the assembly name of the current object
              * Set the assembly_name_ attribute of the current topology assembly
              * @param assembly_name The assembly name attribute of the current object
              */
            void SetAssemblyName(std::string assembly_name);

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
            void Print(std::ostream& out = std::cout);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            std::string assembly_name_;
            TopologyResidueMap residues_;
            TopologyBondMap bonds_;
            TopologyAngleMap angles_;
            TopologyDihedralMap dihedrals_;

    };
}

#endif // TOPOLOGYASSEMBLY_HPP
