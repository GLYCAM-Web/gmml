#ifndef TOPOLOGYRESIDUE_HPP
#define TOPOLOGYRESIDUE_HPP

#include <string>
#include <iostream>
#include <vector>
#include <map>

namespace TopologyFileSpace
{
    class TopologyAtom;
    class TopologyBond;
    class TopologyAngle;
    class TopologyDihedral;

    class TopologyResidue
    {
        public:
            //////////////////////////////////////////////////////////
            //                    TYPE DEFINITION                   //
            //////////////////////////////////////////////////////////
            typedef std::map<std::string, TopologyAtom*> TopologyAtomMap;
            typedef std::map<std::vector<std::string>, TopologyBond*> TopologyBondMap;
            typedef std::map<std::vector<std::string>, TopologyAngle*> TopologyAngleMap;
            typedef std::map<std::vector<std::string>, TopologyDihedral*> TopologyDihedralMap;

            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            TopologyResidue();

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
            /*! \fn
              * An accessor function in order to access to the residue name
              * @return residue_name_ attribute of the current object of this class
              */
            std::string GetResidueName();
            /*! \fn
              * An accessor function in order to access to the atoms
              * @return atoms_ attribute of the current object of this class
              */
            TopologyAtomMap GetAtoms();
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
              * A mutator function in order to set the residue name of the current object
              * Set the residue_name_ attribute of the current topology residue
              * @param residue_name The residue name attribute of the current object
              */
            void SetResidueName(std::string residue_name);

            //////////////////////////////////////////////////////////
            //                        FUNCTIONS                     //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the topology residue contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cout);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            std::string residue_name_;
            TopologyAtomMap atoms_;
            TopologyBondMap bonds_;
            TopologyAngleMap angles_;
            TopologyDihedralMap dihedrals_;

    };
}

#endif // TOPOLOGYRESIDUE_HPP
