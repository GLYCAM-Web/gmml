#ifndef TOPOLOGYANGLE_HPP
#define TOPOLOGYANGLE_HPP

#include <string>
#include <iostream>
#include <vector>

namespace TopologyFileSpace
{
    class TopologyAngleType;

    class TopologyAngle
    {
        public:
            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            TopologyAngle();
            /*! \fn
              * Constructor with required parameters
              * @param angle_atoms Names of atoms involving in an angle (three bonded atom names)
              * @param residue_names Names of residues of the atoms that are involving in an angle
              */
            TopologyAngle(std::vector<std::string> angle_atoms, std::vector<std::string> residue_names);

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
            /*! \fn
              * An accessor function in order to access to the angles
              * @return angles_ attribute of the current object of this class
              */
            std::vector<std::string> GetAngles();
            /*! \fn
              * An accessor function in order to access to the angle type
              * @return angle_type_ attribute of the current object of this class
              */
            TopologyAngleType* GetAngleType();
            /*! \fn
              * An accessor function in order to access to boolean value that indicates if angle includes hydrogen
              * @return including_hydrogen_ attribute of the current object of this class
              */
            bool GetIncludingHydrogen();
            /*! \fn
              * An accessor function in order to access to the residue names
              * @return residue_names_ attribute of the current object of this class
              */
            std::vector<std::string> GetResidueNames();

            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A mutator function in order to set the angles of the current object
              * Set the angles_ attribute of the current topology angle
              * @param angles The angles attribute of the current object
              */
            void SetAngles(std::vector<std::string> angles);
            /*! \fn
              * A mutator function in order to set the angle type of the current object
              * Set the angle_type_ attribute of the current topology angle
              * @param angle_type The angle type attribute of the current object
              */
            void SetAnlgeType(TopologyAngleType* angle_type);
            /*! \fn
              * A mutator function in order to set the boolean including hydrogen attribute of the current object
              * Set the including_hydrogen_ attribute of the current topology angle
              * @param including_hydrogen The including hydrogen attribute of the current object
              */
            void SetIncludingHydrogen(bool including_hydrogen);
            /*! \fn
              * A mutator function in order to set residue names of the current object
              * Set the residue_names_ attribute of the current topology angle
              * @param residue_names The residue names of the current object
              */
            void SetResidueNames(std::vector<std::string> residue_names);

            //////////////////////////////////////////////////////////
            //                        FUNCTIONS                     //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the topology angle contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cerr);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            std::vector<std::string> angles_;               /*!< Name of atoms involving in an angle in a topology file >*/
            TopologyAngleType* angle_type_;                 /*!< Type of the angle corresponding to the angle type definition >*/
            bool including_hydrogen_;                       /*!< Indicates that the angle has a hydrogen atom or not >*/
            std::vector<std::string> residue_names_;        /*!< Names of residues of the atoms that are invovling in an angle >*/

    };
}

#endif // TOPOLOGYANGLE_HPP
