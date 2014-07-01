#ifndef TOPOLOGYANGLE_HPP
#define TOPOLOGYANGLE_HPP

#include <string>
#include <iostream>
#include <vector>

namespace TopologyFileSpace
{
    /*! \enum
      * Topology positional angle type enumerator
      */
        enum TopologyPositionalAngleType
        {
            INTER_RESIDUE_ANGLE_ = 1,
            INTRA_RESIDUE_ANGLE_ = 2
        };

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

            TopologyAngle(std::vector<std::string> angle_atoms);

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
              * An accessor function in order to access to positional angle type
              * @return posotional_type_ attribute of the current object of this class
              */
            TopologyPositionalAngleType GetPositionalType();

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
              * Set the angle_type_ attribute of the current topology angle type
              * @param angle_type The angle type attribute of the current object
              */
            void SetAnlgeType(TopologyAngleType* angle_type);
            /*! \fn
              * A mutator function in order to set the boolean including hydrogen attribute of the current object
              * Set the including_hydrogen_ attribute of the current topology angle type
              * @param including_hydrogen The including hydrogen attribute of the current object
              */
            void SetIncludingHydrogen(bool including_hydrogen);
            /*! \fn
              * A mutator function in order to set the positional type of the current object
              * Set the positional_type_ attribute of the current topology angle type
              * @param positional_type The positional_type attribute of the current object
              */
            void SetPositionalType(TopologyPositionalAngleType positional_type);

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
            void Print(std::ostream& out = std::cout);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            std::vector<std::string> angles_;
            TopologyAngleType* angle_type_;
            bool including_hydrogen_;
            TopologyPositionalAngleType positional_type_;

    };
}

#endif // TOPOLOGYANGLE_HPP
