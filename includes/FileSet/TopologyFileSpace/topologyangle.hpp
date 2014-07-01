#ifndef TOPOLOGYANGLE_HPP
#define TOPOLOGYANGLE_HPP

#include <string>
#include <iostream>
#include <vector>

namespace TopologyFileSpace
{
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
              * An accessor function in order to access to boolean value that indicates if angle includes hydrogen
              * @return including_hydrogen_ attribute of the current object of this class
              */
            bool GetIncludingHydrogen();

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
              * A mutator function in order to set the boolean including hydrogen attribute of the current object
              * Set the including_hydrogen_ attribute of the current topology angle type
              * @param including_hydrogen The including hydrogen attribute of the current object
              */
            void SetIncludingHydrogen(bool including_hydrogen);

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
            bool including_hydrogen_;

    };
}

#endif // TOPOLOGYANGLE_HPP
