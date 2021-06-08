#ifndef GEOMETRY_TOPOLOGY_ROTATION_HPP
#define GEOMETRY_TOPOLOGY_ROTATION_HPP

#include <iostream>
#include "coordinate.hpp"

namespace GeometryTopology
{
    class Rotation
    {
        public:
            //////////////////////////////////////////////////////////
            //                       Constructor                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            Rotation();

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
              * A function to rotate list of coordinate vector
              * @param pivot_point A coordinate point in space around which the coordonates will rotate
              * @param direction_vector A coordinate poimnt in space that will give direction coordinate rotation
              * @param rotation_angle An angle of rotation
              * @param coordinate_set List of coordinates that need to be rotated by an angle
              * @return A list of coordinates rotated by an angle
              */
            GeometryTopology::CoordinateVector RotateCoordinates(GeometryTopology::Coordinate* pivot_point,GeometryTopology::Coordinate* direction_point,double rotation_angle,GeometryTopology::CoordinateVector coordinate_set);

            //////////////////////////////////////////////////////////
            //                     DISPLAY FUNCTIONS                //
            //////////////////////////////////////////////////////////

       private:
            //////////////////////////////////////////////////////////
            //                         ATTRIBUTES                   //
            //////////////////////////////////////////////////////////
        };
}
#endif // GEOMETRY_TOPOLOGY_ROTATION_HPP
