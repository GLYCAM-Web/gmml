#ifndef GEOMETRY_HPP
#define GEOMETRY_HPP

#include <iostream>
#include "../../../includes/GeometryTopology/coordinate.hpp"

namespace GeometryOperation
{
    class Geometry
    {
        public:
            //////////////////////////////////////////////////////////
            //                    TYPE DEFINITION                   //
            //////////////////////////////////////////////////////////

            /*! \typedef
            * List of GeometryTopology::Coordinates
            */
            typedef std::vector<GeometryTopology::Coordinate*> CoordinateVector;

            //////////////////////////////////////////////////////////
            //                       Constructor                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            Geometry();

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
            CoordinateVector RotateCoordinates(GeometryTopology::Coordinate* pivot_point,GeometryTopology::Coordinate* direction_point,double rotation_angle,CoordinateVector coordinate_set);

            //////////////////////////////////////////////////////////
            //                     DISPLAY FUNCTIONS                //
            //////////////////////////////////////////////////////////

       private:
            //////////////////////////////////////////////////////////
            //                         ATTRIBUTES                   //
            //////////////////////////////////////////////////////////
        };
}
#endif // GEOMETRY_HPP
