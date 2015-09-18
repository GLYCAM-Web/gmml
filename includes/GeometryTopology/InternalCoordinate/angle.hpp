#ifndef ANGLE_HPP
#define ANGLE_HPP

#include <string>
#include <iostream>
#include <vector>

namespace GeometryTopology
{
    class Coordinate;
    class Angle
    {
        public:
            //////////////////////////////////////////////////////////
            //                    TYPE DEFINITION                   //
            //////////////////////////////////////////////////////////
            /*! \typedef
              * Vector of coordinates
              */
            typedef std::vector<Coordinate*> CoordinateVector;

            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            Angle();

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
            /*! \fn
              * An accessor function in order to access to the coordinates
              * @return coordinates_ attribute of the current object of this class
              */
            CoordinateVector GetCoordinates();
            /*! \fn
              * An accessor function in order to access to the angle
              * @return angle_ attribute of the current object of this class
              */
            double GetAngle();

            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A mutator function in order to set the coordinates of the current object
              * Set the coordinates_ attribute of the current angle
              * @param cordinates The coordinates attribute of the current object
              */
            void SetCoordinates(CoordinateVector coordinates);
            /*! \fn
              * A function in order to add the coordinate to the current object
              * Set the coordinates_ attribute of the current angle
              * @param coordinate The coordinate of the current object
              */
            void AddCoordinate(Coordinate* coordinate);
            /*! \fn
              * A mutator function in order to set the angle of the current object
              * Set the angle_ attribute of the current angle
              * @param angle The angle attribute of the current object
              */
            void SetAngle(double angle);

            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the angle contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cout);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            CoordinateVector coordinates_;              /*!< Vector of three coordinates representing an angle between three points >*/
            double angle_;                              /*!< Value of angle between three points >*/

    };
}

#endif // ANGLE_HPP
