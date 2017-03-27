#ifndef DISTANCE_HPP
#define DISTANCE_HPP

#include <string>
#include <iostream>
#include <vector>

namespace GeometryTopology
{
    class Coordinate;
    class Distance
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
            Distance();

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
/** \addtogroup Molecular_Data_Structure
               * @{
               */
            /*! \fn
              * An accessor function in order to access to the coordinates
              * @return coordinates_ attribute of the current object of this class
              */
            CoordinateVector GetCoordinates();
            /*! \fn
              * An accessor function in order to access to the length
              * @return length_ attribute of the current object of this class
              */
            double GetLength();
/** @}*/
            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
/** \addtogroup Manipulators
               * @{
               */
            /*! \fn
              * A mutator function in order to set the coordinates of the current object
              * Set the coordinates_ attribute of the current distance
              * @param coordinates The coordinates attribute of the current object
              */
            void SetCoordinates(CoordinateVector coordinates);
            /*! \fn
              * A function in order to add the coordinate to the current object
              * Set the coordinates_ attribute of the current distance
              * @param coordinate The coordinate of the current object
              */
            void AddCoordinate(Coordinate* coordinate);
            /*! \fn
              * A mutator function in order to set the length of the current object
              * Set the length_ attribute of the current distance
              * @param length The length attribute of the current object
              */
            void SetLenght(double length);
/** @}*/
            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the distance contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cout);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            CoordinateVector coordinates_;              /*!< Vector of two coordinates representing a bond between two points >*/
            double length_;                             /*!< Length of the bond between two coordinates >*/

    };
}

#endif // DISTANCE_HPP
