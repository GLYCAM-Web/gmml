#ifndef DIHEDRAL_HPP
#define DIHEDRAL_HPP

#include <string>
#include <iostream>
#include <vector>
#include "coordinate.hpp"

namespace GeometryTopology
{
    // class Coordinate;
    class Dihedral
    {
        public:
            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            Dihedral();

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
            GeometryTopology::Coordinate::CoordinateVector GetCoordinates();
            /*! \fn
              * An accessor function in order to access to the torsion
              * @return torsion_ attribute of the current object of this class
              */
            double GetTorsion();
/** @}*/
            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
/** \addtogroup Manipulators
               * @{
               */
            /*! \fn
              * A mutator function in order to set the coordinates of the current object
              * Set the coordinates_ attribute of the current dihedral
              * @param coordinates The coordinates attribute of the current object
              */
            void SetCoordinates(GeometryTopology::Coordinate::CoordinateVector coordinates);
            /*! \fn
              * A function in order to add the coordinate to the current object
              * Set the coordinates_ attribute of the current dihedral
              * @param coordinate The atom of the current object
              */
            void AddCoordinate(Coordinate* coordinate);
            /*! \fn
              * A mutator function in order to set the torsion of the current object
              * Set the torsion_ attribute of the current dihedral
              * @param torsion The torsion attribute of the current object
              */
            void SetTorsion(double torsion);
/** @}*/
            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the dihedral contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cout);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            GeometryTopology::Coordinate::CoordinateVector coordinates_;          /*!< Vector of four coordinates representing a dihedral between four points >*/
            double torsion_;                        /*!< Torsion between the four points >*/

    };
}

#endif // DIHEDRAL_HPP
