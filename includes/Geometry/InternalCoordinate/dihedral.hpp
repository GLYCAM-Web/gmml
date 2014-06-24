#ifndef DIHEDRAL_HPP
#define DIHEDRAL_HPP

#include <string>
#include <iostream>
#include <vector>

namespace Geometry
{
    class Coordinate;
    class Dihedral
    {
        public:
            //////////////////////////////////////////////////////////
            //                    TYPE DEFINITION                   //
            //////////////////////////////////////////////////////////
            typedef std::vector<Coordinate*> CoordinateVector;

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
            /*! \fn
              * An accessor function in order to access to the coordinates
              * @return coordinates_ attribute of the current object of this class
              */
            CoordinateVector GetCoordinates();
            /*! \fn
              * An accessor function in order to access to the torsion
              * @return torsion_ attribute of the current object of this class
              */
            double GetTorsion();

            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A mutator function in order to set the coordinates of the current object
              * Set the coordinates_ attribute of the current dihedral
              * @param coordinates The coordinates attribute of the current object
              */
            void SetCoordinates(CoordinateVector coordinates);
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
            CoordinateVector coordinates_;
            double torsion_;

    };
}

#endif // DIHEDRAL_HPP
