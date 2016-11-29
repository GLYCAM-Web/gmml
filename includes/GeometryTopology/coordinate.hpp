#ifndef COORDINATE_HPP
#define COORDINATE_HPP

#include <ostream>
#include <iomanip>
#include <iostream>
#include <vector>

namespace GeometryTopology
{
    class Coordinate
    {
        public:
            //////////////////////////////////////////////////////////
            //                    TYPE DEFINITION                   //
            //////////////////////////////////////////////////////////
            typedef std::vector<Coordinate*> CoordinateVector;
            //////////////////////////////////////////////////////////
            //                       Constructor                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            Coordinate();
            /*! \fn
              * Constructor with given three double values
              * @param x A double value on X axis in cartesian coordinate
              * @param y A double value on Y axis in cartesian coordinate
              * @param z A double value on Z axis in cartesian coordinate
              */
            Coordinate(double x, double y, double z);
            /*! \fn
              * Copy constructor
              * @param coordinate A coordinate to be copied to another instance
              */
            Coordinate(const Coordinate& coordinate);            
            Coordinate(Coordinate* coordinate);


            //////////////////////////////////////////////////////////
            //                           ACCESSOR                   //
            //////////////////////////////////////////////////////////
            /*! \fn
              * An accessor function in order to access to X coordinate attribute of the current object
              * The attribute is set by the contents of the given file
              * @return x_ of the current object of this class
              */
            double GetX();
            /*! \fn
              * An accessor function in order to access to Y coordinate attribute of the current object
              * The attribute is set by the contents of the given file
              * @return y_ of the current object of this class
              */
            double GetY();
            /*! \fn
              * An accessor function in order to access to Z coordinate attribute of the current object
              * The attribute is set by the contents of the given file
              * @return z_ of the current object of this class
              */
            double GetZ();

            //////////////////////////////////////////////////////////
            //                           MUTATOR                    //
            //////////////////////////////////////////////////////////
            void SetX(double x);
            void SetY(double y);
            void SetZ(double z);

            //////////////////////////////////////////////////////////
            //                         FUNCTIONS                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Translate a coordinate with respect to given coordinate by x, y and z
              * @param x A double value on X axis for the origin of translation
              * @param y A double value on Y axis for the origin of translation
              * @param z A double value on Z axis for the origin of translation
              */
            void Translate(double x, double y, double z);

            /*! \fn
              * Compare current object with the given coordinate and return true if they are the same
              * @param coordinate A coordinate to be compared with the current object
              * @return A boolean value as a result of comparison
              */
            bool CompareTo(Coordinate coordinate);

            /*! \fn
              * Calculate the distance between current object coordinate and the given coordinate
              * @param coordinate A coordinate to calculate the distance between that and the current object
              * @return A double value as a distance between two coordinates
              */
            double Distance(Coordinate coordinate);

            /*! \fn
              * Calculate the length of the vector
              * @return Length of the vector
              */
            double length();
            /*! \fn
              * Normalize the vector and update it on the current object
              */
            void Normalize();
            /*! \fn
              * Calculate the dot product of current object and the given one
              * @param coordinate Second vector in the dot product operator
              * @return Dot product value of two vectors
              */
            double DotProduct(Coordinate coordinate);
            /*! \fn
              * Calculate the vector from cross producting the current object and the given one
              * @param coordinate Second vector in the cross product operator
              * @return Cross product vector of two vectors
              */
            void CrossProduct(Coordinate coordinate);
            /*! \fn
              * Add current object to the given one and update the current one
              * @param coordinate Second vector in the add operator
              */
            void operator+(Coordinate coordinate);
            void operator +(double addition);
            /*! \fn
              * Subtract current object to the given one and update the current one
              * @param coordinate Second vector in the subtract operator
              */
            void operator-(Coordinate coordinate);
            void operator /(Coordinate coordinate);
            void operator/(double divisor);
            void operator*(double multiplier);
            void TranslateAll(CoordinateVector coordinate_set, double margin, int pos);
            void RotateAngularAll(CoordinateVector coordinate_set, double angle, int pos);
            void RotateTorsionalAll(CoordinateVector coordinate_set, double torsion, int pos);
            //////////////////////////////////////////////////////////
            //                     DISPLAY FUNCTIONS                //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out a coordinate
              * Print out the current coordinate in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cout);

	private:
            //////////////////////////////////////////////////////////
            //                         ATTRIBUTES                   //
            //////////////////////////////////////////////////////////
            double x_;          /*!< x */
            double y_;          /*!< y */
            double z_;          /*!< z */
    };
}

#endif // COORDINATE_HPP
