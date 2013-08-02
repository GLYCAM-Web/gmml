#ifndef COORDINATE_HPP
#define COORDINATE_HPP

#include <ostream>
#include <iomanip>

namespace Geometry
{
    class Coordinate
    {
        public:
            /////////////////////////////// CONSTRUCTOR ////////////////////////////////
            Coordinate();
            Coordinate(double x, double y, double z);
            Coordinate(const Coordinate& coordinate);

            /////////////////////////////// FUNCTIONS //////////////////////////////////
            void Translate(double x, double y, double z);

            ///////////////////////////// DISPLAY FUNCTION /////////////////////////////
            void Print(std::ostream& out);

            /////////////////////////////// ATTRIBUTES /////////////////////////////////
            double x_;
            double y_;
            double z_;
    };
}

#endif // COORDINATE_HPP
