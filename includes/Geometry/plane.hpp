#ifndef PLANE_HPP
#define PLANE_HPP

#include <ostream>
#include <iomanip>
#include <iostream>

#include "coordinate.hpp"
#include "../common.hpp"

namespace Geometry
{
    class Plane
    {
        public:
            //////////////////////////////////////////////////////////
            //                       Constructor                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            Plane();

            //////////////////////////////////////////////////////////
            //                           ACCESSOR                   //
            //////////////////////////////////////////////////////////
            /*! \fn
              * An accessor function in order to access to the first vector of plane
              * @return v1_ first vector of plane
              */
            gmml::Vector GetV1();
            /*! \fn
              * An accessor function in order to access to the second vector of plane
              * @return v2_ second vector of plane
              */
            gmml::Vector GetV2();

            gmml::Vector GetUnitNormalVector();

            //////////////////////////////////////////////////////////
            //                           MUTATOR                    //
            //////////////////////////////////////////////////////////
            void SetV1(gmml::Vector v1);
            void SetV1(double x, double y, double z);
            void SetV2(gmml::Vector v2);
            void SetV2(double x, double y, double z);


            //////////////////////////////////////////////////////////
            //                         FUNCTIONS                    //
            //////////////////////////////////////////////////////////


            //////////////////////////////////////////////////////////
            //                     DISPLAY FUNCTIONS                //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out a plane
              * Print out the current plane in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cout);

    private:
            //////////////////////////////////////////////////////////
            //                         ATTRIBUTES                   //
            //////////////////////////////////////////////////////////
            gmml::Vector v1_;          /*!< first vector of plane*/
            gmml::Vector v2_;          /*!< second vector of plane */

    };
}

#endif // PLANE_HPP
