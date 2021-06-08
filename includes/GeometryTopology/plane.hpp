#ifndef PLANE_HPP
#define PLANE_HPP

#include <ostream>
#include <iomanip>
#include <iostream>

#include "coordinate.hpp"
#include "../common.hpp"

namespace GeometryTopology
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
/** \addtogroup Molecular_Data_Structure
               * @{
               */
            /*! \fn
              * An accessor function in order to access to the first vector of plane
              * @return v1_ first vector of plane
              */
            Coordinate GetV1();
            /*! \fn
              * An accessor function in order to access to the second vector of plane
              * @return v2_ second vector of plane
              */
            Coordinate GetV2();

            Coordinate GetUnitNormalVector();
/** @}*/
            //////////////////////////////////////////////////////////
            //                           MUTATOR                    //
            //////////////////////////////////////////////////////////
/** \addtogroup Manipulators
               * @{
               */
            void SetV1(Coordinate v1);
            void SetV1(double x, double y, double z);
            void SetV2(Coordinate v2);
            void SetV2(double x, double y, double z);
/** @}*/

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
            void Print(std::ostream& out = std::cerr);

    private:
            //////////////////////////////////////////////////////////
            //                         ATTRIBUTES                   //
            //////////////////////////////////////////////////////////
            Coordinate v1_;          /*!< first vector of plane*/
            Coordinate v2_;          /*!< second vector of plane */

    };
}

#endif // PLANE_HPP
