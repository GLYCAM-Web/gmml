// Created by: Delaram Rahbarinia
// Modified by: Alireza Khatamian, Delaram Rahbarinia

#ifndef PDBORIGINXNCARD_HPP
#define PDBORIGINXNCARD_HPP

#include <string>
#include <vector>
#include <sstream>
#include <iostream>

namespace PdbFileSpace
{
    class PdbOriginXn;

    class PdbOriginXnCard
    {
        public:
            //////////////////////////////////////////////////////////
            //                       TYPE DEFINITION                //
            //////////////////////////////////////////////////////////
            /*! \typedef
              * Vector of origins
              */
            typedef std::vector< PdbOriginXn* > OriginXnVector;

            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            PdbOriginXnCard();
            /*! \fn
              * Constructor with required parameters
              * @param stream_block
              */
            PdbOriginXnCard(std::stringstream& stream_block);

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
            /*! \fn
              * An accessor function in order to access to the origin X N attribute in a origin X N card
              * @return origin_x_n_ attribute of the current object of this class
              */
            OriginXnVector GetOriginXN();

            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A mutator function in order to set the origin X N of the current object
              * Set the origin_x_n_ attribute of the current origin X N
              * @param origin_x_n The origin X N of the current object
              */
            void SetOriginXN(OriginXnVector origin_x_n);
            /*! \fn
              * A function in order to add the origin to the current object
              * Set the origin_ attribute of the current origin X N
              * @param origin The origin of the current object
              */
            void AddOriginXN(PdbOriginXn* origin);

            //////////////////////////////////////////////////////////
            //                        FUNCTIONS                     //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the origin xn card contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cout);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            OriginXnVector origin_x_n_;     /*!< Vector of origins >*/

    };
}

#endif // PDBORIGINXNCARD_HPP
