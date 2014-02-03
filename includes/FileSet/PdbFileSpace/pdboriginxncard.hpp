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
            typedef std::vector< PdbOriginXn* > OriginXnVector;

            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            PdbOriginXnCard();
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
            void SetOriginXN(const OriginXnVector origin_x_n);
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
            void Print(std::ostream& out = std::cout);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            OriginXnVector origin_x_n_;

    };
}

#endif // PDBORIGINXNCARD_HPP
