// Created by: Delaram Rahbarinia
// Modified by: Alireza Khatamian, Delaram Rahbarinia, Dave Montgomery

#ifndef PDBORIGINXNSECTION_HPP
#define PDBORIGINXNSECTION_HPP

#include <string>
#include <vector>
#include <sstream>
#include <iostream>

namespace PdbFileSpace
{
    class PdbOriginXnCard;

    class PdbOriginXnSection
    {
        public:
            //////////////////////////////////////////////////////////
            //                       TYPE DEFINITION                //
            //////////////////////////////////////////////////////////
            /*! \typedef
              * Vector of origins
              */
            typedef std::vector< PdbOriginXnCard* > OriginXnCardVector;

            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            PdbOriginXnSection();
            /*! \fn
              * Constructor with required parameters
              * @param stream_block
              */
            PdbOriginXnSection(std::stringstream& stream_block);

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
/** \addtogroup Molecular_Data_Structure
               * @{
               */
            /*! \fn
              * An accessor function in order to access to the origin X N attribute in a origin X N card
              * @return origin_x_n_ attribute of the current object of this class
              */
            OriginXnCardVector GetOriginXN();
/** @}*/
            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
/** \addtogroup Manipulators
               * @{
               */
            /*! \fn
              * A mutator function in order to set the origin X N of the current object
              * Set the origin_x_n_ attribute of the current origin X N
              * @param origin_x_n The origin X N of the current object
              */
            void SetOriginXN(OriginXnCardVector origin_x_n);
            /*! \fn
              * A function in order to add the origin to the current object
              * Set the origin_ attribute of the current origin X N
              * @param origin The origin of the current object
              */
            void AddOriginXN(PdbOriginXnCard* origin);
/** @}*/
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
            void Print(std::ostream& out = std::cerr);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            OriginXnCardVector origin_x_n_;     /*!< Vector of origins >*/

    };
}

#endif // PDBORIGINXNSECTION_HPP
