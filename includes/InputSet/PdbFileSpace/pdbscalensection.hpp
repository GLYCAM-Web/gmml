// Created by: Delaram Rahbarinia
// Modified by: Alireza Khatamian, Delaram Rahbarinia, Dave Montgomery

#ifndef PDBSCALENSECTION_HPP
#define PDBSCALENSECTION_HPP

#include <string>
#include <vector>
#include <sstream>
#include <iostream>

namespace PdbFileSpace
{
    class PdbScaleNCard;

    class PdbScaleNSection
    {
        public:
            //////////////////////////////////////////////////////////
            //                       TYPE DEFINITION                //
            //////////////////////////////////////////////////////////
            /*! \typedef
              * Vector of scales
              */
            typedef std::vector< PdbScaleNCard* > ScaleNCardVector;

            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            PdbScaleNSection();
            /*! \fn
              * Constructor with required parameters
              * @param stream_block
              */
            PdbScaleNSection(std::stringstream& stream_block);

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
/** \addtogroup Molecular_Data_Structure
               * @{
               */
            /*! \fn
              * An accessor function in order to access to the scale n in a scale n card
              * @return scale_n_ attribute of the current object of this class
              */
            ScaleNCardVector GetScaleNCard();
/** @}*/
            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
/** \addtogroup Manipulators
               * @{
               */
            /*! \fn
              * A mutator function in order to set the scale n attribute of the current object
              * Set the scale_n_ attribute of the current scale n card
              * @param scale_n The scale n attribute of the current object
              */
            void SetScaleNCard(ScaleNCardVector scale_n);
            /*! \fn
              * A function in order to add the record name to the current object
              * Set the scale_ attribute of the current scale n card
              * @param scale The scale attribute of the current object
              */
            void AddScaleNCard(PdbScaleNCard* scale);
/** @}*/
            //////////////////////////////////////////////////////////
            //                        FUNCTIONS                     //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the pdb scale n card contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cerr);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            ScaleNCardVector scale_n_;      /*!< Vector of scales are in scalen card of a pdb file >*/

    };
}

#endif // PDBSCALENSECTION_HPP
