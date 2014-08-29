// Created by: Delaram Rahbarinia
// Modified by: Alireza Khatamian, Delaram Rahbarinia

#ifndef PDBSCALENCARD_HPP
#define PDBSCALENCARD_HPP

#include <string>
#include <vector>
#include <sstream>
#include <iostream>

namespace PdbFileSpace
{
    class PdbScaleN;

    class PdbScaleNCard
    {
        public:
            //////////////////////////////////////////////////////////
            //                       TYPE DEFINITION                //
            //////////////////////////////////////////////////////////
            typedef std::vector< PdbScaleN* > ScaleNVector;

            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            PdbScaleNCard();
            /*! \fn
              * Constructor with required parameters
              * @param stream_block
              */
            PdbScaleNCard(std::stringstream& stream_block);

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
            /*! \fn
              * An accessor function in order to access to the scale n in a scale n card
              * @return scale_n_ attribute of the current object of this class
              */
            ScaleNVector GetScaleN();

            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A mutator function in order to set the scale n attribute of the current object
              * Set the scale_n_ attribute of the current scale n card
              * @param scale_n The scale n attribute of the current object
              */
            void SetScaleN(ScaleNVector scale_n);
            /*! \fn
              * A function in order to add the record name to the current object
              * Set the scale_ attribute of the current scale n card
              * @param scale The scale attribute of the current object
              */
            void AddScaleN(PdbScaleN* scale);

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
            void Print(std::ostream& out = std::cout);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            ScaleNVector scale_n_;

    };
}

#endif // PDBSCALENCARD_HPP
