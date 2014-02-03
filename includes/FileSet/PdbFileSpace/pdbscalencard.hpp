// Created by: Delaram Rahbarinia
// Modified by: Alireza Khatamian, Delaram Rahbarinia

#ifndef PDBSCALENCARD_HPP
#define PDBSCALENCARD_HPP

#include <string>
#include <vector>
#include <sstream>

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
            void SetScaleN(const ScaleNVector scale_n);
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

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            ScaleNVector scale_n_;

    };
}

#endif // PDBSCALENCARD_HPP
