// Created by: Dave Montgomery

#ifndef PDBSEQUENCEADVANCEDSECTION_HPP
#define PDBSEQUENCEADVANCEDSECTION_HPP

#include <string>
#include <vector>
#include <iostream>

namespace PdbFileSpace
{
    class PdbSequenceAdvancedCard;

    class PdbSequenceAdvancedSection
    {
        public:
            //////////////////////////////////////////////////////////
            //                       TYPE DEFINITION                //
            //////////////////////////////////////////////////////////
            /*! \typedef
              * List of sequence advanced cards
              */
            typedef std::vector<PdbSequenceAdvancedCard*> SequenceAdvancedCardVector;


            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            PdbSequenceAdvancedSection();
            /*! \fn
              * Constructor with required parameters
              * @param stream_block
              */
            PdbSequenceAdvancedSection(std::stringstream& stream_block);

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
/** \addtogroup Molecular_Data_Structure
              * @{
              */
            /*! \fn
              * An accessor function in order to access to the sequence_advanced in a sequence_advanced card
              * @return sequence_advanced_ attribute of the current object of this class
              */
            SequenceAdvancedCardVector GetSequenceAdvancedCards();
/** @}*/
            //////////////////////////////////////////////////////////
            //                       MUTATOR                 PdbSequenceAdvancedSection       //
            //////////////////////////////////////////////////////////
/** \addtogroup Manipulators
              * @{
              */
            /*! \fn
              * A mutator function in order to set the sequence_advanced attribute of the current object
              * Set the sequence_advanced_ attribute of the current sequence_advanced card
              * @param sequence_advanced The sequence_advanced attribute of the current object
              */
            void SetSequenceAdvancedCards(SequenceAdvancedCardVector sequence_advanced);
            /*! \fn
              * A function in order to add the sequence_advanced attribute to the current object
              * Set the sequence_advanced_ attribute of the current sequence_advanced card
              * @param sequence_advanced The sequence_advanced attribute of the current object
              */
            void AddSequenceAdvancedCards(PdbSequenceAdvancedCard *sequence_advanced);
/** @}*/
            //////////////////////////////////////////////////////////
            //                        FUNCTIONS                     //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the sequence_advanced section contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cerr);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            SequenceAdvancedCardVector sequence_advanced_;      /*!< List of sequence advanced cards >*/

    };
}

#endif // PDBSEQUENCEADVANCEDSECTION_HPP
