// Created by: Dave Montgomery

#ifndef PDBSOURCESECTION_HPP
#define PDBSOURCESECTION_HPP

#include <string>
#include <vector>
#include <iostream>

namespace PdbFileSpace
{
    class PdbSourceCard;

    class PdbSourceSection
    {
        public:
            //////////////////////////////////////////////////////////
            //                       TYPE DEFINITION                //
            //////////////////////////////////////////////////////////
            /*! \typedef
              * List of sources
              */
            typedef std::vector<PdbSourceCard*> SourceCardVector;


            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            PdbSourceSection();
            /*! \fn
              * Constructor with required parameters
              * @param stream_block
              */
            PdbSourceSection(std::stringstream& stream_block);

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
/** \addtogroup Molecular_Data_Structure
              * @{
              */
            /*! \fn
              * An accessor function in order to access to the source in a source card
              * @return source_ attribute of the current object of this class
              */
            SourceCardVector GetSourceCards();
/** @}*/
            //////////////////////////////////////////////////////////
            //                       MUTATOR                 PdbSourceSection       //
            //////////////////////////////////////////////////////////
/** \addtogroup Manipulators
              * @{
              */
            /*! \fn
              * A mutator function in order to set the source attribute of the current object
              * Set the source_ attribute of the current source card
              * @param source The source attribute of the current object
              */
            void SetSourceCards(SourceCardVector source);
            /*! \fn
              * A function in order to add the source attribute to the current object
              * Set the source_ attribute of the current source card
              * @param source The source attribute of the current object
              */
            void AddSourceCards(PdbSourceCard *source);
/** @}*/
            //////////////////////////////////////////////////////////
            //                        FUNCTIONS                     //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the source card contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cerr);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            SourceCardVector source_;      /*!< List of sources >*/

    };
}

#endif // PDBSOURCESECTION_HPP
