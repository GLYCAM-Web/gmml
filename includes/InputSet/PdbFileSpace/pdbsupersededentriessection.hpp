// Created by: Dave Montgomery

#ifndef PDBSUPERSEDEDENTRIESSECTION_HPP
#define PDBSUPERSEDEDENTRIESSECTION_HPP

#include <string>
#include <vector>
#include <iostream>

namespace PdbFileSpace
{
    class PdbSupersededEntriesCard;

    class PdbSupersededEntriesSection
    {
        public:
            //////////////////////////////////////////////////////////
            //                       TYPE DEFINITION                //
            //////////////////////////////////////////////////////////
            /*! \typedef
              * List of superseded_entriess
              */
            typedef std::vector<PdbSupersededEntriesCard*> SupersededEntriesCardVector;


            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            PdbSupersededEntriesSection();
            /*! \fn
              * Constructor with required parameters
              * @param stream_block
              */
            PdbSupersededEntriesSection(std::stringstream& stream_block);

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
/** \addtogroup Molecular_Data_Structure
              * @{
              */
            /*! \fn
              * An accessor function in order to access to the superseded_entries in a superseded_entries card
              * @return superseded_entries_ attribute of the current object of this class
              */
            SupersededEntriesCardVector GetSupersededEntriesCards();
/** @}*/
            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
/** \addtogroup Manipulators
              * @{
              */
            /*! \fn
              * A mutator function in order to set the superseded_entries attribute of the current object
              * Set the superseded_entries_ attribute of the current superseded_entries card
              * @param superseded_entries The superseded_entries attribute of the current object
              */
            void SetSupersededEntriesCards(SupersededEntriesCardVector superseded_entries);
            /*! \fn
              * A function in order to add the superseded_entries attribute to the current object
              * Set the superseded_entries_ attribute of the current superseded_entries card
              * @param superseded_entries The superseded_entries attribute of the current object
              */
            void AddSupersededEntriesCards(PdbSupersededEntriesCard *superseded_entries);
/** @}*/
            //////////////////////////////////////////////////////////
            //                        FUNCTIONS                     //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the superseded_entries card contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cerr);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            SupersededEntriesCardVector superseded_entries_;      /*!< List of superseded_entriess >*/

    };
}

#endif // PDBSUPERSEDEDENTRIESSECTION_HPP
