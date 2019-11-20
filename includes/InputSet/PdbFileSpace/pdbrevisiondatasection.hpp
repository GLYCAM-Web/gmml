// Created by: Dave Montgomery

#ifndef PDBREVISIONDATASECTION_HPP
#define PDBREVISIONDATASECTION_HPP

#include <string>
#include <vector>
#include <iostream>

namespace PdbFileSpace
{
    class PdbRevisionDataCard;

    class PdbRevisionDataSection
    {
        public:
            //////////////////////////////////////////////////////////
            //                       TYPE DEFINITION                //
            //////////////////////////////////////////////////////////
            /*! \typedef
              * List of revision_datas
              */
            typedef std::vector<PdbRevisionDataCard*> RevisionDataCardVector;


            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            PdbRevisionDataSection();
            /*! \fn
              * Constructor with required parameters
              * @param stream_block
              */
            PdbRevisionDataSection(std::stringstream& stream_block);

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
/** \addtogroup Molecular_Data_Structure
              * @{
              */
            /*! \fn
              * An accessor function in order to access to the revision_data in a revision_data card
              * @return revision_data_ attribute of the current object of this class
              */
            RevisionDataCardVector GetRevisionDataCards();
/** @}*/
            //////////////////////////////////////////////////////////
            //                       MUTATOR                 PdbRevisionDataSection       //
            //////////////////////////////////////////////////////////
/** \addtogroup Manipulators
              * @{
              */
            /*! \fn
              * A mutator function in order to set the revision_data attribute of the current object
              * Set the revision_data_ attribute of the current revision_data card
              * @param revision_data The revision_data attribute of the current object
              */
            void SetRevisionDataCards(RevisionDataCardVector revision_data);
            /*! \fn
              * A function in order to add the revision_data attribute to the current object
              * Set the revision_data_ attribute of the current revision_data card
              * @param revision_data The revision_data attribute of the current object
              */
            void AddRevisionDataCards(PdbRevisionDataCard *revision_data);
/** @}*/
            //////////////////////////////////////////////////////////
            //                        FUNCTIONS                     //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the revision_data card contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cerr);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            RevisionDataCardVector revision_data_;      /*!< List of revision_datas >*/

    };
}

#endif // PDBREVISIONDATASECTION_HPP
