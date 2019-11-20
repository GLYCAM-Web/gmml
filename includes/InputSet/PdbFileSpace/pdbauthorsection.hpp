// Created by: Dave Montgomery

#ifndef PDBAUTHORSECTION_HPP
#define PDBAUTHORSECTION_HPP

#include <string>
#include <sstream>
#include <iostream>

namespace PdbFileSpace
{
    class PdbAuthorSection
    {
        public:
            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            PdbAuthorSection();
            /*! \fn
              * Constructor with required parameters
              * @param record_name Name for a author card record which appears in the first column of each line in a pdb file
              * @param author Author of a pdb file
              */
            PdbAuthorSection(const std::string& record_name, const std::string& author);
            /*! \fn
              * Constructor with required parameters
              * @param stream_block
              */
            PdbAuthorSection(std::stringstream& stream_block);

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
/** \addtogroup Molecular_Data_Structure
               * @{
               */
            /*! \fn
              * An accessor function in order to access to the record name in a author card
              * @return record_name_ attribute of the current object of this class
              */
            std::string GetRecordName();
            /*! \fn
              * An accessor function in order to access to the author in a author card
              * @return author_ attribute of the current object of this class
              */
            std::string GetAuthor();
/** @}*/
            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
/** \addtogroup Manipulators
               * @{
               */
            /*! \fn
              * A mutator function in order to set the record name of the current object
              * Set the record_name_ attribute of the current author card
              * @param record_name The record name of the current object
              */
            void SetRecordName(const std::string record_name);
            /*! \fn
              * A mutator function in order to set the author of the current object
              * Set the author_ attribute of the current author card
              * @param author The author of the current object
              */
            void SetAuthor(const std::string author);
/** @}*/
            //////////////////////////////////////////////////////////
            //                        FUNCTIONS                     //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the pdb author card contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cerr);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            std::string record_name_;           /*!< Record name which appears in the first column of each line in a pdb file >*/
            std::string author_;              /*!< Author that appears in KEYWORD record of a pdb file >*/
    };
}

#endif // PDBAUTHORSECTION_HPP
