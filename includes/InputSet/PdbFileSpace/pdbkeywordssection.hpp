// Created by: Dave Montgomery

#ifndef PDBKEYWORDSSECTION_HPP
#define PDBKEYWORDSSECTION_HPP

#include <string>
#include <sstream>
#include <iostream>

namespace PdbFileSpace
{
    class PdbKeywordsSection
    {
        public:
            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            PdbKeywordsSection();
            /*! \fn
              * Constructor with required parameters
              * @param record_name Name for a keywords card record which appears in the first column of each line in a pdb file
              * @param keywords Keywords of a pdb file
              */
            PdbKeywordsSection(const std::string& record_name, const std::string& keywords);
            /*! \fn
              * Constructor with required parameters
              * @param stream_block
              */
            PdbKeywordsSection(std::stringstream& stream_block);

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
/** \addtogroup Molecular_Data_Structure
               * @{
               */
            /*! \fn
              * An accessor function in order to access to the record name in a keywords card
              * @return record_name_ attribute of the current object of this class
              */
            std::string GetRecordName();
            /*! \fn
              * An accessor function in order to access to the keywords in a keywords card
              * @return keywords_ attribute of the current object of this class
              */
            std::string GetKeywords();
/** @}*/
            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
/** \addtogroup Manipulators
               * @{
               */
            /*! \fn
              * A mutator function in order to set the record name of the current object
              * Set the record_name_ attribute of the current keywords card
              * @param record_name The record name of the current object
              */
            void SetRecordName(const std::string record_name);
            /*! \fn
              * A mutator function in order to set the keywords of the current object
              * Set the keywords_ attribute of the current keywords card
              * @param keywords The keywords of the current object
              */
            void SetKeywords(const std::string keywords);
/** @}*/
            //////////////////////////////////////////////////////////
            //                        FUNCTIONS                     //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the pdb keywords card contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cerr);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            std::string record_name_;           /*!< Record name which appears in the first column of each line in a pdb file >*/
            std::string keywords_;              /*!< Keywords that appears in KEYWORD record of a pdb file >*/
    };
}

#endif // PDBKEYWORDSSECTION_HPP
