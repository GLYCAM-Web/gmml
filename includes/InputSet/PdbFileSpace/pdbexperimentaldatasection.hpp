// Created by: Dave Montgomery

#ifndef PDBEXPERIMENTALDATASECTION_HPP
#define PDBEXPERIMENTALDATASECTION_HPP

#include <string>
#include <sstream>
#include <iostream>

namespace PdbFileSpace
{
    class PdbExperimentalDataSection
    {
        public:
            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            PdbExperimentalDataSection();
            /*! \fn
              * Constructor with required parameters
              * @param record_name Name for a experimental_data card record which appears in the first column of each line in a pdb file
              * @param experimental_data ExperimentalData of a pdb file
              */
            PdbExperimentalDataSection(const std::string& record_name, const std::string& experimental_data);
            /*! \fn
              * Constructor with required parameters
              * @param stream_block
              */
            PdbExperimentalDataSection(std::stringstream& stream_block);

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
/** \addtogroup Molecular_Data_Structure
               * @{
               */
            /*! \fn
              * An accessor function in order to access to the record name in a experimental_data card
              * @return record_name_ attribute of the current object of this class
              */
            std::string GetRecordName();
            /*! \fn
              * An accessor function in order to access to the experimental_data in a experimental_data card
              * @return experimental_data_ attribute of the current object of this class
              */
            std::string GetExperimentalData();
/** @}*/
            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
/** \addtogroup Manipulators
               * @{
               */
            /*! \fn
              * A mutator function in order to set the record name of the current object
              * Set the record_name_ attribute of the current experimental_data card
              * @param record_name The record name of the current object
              */
            void SetRecordName(const std::string record_name);
            /*! \fn
              * A mutator function in order to set the experimental_data of the current object
              * Set the experimental_data_ attribute of the current experimental_data card
              * @param experimental_data The experimental_data of the current object
              */
            void SetExperimentalData(const std::string experimental_data);
/** @}*/
            //////////////////////////////////////////////////////////
            //                        FUNCTIONS                     //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the pdb experimental_data card contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cerr);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            std::string record_name_;           /*!< Record name which appears in the first column of each line in a pdb file >*/
            std::string experimental_data_;              /*!< ExperimentalData that appears in KEYWORD record of a pdb file >*/
    };
}

#endif // PDBEXPERIMENTALDATASECTION_HPP
