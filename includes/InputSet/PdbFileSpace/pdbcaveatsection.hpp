// Created by: Dave Montgomery

#ifndef PDBCAVEATSECTION_HPP
#define PDBCAVEATSECTION_HPP

#include <string>
#include <sstream>
#include <iostream>

namespace PdbFileSpace
{
    class PdbCaveatSection
    {
        public:
            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            PdbCaveatSection();
            /*! \fn
              * Constructor with required parameters
              * @param record_name Name for a caveat card record which appears in the first column of each line in a pdb file
              * @param caveat Caveat of a pdb file
              */
            PdbCaveatSection(const std::string& record_name, const std::string& caveat);
            /*! \fn
              * Constructor with required parameters
              * @param stream_block
              */
            PdbCaveatSection(std::stringstream& stream_block);

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
/** \addtogroup Molecular_Data_Structure
               * @{
               */
            /*! \fn
              * An accessor function in order to access to the record name in a caveat card
              * @return record_name_ attribute of the current object of this class
              */
            std::string GetRecordName();
            /*! \fn
              * An accessor function in order to access to the caveat in a caveat card
              * @return caveat_ attribute of the current object of this class
              */
            std::string GetCaveat();
/** @}*/
            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
/** \addtogroup Manipulators
               * @{
               */
            /*! \fn
              * A mutator function in order to set the record name of the current object
              * Set the record_name_ attribute of the current caveat card
              * @param record_name The record name of the current object
              */
            void SetRecordName(const std::string record_name);
            /*! \fn
              * A mutator function in order to set the caveat of the current object
              * Set the caveat_ attribute of the current caveat card
              * @param caveat The caveat of the current object
              */
            void SetCaveat(const std::string caveat);
/** @}*/
            //////////////////////////////////////////////////////////
            //                        FUNCTIONS                     //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the pdb caveat card contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cerr);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            std::string record_name_;           /*!< Record name which appears in the first column of each line in a pdb file >*/
            std::string caveat_;                 /*!< Caveat that appears in CAVEAT record of a pdb file >*/
    };
}

#endif // PDBCAVEATSECTION_HPP
