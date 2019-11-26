// Created by: Dave Montgomery

#ifndef PDBSPLITSECTION_HPP
#define PDBSPLITSECTION_HPP

#include <string>
#include <sstream>
#include <iostream>

namespace PdbFileSpace
{
    class PdbSplitSection
    {
        public:
            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            PdbSplitSection();
            /*! \fn
              * Constructor with required parameters
              * @param record_name Name for a split card record which appears in the first column of each line in a pdb file
              * @param split Split of a pdb file
              */
            PdbSplitSection(const std::string& record_name, const std::string& split);
            /*! \fn
              * Constructor with required parameters
              * @param stream_block
              */
            PdbSplitSection(std::stringstream& stream_block);

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
/** \addtogroup Molecular_Data_Structure
               * @{
               */
            /*! \fn
              * An accessor function in order to access to the record name in a split card
              * @return record_name_ attribute of the current object of this class
              */
            std::string GetRecordName();
            /*! \fn
              * An accessor function in order to access to the split in a split card
              * @return split_ attribute of the current object of this class
              */
            std::string GetSplit();
/** @}*/
            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
/** \addtogroup Manipulators
               * @{
               */
            /*! \fn
              * A mutator function in order to set the record name of the current object
              * Set the record_name_ attribute of the current split card
              * @param record_name The record name of the current object
              */
            void SetRecordName(const std::string record_name);
            /*! \fn
              * A mutator function in order to set the split of the current object
              * Set the split_ attribute of the current split card
              * @param split The split of the current object
              */
            void SetSplit(const std::string split);
/** @}*/
            //////////////////////////////////////////////////////////
            //                        FUNCTIONS                     //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the pdb split card contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cerr);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            std::string record_name_;           /*!< Record name which appears in the first column of each line in a pdb file >*/
            std::string split_;                 /*!< Split that appears in TITLE record of a pdb file >*/
    };
}

#endif // PDBSPLITSECTION_HPP
