// Created by: Alireza Khatamian
// Modified by: Alireza Khatamian, Delaram Rahbarinia

#ifndef PDBSHEETCARD_HPP
#define PDBSHEETCARD_HPP

#include <string>
#include <map>
#include <sstream>
#include <iostream>

namespace PdbFileSpace
{
    class PdbSheet;
    class PdbSheetCard
    {
        public:
            //////////////////////////////////////////////////////////
            //                    TYPE DEFINITION                   //
            //////////////////////////////////////////////////////////
            /*! \typedef
              * Mapping between sheet id and the sheet itself
              */
            typedef std::map<std::string, PdbSheet*> SheetMap;

            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            PdbSheetCard();
            /*! \fn
              * Constructor with required parameters
              * @param record_name
              */
            PdbSheetCard(const std::string& record_name);
            /*! \fn
              * Constructor with required parameters
              * @param stream_block
              */
            PdbSheetCard(std::stringstream& stream_block);

            //////////////////////////////////////////////////////////
            //                         ACCESSOR                     //
            //////////////////////////////////////////////////////////
/** \addtogroup Molecular_Data_Structure
               * @{
               */
            /*! \fn
              * An accessor function in order to access to the record name in a sheet card
              * @return record_name_ attribute of the current object of this class
              */
            std::string GetRecordName();
            /*! \fn
              * An accessor function in order to access to the sheets in a sheet card
              * @return sheets_ attribute of the current object of this class
              */
            SheetMap GetSheets();
/** @}*/
            //////////////////////////////////////////////////////////
            //                          MUTATOR                     //
            //////////////////////////////////////////////////////////
/** \addtogroup Manipulators
               * @{
               */
            /*! \fn
              * A mutator function in order to set the record name of the current object
              * Set the record_name_ attribute of the current sheet card
              * @param record_name The record name of the current object
              */
            void SetRecordName(const std::string record_name);
/** @}*/
            //////////////////////////////////////////////////////////
            //                        FUNCTIONS                     //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                      DISPLAY FUNCTION                //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the pdb sheet card contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cout);

        private:
            //////////////////////////////////////////////////////////
            //                        ATTRIBUTES                    //
            //////////////////////////////////////////////////////////
            std::string record_name_;       /*!< Name of the sheet card record that appears in the first column of each line of a pdb file >*/
            SheetMap sheets_;               /*!< All sheets that are in a sheet card of a pdb file in a map data structure >*/
    };
}

#endif // PDBSHEETCARD_HPP
