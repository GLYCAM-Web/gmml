// Created by: Alireza Khatamian
// Modified by: Alireza Khatamian, Delaram Rahbarinia, Dave Montgomery
#ifndef PDBSHEETCARD_HPP
#define PDBSHEETCARD_HPP

#include <string>
#include <vector>
#include <sstream>
#include <iostream>

namespace PdbFileSpace
{
    class PdbSheetStrand;
    class PdbSheetCard
    {
        public:
            //////////////////////////////////////////////////////////
            //                    TYPE DEFINITION                   //
            //////////////////////////////////////////////////////////
            /*! \typedef
              * Vector of sheet strands
              */
            typedef std::vector<PdbSheetStrand*> SheetStrandVector;

            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            PdbSheetCard();
            /*! \fn
              * Constructor with required parameters
              * @param sheet_id Sheet identifier of a sheet in a pdb file
              * @param number_of_strands Number of strands that are involved in a sheet in a pdb file
              * @param strands Vector of involving strands in the current sheet
              */
            PdbSheetCard(const std::string& sheet_id, int number_of_strands, const SheetStrandVector strands);
            /*! \fn
              * Constructor with required parameters
              * @param sheet_block
              */
            PdbSheetCard(std::stringstream& sheet_block);

            //////////////////////////////////////////////////////////
            //                         ACCESSOR                     //
            //////////////////////////////////////////////////////////
/** \addtogroup Molecular_Data_Structure
               * @{
               */
            /*! \fn
              * An accessor function in order to access to the sheet id in a pdb sheet
              * @return sheet_id_ attribute of the current object of this class
              */
            std::string GetSheetId();
            /*! \fn
              * An accessor function in order to access to the number of strands in a pdb sheet
              * @return number_of_strands_ attribute of the current object of this class
              */
            int GetNumberOfStrands();
            /*! \fn
              * An accessor function in order to access to the strands in a pdb sheet
              * @return strands_ attribute of the current object of this class
              */
            SheetStrandVector GetStrands();
/** @}*/
            //////////////////////////////////////////////////////////
            //                          MUTATOR                     //
            //////////////////////////////////////////////////////////
/** \addtogroup Manipulators
               * @{
               */
            /*! \fn
              * A mutator function in order to set the sheet id of the current object
              * Set the sheet_id_ attribute of the current sheet
              * @param sheet_id The sheet id of the current object
              */
            void SetSheetId(const std::string sheet_id);
            /*! \fn
              * A mutator function in order to set the number of strands of the current object
              * Set the number_of_strands_ attribute of the current sheet
              * @param number_of_strands The number of strands of the current object
              */
            void SetNumberOfStrands(int number_of_strands);
            /*! \fn
              * A mutator function in order to set the strands of the current object
              * Set the strands_ attribute of the current sheet
              * @param strands The strands of the current object
              */
            void SetStrands(const SheetStrandVector strands);
            /*! \fn
              * A function in order to add the strand to the current object
              * Set the strand_ attribute of the current sheet
              * @param strand The strand of the current object
              */
            void AddStrand(PdbSheetStrand* strand);
/** @}*/
            //////////////////////////////////////////////////////////
            //                        FUNCTIONS                     //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                      DISPLAY FUNCTION                //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the pdb sheet contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cerr);

        private:
            //////////////////////////////////////////////////////////
            //                        ATTRIBUTES                    //
            //////////////////////////////////////////////////////////
            std::string sheet_id_;          /*!< Sheet identifier >*/
            int number_of_strands_;         /*!< Number of strands involving in a sheet >*/
            SheetStrandVector strands_;     /*!< Strands involving in a sheet >*/
    };
}

#endif // PDBSHEETCARD_HPP
