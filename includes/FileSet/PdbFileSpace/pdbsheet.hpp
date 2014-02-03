// Created by: Alireza Khatamian
// Modified by: Alireza Khatamian, Delaram Rahbarinia

#ifndef PDBSHEET_HPP
#define PDBSHEET_HPP

#include <string>
#include <vector>
#include <sstream>
#include <iostream>

namespace PdbFileSpace
{
    class PdbSheetStrand;
    class PdbSheet
    {
        public:
            //////////////////////////////////////////////////////////
            //                    TYPE DEFINITION                   //
            //////////////////////////////////////////////////////////
            typedef std::vector<PdbSheetStrand*> SheetStrandVector;

            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            PdbSheet();
            /*! \fn
              * Constructor with required parameters
              * @param sheet_id
              * @param number_of_strands
              * @param strands
              */
            PdbSheet(const std::string& sheet_id, int number_of_strands, const SheetStrandVector strands);
            PdbSheet(std::stringstream& sheet_block);

            //////////////////////////////////////////////////////////
            //                         ACCESSOR                     //
            //////////////////////////////////////////////////////////
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

            //////////////////////////////////////////////////////////
            //                          MUTATOR                     //
            //////////////////////////////////////////////////////////
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

            //////////////////////////////////////////////////////////
            //                        FUNCTIONS                     //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                      DISPLAY FUNCTION                //
            //////////////////////////////////////////////////////////
            void Print(std::ostream& out = std::cout);

        private:
            //////////////////////////////////////////////////////////
            //                        ATTRIBUTES                    //
            //////////////////////////////////////////////////////////
            std::string sheet_id_;
            int number_of_strands_;
            SheetStrandVector strands_;
    };
}

#endif // PDBSHEET_HPP
