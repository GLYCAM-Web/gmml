// Author: Alireza Khatamian

#ifndef PDBSHEET_HPP
#define PDBSHEET_HPP

#include <string>
#include <vector>
#include<sstream>

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
            PdbSheet();
            PdbSheet(const std::string& sheet_id, int number_of_strands, const SheetStrandVector strands);
            PdbSheet(std::stringstream& sheet_block);

            //////////////////////////////////////////////////////////
            //                         ACCESSOR                     //
            //////////////////////////////////////////////////////////
            std::string GetSheetId();
            int GetNumberOfStrands();
            SheetStrandVector GetStrands();

            //////////////////////////////////////////////////////////
            //                          MUTATOR                     //
            //////////////////////////////////////////////////////////
            void SetSheetId(const std::string sheet_id);
            void SetNumberOfStrands(int number_of_strands);
            void SetStrands(const SheetStrandVector strands);
            void AddStrand(PdbSheetStrand* strand);

            //////////////////////////////////////////////////////////
            //                        FUNCTIONS                     //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                      DISPLAY FUNCTION                //
            //////////////////////////////////////////////////////////

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
