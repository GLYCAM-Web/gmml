#ifndef GMML_INCLUDES_INPUTSET_PDBFILE_CONECTRECORD_HPP
#define GMML_INCLUDES_INPUTSET_PDBFILE_CONECTRECORD_HPP

#include <string>
#include <iostream>

#include "includes/InputSet/PdbFile/atomRecord.hpp"

namespace pdb
{
    class ConectRecord
    {
        // This is passed in as serial numbers, but they can/will change so store pointers to the atoms that are bonded instead.
        // When asked to print or output to file, get the serial numbers via the pointers.
        public:
            //////////////////////////////////////////////////////////
            //                    CONSTRUCTOR                       //
            //////////////////////////////////////////////////////////
            ConectRecord(std::string& line, int modelNumber = 1);
            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            void Print(std::ostream& out = std::cerr);
        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            std::pair<AtomRecord*, AtomRecord*> connectRecord_;


    };
}

#endif// GMML_INCLUDES_INPUTSET_PDBFILE_CONECTRECORD_HPP
