#ifndef GMML_INCLUDES_INPUTSET_PDBFILE_CONECTRECORDS_HPP
#define GMML_INCLUDES_INPUTSET_PDBFILE_CONECTRECORDS_HPP

#include <string>
#include <iostream>

#include "includes/InputSet/PdbFile/atomRecord.hpp"
#include "pdbAssembly.hpp"

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
            ConectRecord(std::string &line, PdbAssembly& coordinateSection);
            ConectRecord(std::vector<AtomRecord*> atoms );
            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            void Print(std::ostream& out = std::cerr) const;
            void Write(std::ostream& stream) const;
        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            std::vector<AtomRecord*> atomRecordPtrs_;
    };
}

#endif// GMML_INCLUDES_INPUTSET_PDBFILE_CONECTRECORDS_HPP
