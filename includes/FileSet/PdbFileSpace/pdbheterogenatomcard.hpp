// Created by: Delaram Rahbarinia
// Modified by: Alireza Khatamian, Delaram Rahbarinia

#ifndef PDBHETEROGENATOMCARD_HPP
#define PDBHETEROGENATOMCARD_HPP

#include <string>
#include <map>
#include <iostream>

namespace PdbFileSpace
{
    class PdbAtom;

    class PdbHeterogenAtomCard
    {
        public:
            //////////////////////////////////////////////////////////
            //                       TYPE DEFINITION                //
            //////////////////////////////////////////////////////////
            typedef PdbAtom PdbHeterogenAtom;
            typedef std::map<int, PdbHeterogenAtom*> PdbHeterogenAtomMap;            

            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            PdbHeterogenAtomCard();
            /*! \fn
              * Constructor with required parameters
              * @param stream_block
              * @param index
              */
            PdbHeterogenAtomCard(std::stringstream& stream_block, std::string index);

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
            /*! \fn
              * An accessor function in order to access to the record name in a heterogen atom card
              * @return record_name_ attribute of the current object of this class
              */
            std::string GetRecordName();
            /*! \fn
              * An accessor function in order to access to the heterogen atoms in a heterogen atom card
              * @return heterogen_atoms_ attribute of the current object of this class
              */
            PdbHeterogenAtomMap GetHeterogenAtoms();

            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A mutator function in order to set the record name of the current object
              * Set the record_name_ attribute of the current heterogen atom card
              * @param record_name The record name of the current object
              */
            void SetRecordName(const std::string record_name);

            /*! \fn
              * A mutator function in order to set the atoms of the current object
              * Set the heterogen_atoms_ attribute of the current atom card
              * @param heterogen_atoms Heterogen atoms attribute of the current object
              */
            void SetHeterogenAtoms(PdbHeterogenAtomMap heterogen_atoms);

            //////////////////////////////////////////////////////////
            //                        FUNCTIONS                     //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the heterogen name card contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cout);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            std::string record_name_;
            PdbHeterogenAtomMap heterogen_atoms_;

    };
}

#endif // PDBHETEROGENATOMCARD_HPP
