// Created by: Delaram Rahbarinia
// Modified by: Alireza Khatamian, Delaram Rahbarinia

#ifndef PDBCONNECTCARD_HPP
#define PDBCONNECTCARD_HPP

#include <string>
#include <map>
#include <vector>
#include <sstream>
#include <iostream>

namespace PdbFileSpace
{
    class PdbConnectCard
    {
        public:

            //////////////////////////////////////////////////////////
            //                       TYPE DEFINITION                //
            //////////////////////////////////////////////////////////
            /*! \typedef
              * A mapping between a specific atom serial number and its bonded atoms serial numbers in a connect card of a pdb file
              */
            typedef std::map<int, std::vector<int> > BondedAtomsSerialNumbersMap;

            //////////////////////////////////////////////////////////
            //                       Constructor                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            PdbConnectCard();
            /*! \fn
              * A constructor that get a stream block of connect card and parse the whole block to fill the related fields
              * @param stream_block A whole block of connect card in a pdb file
              */
            PdbConnectCard(std::stringstream& stream_block);

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
            /*! \fn
              * An accessor function in order to access to the record name in a connect card
              * @return record_name_ attribute of the current object of this class
              */
            std::string GetRecordName();            
            /*! \fn
              * An accessor function in order to access to the bonded atoms serial numbers in a connect card
              * @return bonded_atoms_serial_numbers_ attribute of the current object of this class
              */
            BondedAtomsSerialNumbersMap GetBondedAtomsSerialNumbers();

            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A mutator function in order to set the record name of the current object
              * Set the record_name_ attribute of the current connect card
              * @param record_name The atom serial number of the current object
              */
            void SetRecordName(const std::string record_name);            

            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the connect card contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cout);


        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            std::string record_name_;                                   /*!< Record name of connect card in a pdb file: "CONECT" */
            BondedAtomsSerialNumbersMap bonded_atom_serial_numbers_;    /*!< Map of bonded atoms to a specific atom by its serial number */

    };
}

#endif // PDBCONNECTCARD_HPP
