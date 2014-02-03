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
            typedef std::map<int, std::vector<int> > BondedAtomsSerialNumbersMap;

            //////////////////////////////////////////////////////////
            //                       Constructor                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            PdbConnectCard();
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
            void Print(std::ostream& out = std::cout);


        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            std::string record_name_;
            BondedAtomsSerialNumbersMap bonded_atom_serial_numbers_;

    };
}

#endif // PDBCONNECTCARD_HPP
