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
            PdbConnectCard();
            PdbConnectCard(std::stringstream& stream_block);

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
            std::string GetRecordName();
            BondedAtomsSerialNumbersMap GetBondedAtomsSerialNumbers();

            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
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
