#ifndef PDBCONNECTCARD_H
#define PDBCONNECTCARD_H

#include <string>
#include <map>
#include <vector>

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

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
            std::string GetRecordName();
            BondedAtomsSerialNumbersMap GetBondedAtomsSerialNumbers();

            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
            void SetRecordName(std::string record_name);

            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////



        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            std::string record_name_;
            BondedAtomsSerialNumbersMap bonded_atom_serial_numbers_;

    };
}

#endif // PDBCONNECTCARD_H
