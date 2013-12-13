#ifndef PDBATOM_HPP
#define PDBATOM_HPP

#include <string>

namespace Geometry
{
    class Coordinate;
}

namespace PdbFileSpace
{
    class PdbAtom
    {
        public:

            //////////////////////////////////////////////////////////
            //                       Constructor                    //
            //////////////////////////////////////////////////////////
            PdbAtomCard();

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
            int GetAtomSerialNumber();
            std::string GetName();
            char GetAtomAlternateLocation();
            std::string GetAtomResidueName();
            char GetAtomChainId();
            int GetAtomResidueSequenceNumber();
            char GetAtomInsertionCode();
            Geometry::Coordinate GetAtomOrthogonalCoordinate();
            double GetAtomOccupancy();
            double GetAtomTempretureFactor();
            std::string GetAtomElementFactor();

            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
            void SetRecordName(const std::string record_name);

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

#endif // PDBATOM_HPP
