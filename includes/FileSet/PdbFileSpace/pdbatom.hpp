#ifndef PDBATOM_HPP
#define PDBATOM_HPP

#include <string>
#include "../../../includes/Geometry/coordinate.hpp"

namespace PdbFileSpace
{
    class PdbAtom
    {
        public:

            //////////////////////////////////////////////////////////
            //                       Constructor                    //
            //////////////////////////////////////////////////////////
            PdbAtom();

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
            int GetAtomSerialNumber();
            std::string GetAtomName();
            char GetAtomAlternateLocation();
            std::string GetAtomResidueName();
            char GetAtomChainId();
            int GetAtomResidueSequenceNumber();
            char GetAtomInsertionCode();
            Geometry::Coordinate GetAtomOrthogonalCoordinate();
            double GetAtomOccupancy();
            double GetAtomTempretureFactor();
            std::string GetAtomElementSymbol();
            std::string GetAtomCharge();

            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
            void SetAtomSerialNumber(int atom_serial_number);
            void SetAtomName(const std::string atom_name);
            void SetAtomAlternateLocation(char atom_alternate_location);
            void SetAtomResidueName(const std::string atom_residue_name);
            void SetAtomChainId(char atom_chain_id);
            void SetAtomResidueSequenceNumber(int atom_residue_sequence_number);
            void SetAtomInsertionCode(char atom_insertion_code);
            void SetAtomOrthogonalCoordinate(Geometry::Coordinate atom_orthogonal_coordinate);
            void SetAtomOccupancy(double atom_occupancy);
            void SetAtomTempretureFactor(double atom_tempreture_factor);
            void SetAtomElementSymbol(const std::string atom_element_symbol);
            void SetAtomCharge(const std::string atom_charge);

            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////



        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            int atom_serial_number_;
            std::string atom_name_;
            char atom_alternate_location_;
            std::string atom_residue_name_;
            char atom_chain_id_;
            int atom_residue_sequence_number_;
            char atom_insertion_code_;
            Geometry::Coordinate atom_orthogonal_coordinate_;
            double atom_occupancy_;
            double atom_tempreture_factor_;
            std::string atom_element_symbol_;
            std::string atom_charge_;

    };
}

#endif // PDBATOM_HPP
