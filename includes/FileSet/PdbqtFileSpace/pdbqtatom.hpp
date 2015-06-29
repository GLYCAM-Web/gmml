#ifndef PDBQTATOM_HPP
#define PDBQTATOM_HPP

#include <string>
#include "../../Geometry/coordinate.hpp"

namespace PdbqtFileSpace
{
    class PdbqtAtom
    {
        public:

            //////////////////////////////////////////////////////////
            //                       Constructor                    //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////


            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////


        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            std::string type_;
            int serial_number_;
            std::string atom_name;
            char atom_alternate_location_;
            std::string atom_residue_name_;
            char atom_chain_id_;
            int atom_residue_sequence_number_;
            char atom_insertion_code_;
            Geometry::Coordinate atom_orthogonal_coordinate_;
            double atom_occupancy_;
            double atom_temperature_factor_;
            double atom_charge_;
            std::string atom_type;

    };
}

#endif // PDBQTATOM_HPP
