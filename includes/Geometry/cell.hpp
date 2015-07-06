#ifndef CELL_HPP
#define CELL_HPP

#include <iostream>
#include "coordinate.hpp"

namespace Geometry
{
    class Cell
    {
        public:
            //////////////////////////////////////////////////////////
            //                       Constructor                    //
            //////////////////////////////////////////////////////////
            Cell();
            Cell(Geometry::Coordinate min, Geometry::Coordinate max);
            Cell(Geometry::Coordinate min, Geometry::Coordinate max, double charge, double potential_energy);

            //////////////////////////////////////////////////////////
            //                           ACCESSOR                   //
            //////////////////////////////////////////////////////////
            Geometry::Coordinate GetMinCorner();
            Geometry::Coordinate GetMaxCorner();
            double GetCellCharge();
            double GetCellPotentialEnergy();
            double GetCellLength();
            double GetCellWidth();
            double GetCellHeight();

            //////////////////////////////////////////////////////////
            //                           MUTATOR                    //
            //////////////////////////////////////////////////////////
            void SetMinCorner(Geometry::Coordinate min);
            void SetMaxCorner(Geometry::Coordinate max);
            void SetCellCharge(double charge);
            void SetCellPotentialEnergy(double potential_energy);

            //////////////////////////////////////////////////////////+
            //                         FUNCTIONS                    //
            //////////////////////////////////////////////////////////
            Geometry::Coordinate GetCellCenter();

            //////////////////////////////////////////////////////////
            //                     DISPLAY FUNCTIONS                //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out a cell
              * Print out the current cell in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cout);

    private:
            //////////////////////////////////////////////////////////
            //                         ATTRIBUTES                   //
            //////////////////////////////////////////////////////////
            Geometry::Coordinate min_corner_;
            Geometry::Coordinate max_corner_;
            double cell_charge_;
            double cell_potential_energy_;
    };
}
#endif // CELL_HPP
