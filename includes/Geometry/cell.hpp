#ifndef CELL_HPP
#define CELL_HPP

#include <ostream>
#include "../MolecularModeling/assembly.hpp"

namespace Geometry
{
    class Coordinate;

    class Cell
    {
        public:
            //////////////////////////////////////////////////////////
            //                       Constructor                    //
            //////////////////////////////////////////////////////////
            Cell();
            Cell(Coordinate min, Coordinate max);
            Cell(Coordinate min, Coordinate max, double charge, double potential_energy);

            //////////////////////////////////////////////////////////
            //                           ACCESSOR                   //
            //////////////////////////////////////////////////////////
            Coordinate GetMinCorner();
            Coordinate GetMaxCorner();
            double GetCellCharge();
            double GetCellPotentialEnergy();
            double GetCellLength();
            double GetCellWidth();
            double GetCellHeight();

            //////////////////////////////////////////////////////////
            //                           MUTATOR                    //
            //////////////////////////////////////////////////////////
            void SetMinCorner(Coordinate min);
            void SetMaxCorner(Coordinate max);
            void SetCellCharge(double charge);
            void SetCellPotentialEnergy(double potential_energy);

            //////////////////////////////////////////////////////////+
            //                         FUNCTIONS                    //
            //////////////////////////////////////////////////////////
            Coordinate GetCellCenter();
            void CalculateCellCharge(MolecularModeling::Assembly* assembly);
            void CalculateCellPotentialEnergy(MolecularModeling::Assembly* assembly, double radius);
            void CalculateBoxCharge(Assembly* assembly);
            void CalculateBoxPotentialEnergy(Assembly* assembly);

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
            Coordinate min_corner_;
            Coordinate max_corner_;
            double cell_charge_;
            double cell_potential_energy_;
    };
}
#endif // CELL_HPP
