#ifndef GRID_HPP
#define GRID_HPP

#include <vector>
#include "../MolecularModeling/assembly.hpp"

namespace GeometryTopology
{
    class Coordinate;
    class Cell;
    class Grid
    {
        public:
            //////////////////////////////////////////////////////////
            //                    TYPE DEFINITION                   //
            //////////////////////////////////////////////////////////
            typedef std::vector<Cell*> CellVector;

            //////////////////////////////////////////////////////////
            //                       Constructor                    //
            //////////////////////////////////////////////////////////
            Grid();
            Grid(MolecularModeling::Assembly* assembly, Coordinate* min, Coordinate* max, double ion_radius, double ion_charge);

            //////////////////////////////////////////////////////////
            //                           ACCESSOR                   //
            //////////////////////////////////////////////////////////
            Coordinate* GetMinCorner();
            Coordinate* GetMaxCorner();
            CellVector GetCells();
            MolecularModeling::Assembly* GetAssembly();

            //////////////////////////////////////////////////////////
            //                           MUTATOR                    //
            //////////////////////////////////////////////////////////
            void SetMinCorner(Coordinate* min);
            void SetMaxCorner(Coordinate* max);
            void SetCells(CellVector cells);
            void SetAssembly(MolecularModeling::Assembly* assembly);

            //////////////////////////////////////////////////////////
            //                         FUNCTIONS                    //
            //////////////////////////////////////////////////////////
            void UpdateGrid(double ion_charge);
            void CalculateCellsCharge();
            void CalculateCellsPotentialEnergy(double ion_radius);
            void CalculateBoxCharge();
            void CalculateBoxPotentialEnergy();
            Cell* GetBestBox(Grid* grid, double ion_charge);
            std::vector<Coordinate*> GetBestPositions(double ion_charge);

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
            Coordinate* min_corner_;
            Coordinate* max_corner_;
            CellVector cells_;
            MolecularModeling::Assembly* assembly_;

    };
}

#endif // GRID_HPP
