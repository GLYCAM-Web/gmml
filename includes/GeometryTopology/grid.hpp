#ifndef GRID_HPP
#define GRID_HPP

#include <vector>
#include "../MolecularModeling/assembly.hpp"

namespace GeometryTopology
{
    // class Coordinate;
    class Cell;
    class Grid
    {
        public:
            //////////////////////////////////////////////////////////
            //                    TYPE DEFINITION                   //
            //////////////////////////////////////////////////////////
            typedef std::vector<GeometryTopology::Cell*> CellVector;

            //////////////////////////////////////////////////////////
            //                       Constructor                    //
            //////////////////////////////////////////////////////////
            Grid();
           /*! \todo Figure out if this function really should use ion_radius.  */
            Grid(MolecularModeling::Assembly* assembly, Coordinate* min, Coordinate* max, double ion_radius, double ion_charge);
            Grid(MolecularModeling::Assembly* assembly, Coordinate* min, Coordinate* max, double cell_length, double cell_width, double cell_height);
            Grid(Grid &grid);

            //////////////////////////////////////////////////////////
            //                           ACCESSOR                   //
            //////////////////////////////////////////////////////////
/** \addtogroup Molecular_Data_Structure
               * @{
               */
            Coordinate* GetMinCorner();
            Coordinate* GetMaxCorner();
            CellVector GetCells();
            MolecularModeling::Assembly* GetAssembly();
/** @}*/
            //////////////////////////////////////////////////////////
            //                           MUTATOR                    //
            //////////////////////////////////////////////////////////
/** \addtogroup Manipulators
               * @{
               */
            void SetMinCorner(Coordinate* min);
            void SetMaxCorner(Coordinate* max);
            void SetCells(CellVector cells);
            void SetAssembly(MolecularModeling::Assembly* assembly);

            //////////////////////////////////////////////////////////
            //                         FUNCTIONS                    //
            //////////////////////////////////////////////////////////
            void UpdateGrid(double lenght, double width, double height);
            void UpdateGrid(double ion_charge);
            void CalculateCellsCharge();
            void CalculateCellsPotentialEnergy(double ion_radius);
            void CalculateBoxCharge();
            void CalculateBoxPotentialEnergy();
            GeometryTopology::Cell* GetBestBox(Grid* grid, double ion_charge);
            std::vector<Coordinate*> GetBestPositions(double ion_charge);
/** @}*/
            //////////////////////////////////////////////////////////
            //                     DISPLAY FUNCTIONS                //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out a cell
              * Print out the current cell in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cerr);

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
