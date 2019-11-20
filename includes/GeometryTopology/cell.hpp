#ifndef CELL_HPP
#define CELL_HPP

#include <iostream>

#include "../../includes/common.hpp"
#include "coordinate.hpp"
#include "grid.hpp"

namespace GeometryTopology
{
    class Cell
    {
        public:
            //////////////////////////////////////////////////////////
            //                       Constructor                    //
            //////////////////////////////////////////////////////////
            Cell();
            Cell(GeometryTopology::Coordinate* min, GeometryTopology::Coordinate* max);
            Cell(Grid* grid, GeometryTopology::Coordinate* min, GeometryTopology::Coordinate* max);
            Cell(GeometryTopology::Coordinate* min, GeometryTopology::Coordinate* max, double charge, double potential_energy);
            Cell(Grid* grid, GeometryTopology::Coordinate* min, GeometryTopology::Coordinate* max, double charge, double potential_energy);

            //////////////////////////////////////////////////////////
            //                           ACCESSOR                   //
            //////////////////////////////////////////////////////////
/** \addtogroup Molecular_Data_Structure
               * @{
               */
            GeometryTopology::Coordinate* GetMinCorner();
            GeometryTopology::Coordinate* GetMaxCorner();
            double GetCellCharge();
            double GetCellPotentialEnergy();
            double GetCellLength();
            double GetCellWidth();
            double GetCellHeight();
            Grid* GetGrid();
/** @}*/
            //////////////////////////////////////////////////////////
            //                           MUTATOR                    //
            //////////////////////////////////////////////////////////
/** \addtogroup Manipulators
               * @{
               */
            void SetMinCorner(GeometryTopology::Coordinate* min);
            void SetMaxCorner(GeometryTopology::Coordinate* max);
            void SetCellCharge(double charge);
            void SetCellPotentialEnergy(double potential_energy);
            void SetGrid(Grid* grid);
/** @}*/
            //////////////////////////////////////////////////////////+
            //                         FUNCTIONS                    //
            //////////////////////////////////////////////////////////
/** \addtogroup Molecular_Data_Structure
               * @{
               */
            GeometryTopology::Coordinate* GetCellCenter();
            void CalculateCellCharge();
            void CalculateCellPotentialEnergy(double ion_radius);
            void CalculateBoxCharge();
            void CalculateBoxPotentialEnergy();
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
            GeometryTopology::Coordinate* min_corner_;
            GeometryTopology::Coordinate* max_corner_;
            double cell_charge_;
            double cell_potential_energy_;
            Grid* grid_;
    };
}

#endif // CELL_HPP
