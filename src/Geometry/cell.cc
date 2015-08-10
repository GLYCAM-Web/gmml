#include "../../includes/common.hpp"
#include "../../includes/Geometry/cell.hpp"
#include "../../includes/Geometry/coordinate.hpp"
#include "../../includes/MolecularModeling/assembly.hpp"

using namespace gmml;
using namespace Geometry;
using namespace MolecularModeling;
//////////////////////////////////////////////////////////
//                       Constructor                    //
//////////////////////////////////////////////////////////
Cell::Cell()
{
    min_corner_ = new Coordinate();
    max_corner_ = new Coordinate();
    cell_charge_ = 0.0;
    cell_potential_energy_ = 0.0;
    grid_ = NULL;
}

Cell::Cell(Coordinate* min, Coordinate* max)
{
    min_corner_ = new Coordinate(min->GetX(), min->GetY(), min->GetZ());
    max_corner_ = new Coordinate(max->GetX(), max->GetY(), max->GetZ());
    cell_charge_ = 0.0;
    cell_potential_energy_ = 0.0;
    grid_ = NULL;
}

Cell::Cell(Grid *grid, Coordinate *min, Coordinate *max)
{
    min_corner_ = new Coordinate(min->GetX(), min->GetY(), min->GetZ());
    max_corner_ = new Coordinate(max->GetX(), max->GetY(), max->GetZ());
    cell_charge_ = 0.0;
    cell_potential_energy_ = 0.0;
    grid_ = grid;
}

Cell::Cell(Coordinate* min, Coordinate* max, double charge, double potential_energy)
{
    min_corner_ = new Coordinate(min->GetX(), min->GetY(), min->GetZ());
    max_corner_ = new Coordinate(max->GetX(), max->GetY(), max->GetZ());
    cell_charge_ = charge;
    cell_potential_energy_ = potential_energy;
    grid_ = NULL;
}

Cell::Cell(Grid *grid, Coordinate *min, Coordinate *max, double charge, double potential_energy)
{
    min_corner_ = new Coordinate(min->GetX(), min->GetY(), min->GetZ());
    max_corner_ = new Coordinate(max->GetX(), max->GetY(), max->GetZ());
    cell_charge_ = charge;
    cell_potential_energy_ = potential_energy;
    grid_ = grid;
}

//////////////////////////////////////////////////////////
//                           ACCESSOR                   //
//////////////////////////////////////////////////////////
Geometry::Coordinate* Cell::GetMinCorner()
{
    return min_corner_;
}

Geometry::Coordinate* Cell::GetMaxCorner()
{
    return max_corner_;
}

double Cell::GetCellCharge()
{
    return cell_charge_;
}

double Cell::GetCellPotentialEnergy()
{
    return cell_potential_energy_;
}

double Cell::GetCellLength()
{
    return max_corner_->GetX() - min_corner_->GetX();
}

double Cell::GetCellWidth()
{
    max_corner_->GetY() - min_corner_->GetY();
}

double Cell::GetCellHeight()
{
    max_corner_->GetZ() - min_corner_->GetZ();
}

Grid* Cell::GetGrid()
{
    return grid_;
}

//////////////////////////////////////////////////////////
//                           MUTATOR                    //
//////////////////////////////////////////////////////////
void Cell::SetMinCorner(Geometry::Coordinate* min)
{
    min_corner_->SetX(min->GetX());
    min_corner_->SetY(min->GetY());
    min_corner_->SetZ(min->GetZ());
}

void Cell::SetMaxCorner(Geometry::Coordinate* max)
{
    max_corner_->SetX(max->GetX());
    max_corner_->SetY(max->GetY());
    max_corner_->SetZ(max->GetZ());
}

void Cell::SetCellCharge(double charge)
{
    cell_charge_ = charge;
}

void Cell::SetCellPotentialEnergy(double potential_energy)
{
    cell_potential_energy_ = potential_energy;
}

void Cell::SetGrid(Grid *grid)
{
    grid_ = grid;
}

//////////////////////////////////////////////////////////+
//                         FUNCTIONS                    //
//////////////////////////////////////////////////////////
Geometry::Coordinate* Cell::GetCellCenter()
{
    return new Coordinate(min_corner_->GetX() + (max_corner_->GetX() - min_corner_->GetX()) / 2,
                          min_corner_->GetY() + (max_corner_->GetY() - min_corner_->GetY()) / 2,
                          min_corner_->GetZ() + (max_corner_->GetZ() - min_corner_->GetZ()) / 2);
}

void Cell::CalculateCellCharge()
{
    double charge = 0.0;
    Assembly::AtomVector all_atoms = this->grid_->GetAssembly()->GetAllAtomsOfAssembly();
    for(Assembly::AtomVector::iterator it = all_atoms.begin(); it != all_atoms.end(); it++)
    {
        Atom* atom = *it;
        if(atom->GetCoordinates().at(this->grid_->GetAssembly()->GetModelIndex())->GetX() <= this->GetMaxCorner()->GetX() &&
            atom->GetCoordinates().at(this->grid_->GetAssembly()->GetModelIndex())->GetY() <= this->GetMaxCorner()->GetY() &&
            atom->GetCoordinates().at(this->grid_->GetAssembly()->GetModelIndex())->GetZ() <= this->GetMaxCorner()->GetZ() &&
                atom->GetCoordinates().at(this->grid_->GetAssembly()->GetModelIndex())->GetX() > this->GetMinCorner()->GetX() &&
                atom->GetCoordinates().at(this->grid_->GetAssembly()->GetModelIndex())->GetY() > this->GetMinCorner()->GetY() &&
                atom->GetCoordinates().at(this->grid_->GetAssembly()->GetModelIndex())->GetZ() > this->GetMinCorner()->GetZ())
            charge += (atom->MolecularDynamicAtom::GetCharge() != dNotSet) ? atom->MolecularDynamicAtom::GetCharge() : 0.0;
    }
    cell_charge_ = charge;
}

void Cell::CalculateCellPotentialEnergy(double ion_radius)
{
    double potential_energy = 0.0;
    Assembly::AtomVector all_atoms = this->grid_->GetAssembly()->GetAllAtomsOfAssembly();
    Coordinate* center_of_cell = this->GetCellCenter();
    for(Assembly::AtomVector::iterator it = all_atoms.begin(); it != all_atoms.end(); it++)
    {
        Atom* atom = *it;
        double dist = sqrt((center_of_cell->GetX() - atom->GetCoordinates().at(this->grid_->GetAssembly()->GetModelIndex())->GetX()) *
                           (center_of_cell->GetX() - atom->GetCoordinates().at(this->grid_->GetAssembly()->GetModelIndex())->GetX()) +
                           (center_of_cell->GetY() - atom->GetCoordinates().at(this->grid_->GetAssembly()->GetModelIndex())->GetY()) *
                           (center_of_cell->GetY() - atom->GetCoordinates().at(this->grid_->GetAssembly()->GetModelIndex())->GetY()) +
                           (center_of_cell->GetZ() - atom->GetCoordinates().at(this->grid_->GetAssembly()->GetModelIndex())->GetZ()) *
                           (center_of_cell->GetZ() - atom->GetCoordinates().at(this->grid_->GetAssembly()->GetModelIndex())->GetZ()));

        double radius = (atom->MolecularDynamicAtom::GetRadius() != dNotSet) ? atom->MolecularDynamicAtom::GetRadius() : MINIMUM_RADIUS;
        if(dist < radius + GRID_OFFSET + ion_radius)
        {
            potential_energy = INFINITY;
            break;
        }
        else if(potential_energy == INFINITY)
            break;
        else
            potential_energy += ((atom->MolecularDynamicAtom::GetCharge() != dNotSet) ? atom->MolecularDynamicAtom::GetCharge() / dist : 0.0);
    }
    cell_potential_energy_ = potential_energy;
}

void Cell::CalculateBoxCharge()
{
    double charge = 0.0;
    Assembly::AtomVector all_atoms = this->grid_->GetAssembly()->GetAllAtomsOfAssembly();
    for(Assembly::AtomVector::iterator it = all_atoms.begin(); it != all_atoms.end(); it++)
    {
        Atom* atom = *it;
        if(atom->GetCoordinates().at(this->grid_->GetAssembly()->GetModelIndex())->GetX() <= this->GetMaxCorner()->GetX() &&
                atom->GetCoordinates().at(this->grid_->GetAssembly()->GetModelIndex())->GetY() <= this->GetMaxCorner()->GetY() &&
                atom->GetCoordinates().at(this->grid_->GetAssembly()->GetModelIndex())->GetZ() <= this->GetMaxCorner()->GetZ() &&
                atom->GetCoordinates().at(this->grid_->GetAssembly()->GetModelIndex())->GetX() > this->GetMinCorner()->GetX() &&
                atom->GetCoordinates().at(this->grid_->GetAssembly()->GetModelIndex())->GetY() > this->GetMinCorner()->GetY() &&
                atom->GetCoordinates().at(this->grid_->GetAssembly()->GetModelIndex())->GetZ() > this->GetMinCorner()->GetZ())
        charge += (atom->MolecularDynamicAtom::GetCharge() != dNotSet) ? atom->MolecularDynamicAtom::GetCharge() : 0.0;
    }
    cell_charge_ = charge;
}

void Cell::CalculateBoxPotentialEnergy()
{
    double potential_energy = 0.0;
    Assembly::AtomVector all_atoms = this->grid_->GetAssembly()->GetAllAtomsOfAssembly();
    Coordinate* center_of_cell = this->GetCellCenter();
    for(Assembly::AtomVector::iterator it = all_atoms.begin(); it != all_atoms.end(); it++)
    {
        Atom* atom = *it;
        double dist = sqrt((center_of_cell->GetX() - atom->GetCoordinates().at(this->grid_->GetAssembly()->GetModelIndex())->GetX()) *
                           (center_of_cell->GetX() - atom->GetCoordinates().at(this->grid_->GetAssembly()->GetModelIndex())->GetX()) +
                           (center_of_cell->GetY() - atom->GetCoordinates().at(this->grid_->GetAssembly()->GetModelIndex())->GetY()) *
                           (center_of_cell->GetY() - atom->GetCoordinates().at(this->grid_->GetAssembly()->GetModelIndex())->GetY()) +
                           (center_of_cell->GetZ() - atom->GetCoordinates().at(this->grid_->GetAssembly()->GetModelIndex())->GetZ()) *
                           (center_of_cell->GetZ() - atom->GetCoordinates().at(this->grid_->GetAssembly()->GetModelIndex())->GetZ()));
        if(dist == 0.0)
            dist = DIST_EPSILON;
        else
            potential_energy += ((atom->MolecularDynamicAtom::GetCharge() != dNotSet) ? atom->MolecularDynamicAtom::GetCharge() / dist : 0.0);
    }
    cell_potential_energy_ = potential_energy;
}

//////////////////////////////////////////////////////////
//                     DISPLAY FUNCTIONS                //
//////////////////////////////////////////////////////////
