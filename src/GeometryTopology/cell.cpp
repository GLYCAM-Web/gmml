#include "../../includes/GeometryTopology/cell.hpp"
#include "../../includes/GeometryTopology/coordinate.hpp"
#include "../../includes/MolecularModeling/assembly.hpp"

using GeometryTopology::Cell;	// Do we want to completely exclude "using" statements?

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

Cell::Cell(GeometryTopology::Coordinate* min, GeometryTopology::Coordinate* max)
{
    min_corner_ = new Coordinate(min->GetX(), min->GetY(), min->GetZ());
    max_corner_ = new Coordinate(max->GetX(), max->GetY(), max->GetZ());
    cell_charge_ = 0.0;
    cell_potential_energy_ = 0.0;
    grid_ = NULL;
}

Cell::Cell(GeometryTopology::Grid* grid, GeometryTopology::Coordinate *min, GeometryTopology::Coordinate *max)
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

Cell::Cell(GeometryTopology::Grid* grid, Coordinate *min, Coordinate *max, double charge, double potential_energy)
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
GeometryTopology::Coordinate* Cell::GetMinCorner()
{
    return min_corner_;
}

GeometryTopology::Coordinate* Cell::GetMaxCorner()
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
    return max_corner_->GetY() - min_corner_->GetY();
}

double Cell::GetCellHeight()
{
    return max_corner_->GetZ() - min_corner_->GetZ();
}

GeometryTopology::Grid* Cell::GetGrid()
{
    return grid_;
}

//////////////////////////////////////////////////////////
//                           MUTATOR                    //
//////////////////////////////////////////////////////////
void Cell::SetMinCorner(GeometryTopology::Coordinate* min)
{
    min_corner_->SetX(min->GetX());
    min_corner_->SetY(min->GetY());
    min_corner_->SetZ(min->GetZ());
}

void Cell::SetMaxCorner(GeometryTopology::Coordinate* max)
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

void Cell::SetGrid(GeometryTopology::Grid* grid)
{
    grid_ = grid;
}

//////////////////////////////////////////////////////////+
//                         FUNCTIONS                    //
//////////////////////////////////////////////////////////
GeometryTopology::Coordinate* Cell::GetCellCenter()
{
    return new Coordinate(min_corner_->GetX() + (max_corner_->GetX() - min_corner_->GetX()) / 2,
                          min_corner_->GetY() + (max_corner_->GetY() - min_corner_->GetY()) / 2,
                          min_corner_->GetZ() + (max_corner_->GetZ() - min_corner_->GetZ()) / 2);
}

void Cell::CalculateCellCharge()
{
    double charge = 0.0;
    MolecularModeling::AtomVector all_atoms = this->grid_->GetAssembly()->GetAllAtomsOfAssembly();
    for(MolecularModeling::AtomVector::iterator it = all_atoms.begin(); it != all_atoms.end(); it++)
    {
        MolecularModeling::Atom* atom = *it;
        if(atom->GetCoordinates().at(this->grid_->GetAssembly()->GetModelIndex())->GetX() <= this->GetMaxCorner()->GetX() &&
            atom->GetCoordinates().at(this->grid_->GetAssembly()->GetModelIndex())->GetY() <= this->GetMaxCorner()->GetY() &&
            atom->GetCoordinates().at(this->grid_->GetAssembly()->GetModelIndex())->GetZ() <= this->GetMaxCorner()->GetZ() &&
                atom->GetCoordinates().at(this->grid_->GetAssembly()->GetModelIndex())->GetX() > this->GetMinCorner()->GetX() &&
                atom->GetCoordinates().at(this->grid_->GetAssembly()->GetModelIndex())->GetY() > this->GetMinCorner()->GetY() &&
                atom->GetCoordinates().at(this->grid_->GetAssembly()->GetModelIndex())->GetZ() > this->GetMinCorner()->GetZ())
            charge += (atom->MolecularDynamicAtom::GetCharge() != gmml::dNotSet) ? atom->MolecularDynamicAtom::GetCharge() : 0.0;
    }
    cell_charge_ = charge;
}

void Cell::CalculateCellPotentialEnergy(double ion_radius)
{
    double potential_energy = 0.0;
    MolecularModeling::AtomVector all_atoms = this->grid_->GetAssembly()->GetAllAtomsOfAssembly();
    Coordinate* center_of_cell = this->GetCellCenter();
    for(MolecularModeling::AtomVector::iterator it = all_atoms.begin(); it != all_atoms.end(); it++)
    {
        MolecularModeling::Atom* atom = *it;
        double dist = sqrt((center_of_cell->GetX() - atom->GetCoordinates().at(this->grid_->GetAssembly()->GetModelIndex())->GetX()) *
                           (center_of_cell->GetX() - atom->GetCoordinates().at(this->grid_->GetAssembly()->GetModelIndex())->GetX()) +
                           (center_of_cell->GetY() - atom->GetCoordinates().at(this->grid_->GetAssembly()->GetModelIndex())->GetY()) *
                           (center_of_cell->GetY() - atom->GetCoordinates().at(this->grid_->GetAssembly()->GetModelIndex())->GetY()) +
                           (center_of_cell->GetZ() - atom->GetCoordinates().at(this->grid_->GetAssembly()->GetModelIndex())->GetZ()) *
                           (center_of_cell->GetZ() - atom->GetCoordinates().at(this->grid_->GetAssembly()->GetModelIndex())->GetZ()));

        double radius = (atom->MolecularDynamicAtom::GetRadius() != gmml::dNotSet) ? atom->MolecularDynamicAtom::GetRadius() : gmml::MINIMUM_RADIUS;
        if(dist < radius + gmml::GRID_OFFSET + ion_radius)
        {
            potential_energy = INFINITY;
            break;
        }
        else if(potential_energy == INFINITY)
            break;
        else
            potential_energy += ((atom->MolecularDynamicAtom::GetCharge() != gmml::dNotSet) ? atom->MolecularDynamicAtom::GetCharge() / dist : 0.0);
    }
    cell_potential_energy_ = potential_energy;
}

void Cell::CalculateBoxCharge()
{
    double charge = 0.0;
    MolecularModeling::AtomVector all_atoms = this->grid_->GetAssembly()->GetAllAtomsOfAssembly();
    for(MolecularModeling::AtomVector::iterator it = all_atoms.begin(); it != all_atoms.end(); it++)
    {
        MolecularModeling::Atom* atom = *it;
        if(atom->GetCoordinates().at(this->grid_->GetAssembly()->GetModelIndex())->GetX() <= this->GetMaxCorner()->GetX() &&
                atom->GetCoordinates().at(this->grid_->GetAssembly()->GetModelIndex())->GetY() <= this->GetMaxCorner()->GetY() &&
                atom->GetCoordinates().at(this->grid_->GetAssembly()->GetModelIndex())->GetZ() <= this->GetMaxCorner()->GetZ() &&
                atom->GetCoordinates().at(this->grid_->GetAssembly()->GetModelIndex())->GetX() > this->GetMinCorner()->GetX() &&
                atom->GetCoordinates().at(this->grid_->GetAssembly()->GetModelIndex())->GetY() > this->GetMinCorner()->GetY() &&
                atom->GetCoordinates().at(this->grid_->GetAssembly()->GetModelIndex())->GetZ() > this->GetMinCorner()->GetZ())
        charge += (atom->MolecularDynamicAtom::GetCharge() != gmml::dNotSet) ? atom->MolecularDynamicAtom::GetCharge() : 0.0;
    }
    cell_charge_ = charge;
}

void Cell::CalculateBoxPotentialEnergy()
{
    double potential_energy = 0.0;
    MolecularModeling::AtomVector all_atoms = this->grid_->GetAssembly()->GetAllAtomsOfAssembly();
    Coordinate* center_of_cell = this->GetCellCenter();
    for(MolecularModeling::AtomVector::iterator it = all_atoms.begin(); it != all_atoms.end(); it++)
    {
        MolecularModeling::Atom* atom = *it;
        double dist = sqrt((center_of_cell->GetX() - atom->GetCoordinates().at(this->grid_->GetAssembly()->GetModelIndex())->GetX()) *
                           (center_of_cell->GetX() - atom->GetCoordinates().at(this->grid_->GetAssembly()->GetModelIndex())->GetX()) +
                           (center_of_cell->GetY() - atom->GetCoordinates().at(this->grid_->GetAssembly()->GetModelIndex())->GetY()) *
                           (center_of_cell->GetY() - atom->GetCoordinates().at(this->grid_->GetAssembly()->GetModelIndex())->GetY()) +
                           (center_of_cell->GetZ() - atom->GetCoordinates().at(this->grid_->GetAssembly()->GetModelIndex())->GetZ()) *
                           (center_of_cell->GetZ() - atom->GetCoordinates().at(this->grid_->GetAssembly()->GetModelIndex())->GetZ()));
        if(dist == 0.0)
            dist = gmml::DIST_EPSILON;
        else
            potential_energy += ((atom->MolecularDynamicAtom::GetCharge() != gmml::dNotSet) ? atom->MolecularDynamicAtom::GetCharge() / dist : 0.0);
    }
    cell_potential_energy_ = potential_energy;
}

//////////////////////////////////////////////////////////
//                     DISPLAY FUNCTIONS                //
//////////////////////////////////////////////////////////
void Cell::Print(std::ostream &out)
{
	out << "";
}
