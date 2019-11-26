#include <math.h>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <set>
#include <queue>
#include <stack>

#include "../../../includes/MolecularModeling/assembly.hpp"
#include "../../../includes/MolecularModeling/residue.hpp"
#include "../../../includes/MolecularModeling/atom.hpp"
#include "../../../includes/MolecularModeling/atomnode.hpp"
#include "../../../includes/InputSet/CondensedSequenceSpace/condensedsequence.hpp"
#include "../../../includes/InputSet/CondensedSequenceSpace/condensedsequenceresidue.hpp"
#include "../../../includes/InputSet/CondensedSequenceSpace/condensedsequenceglycam06residue.hpp"
#include "../../../includes/InputSet/TopologyFileSpace/topologyfile.hpp"
#include "../../../includes/InputSet/TopologyFileSpace/topologyassembly.hpp"
#include "../../../includes/InputSet/TopologyFileSpace/topologyresidue.hpp"
#include "../../../includes/InputSet/TopologyFileSpace/topologyatom.hpp"
#include "../../../includes/InputSet/TopologyFileSpace/topologybond.hpp"
#include "../../../includes/InputSet/TopologyFileSpace/topologybondtype.hpp"
#include "../../../includes/InputSet/TopologyFileSpace/topologyangle.hpp"
#include "../../../includes/InputSet/TopologyFileSpace/topologyangletype.hpp"
#include "../../../includes/InputSet/TopologyFileSpace/topologydihedral.hpp"
#include "../../../includes/InputSet/TopologyFileSpace/topologydihedraltype.hpp"
#include "../../../includes/InputSet/TopologyFileSpace/topologyatompair.hpp"
#include "../../../includes/InputSet/CoordinateFileSpace/coordinatefile.hpp"
#include "../../../includes/ParameterSet/PrepFileSpace/prepfile.hpp"
#include "../../../includes/ParameterSet/PrepFileSpace/prepfileresidue.hpp"
#include "../../../includes/ParameterSet/PrepFileSpace/prepfileatom.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbfile.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbtitlesection.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbmodelcard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbmodelsection.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbmodelresidueset.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbatomcard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbheterogenatomsection.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbatomsection.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbconnectsection.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdblinkcard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdblinksection.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdblinkcardresidue.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbfileprocessingexception.hpp"
#include "../../../includes/InputSet/PdbqtFileSpace/pdbqtfile.hpp"
#include "../../../includes/InputSet/PdbqtFileSpace/pdbqtatom.hpp"
#include "../../../includes/InputSet/PdbqtFileSpace/pdbqtmodel.hpp"
#include "../../../includes/InputSet/PdbqtFileSpace/pdbqtmodelcard.hpp"
#include "../../../includes/InputSet/PdbqtFileSpace/pdbqtatomcard.hpp"
#include "../../../includes/InputSet/PdbqtFileSpace/pdbqtmodelresidueset.hpp"
#include "../../../includes/ParameterSet/LibraryFileSpace/libraryfile.hpp"
#include "../../../includes/ParameterSet/LibraryFileSpace/libraryfileatom.hpp"
#include "../../../includes/ParameterSet/LibraryFileSpace/libraryfileresidue.hpp"
#include "../../../includes/ParameterSet/ParameterFileSpace/parameterfile.hpp"
#include "../../../includes/ParameterSet/ParameterFileSpace/parameterfilebond.hpp"
#include "../../../includes/ParameterSet/ParameterFileSpace/parameterfileangle.hpp"
#include "../../../includes/ParameterSet/ParameterFileSpace/parameterfiledihedral.hpp"
#include "../../../includes/ParameterSet/ParameterFileSpace/parameterfiledihedralterm.hpp"
#include "../../../includes/ParameterSet/ParameterFileSpace/parameterfileatom.hpp"
#include "../../../includes/ParameterSet/PrepFileSpace/prepfile.hpp"
#include "../../../includes/ParameterSet/PrepFileSpace/prepfileresidue.hpp"
#include "../../../includes/ParameterSet/PrepFileSpace/prepfileatom.hpp"
#include "../../../includes/utils.hpp"
#include "../../../includes/common.hpp"
#include "../../../includes/GeometryTopology/grid.hpp"
#include "../../../includes/GeometryTopology/cell.hpp"

#include <unistd.h>
#include <errno.h>
#include <string.h>

using MolecularModeling::Assembly;

//////////////////////////////////////////////////////////
//                       FUNCTIONS                      //
//////////////////////////////////////////////////////////
void Assembly::AddIon(std::string ion_name, std::string lib_file, std::string parameter_file, int ion_count)
{
    if(ion_count == 0)
    {
        gmml::log(__LINE__, __FILE__,  gmml::INF, "Neutralizing .......");
//        std::cout << "Neutralizing ......." << std::endl;
        LibraryFileSpace::LibraryFile* lib = new LibraryFileSpace::LibraryFile(lib_file);
        ParameterFileSpace::ParameterFile* param = new ParameterFileSpace::ParameterFile(parameter_file, gmml::IONICMOD);
        double charge = this->GetTotalCharge();
        if(fabs(charge) < gmml::CHARGE_TOLERANCE)
        {
            gmml::log(__LINE__, __FILE__,  gmml::INF, "The assembly has 0 charge and is neutral.");
//            std::cout << "The assembly has 0 charge and is neutral." << std::endl;
            return;
        }
        else
        {
            std::stringstream ss;
            ss << "Total charge of the assembly is " << charge;
            gmml::log(__LINE__, __FILE__,  gmml::INF, ss.str());
//            std::cout << ss.str() << std::endl;
        }
        double ion_charge = 0;
        std::string ion_residue_name = "";
        std::vector<std::string> ion_list = lib->GetAllResidueNames();
        if(find(ion_list.begin(), ion_list.end(), ion_name) != ion_list.end())
        {
            LibraryFileSpace::LibraryFileResidue* lib_ion_residue = lib->GetLibraryResidueByResidueName(ion_name);
            ion_charge = lib_ion_residue->GetLibraryAtomByAtomName(ion_name)->GetCharge();
            ion_residue_name = lib_ion_residue->GetName();

            if(ion_charge == 0)
            {
                gmml::log(__LINE__, __FILE__,  gmml::INF, "The ion has 0 charge");
//                std::cout << "The ion has 0 charge" << std::endl;
                return;
            }
            else if(ion_charge > 0 && charge > 0)
            {
                gmml::log(__LINE__, __FILE__,  gmml::ERR, "The assembly and the given ion have positive charges, neutralizing process is aborted.");
//                std::cout << "The assembly and the given ion have positive charges, neutralizing process is aborted." << std::endl;
                return;
            }
            else if(ion_charge < 0 && charge < 0)
            {
                gmml::log(__LINE__, __FILE__,  gmml::ERR, "The assembly and the given ion have positive charges, neutralizing process is aborted.");
//                std::cout << "The assembly and the given ion have negative charges, neutralizing process is aborted." << std::endl;
                return;
            }
            else
            {
                int number_of_neutralizing_ion = (int)(fabs(charge) + gmml::CHARGE_TOLERANCE) / (int)(fabs(ion_charge) + gmml::CHARGE_TOLERANCE);
                std::stringstream ss;
                ss << "The assembly will be neutralized by " << number_of_neutralizing_ion << " ion(s)";
                gmml::log(__LINE__, __FILE__,  gmml::INF, ss.str());

                ParameterFileSpace::ParameterFile::AtomTypeMap atom_type_map = param->GetAtomTypes();
                double ion_radius = gmml::MINIMUM_RADIUS;
                double ion_mass = gmml::dNotSet;
                if(atom_type_map.find(ion_name) != atom_type_map.end())
                {
                    ion_radius = atom_type_map[ion_name]->GetRadius();
                    ion_mass = atom_type_map[ion_name]->GetMass();
                }
                GeometryTopology::Coordinate* minimum_boundary = new GeometryTopology::Coordinate();
                GeometryTopology::Coordinate* maximum_boundary = new GeometryTopology::Coordinate();
                this->GetBoundary(minimum_boundary, maximum_boundary);
                if(minimum_boundary->GetX() == INFINITY || minimum_boundary->GetY() == INFINITY || minimum_boundary->GetZ() == INFINITY ||
                        maximum_boundary->GetX() == -INFINITY || maximum_boundary->GetY() == -INFINITY || maximum_boundary->GetZ() == -INFINITY)
                    return;

                minimum_boundary->operator +(-gmml::GRID_OFFSET - 2 * ion_radius - gmml::MARGIN);
                maximum_boundary->operator +(gmml::GRID_OFFSET + 2 * ion_radius + gmml::MARGIN);

                for(int i = 0; i < number_of_neutralizing_ion; i++)
                {
                    GeometryTopology::Grid* grid = new GeometryTopology::Grid(this, minimum_boundary, maximum_boundary, ion_radius, ion_charge);
                    grid->CalculateCellsCharge();
                    grid->CalculateCellsPotentialEnergy(ion_radius);
                    GeometryTopology::CoordinateVector best_positions = grid->GetBestPositions(ion_charge);

                    if(best_positions.size() == 0)
                    {
                        gmml::log(__LINE__, __FILE__,  gmml::ERR, "There is no optimum position to place the ion");
//                        std::cout << "There is no optimum position to place the ion" << std::endl;
                        return;
                    }
                    else
                    {
                        int index = rand() % best_positions.size();
                        GeometryTopology::Coordinate* best_position = new GeometryTopology::Coordinate(best_positions.at(index)->GetX(),
                                                                   best_positions.at(index)->GetY(), best_positions.at(index)->GetZ());
                        GeometryTopology::Grid::CellVector cells = grid->GetCells();
                        for(GeometryTopology::Grid::CellVector::iterator it = cells.begin(); it != cells.end(); it++)
                        {
                            if(best_position->GetX() + gmml::CRITICAL_RADIOUS * ion_radius + gmml::GRID_OFFSET > (*it)->GetCellCenter()->GetX() &&
                                    best_position->GetY() + gmml::CRITICAL_RADIOUS * ion_radius + gmml::GRID_OFFSET > (*it)->GetCellCenter()->GetY() &&
                                    best_position->GetZ() + gmml::CRITICAL_RADIOUS * ion_radius + gmml::GRID_OFFSET > (*it)->GetCellCenter()->GetZ() &&
                                    best_position->GetX() - gmml::CRITICAL_RADIOUS * ion_radius - gmml::GRID_OFFSET < (*it)->GetCellCenter()->GetX() &&
                                    best_position->GetY() - gmml::CRITICAL_RADIOUS * ion_radius - gmml::GRID_OFFSET < (*it)->GetCellCenter()->GetY() &&
                                    best_position->GetZ() - gmml::CRITICAL_RADIOUS * ion_radius - gmml::GRID_OFFSET < (*it)->GetCellCenter()->GetZ())
                            {
                                (*it)->SetCellPotentialEnergy(INFINITY);
                            }
                        }

                        Residue* ion = new Residue(this, ion_residue_name);
                        AtomVector atoms = AtomVector();
                        std::stringstream residue_id;
                        residue_id << ion->GetName() << "_" << gmml::BLANK_SPACE << "_" << (i+1) << "_" << gmml::BLANK_SPACE << "_" << gmml::BLANK_SPACE << "_" << id_;
                        ion->SetId(residue_id.str());

                        GeometryTopology::CoordinateVector atom_coordinates;
                        atom_coordinates.push_back(best_position);
                        Atom* ion_atom = new Atom(ion, lib_ion_residue->GetLibraryAtomByAtomName(ion_name)->GetName(), atom_coordinates);
                        ion_atom->MolecularDynamicAtom::SetAtomType(lib_ion_residue->GetLibraryAtomByAtomName(ion_name)->GetType());
                        ion_atom->MolecularDynamicAtom::SetCharge(lib_ion_residue->GetLibraryAtomByAtomName(ion_name)->GetCharge());
                        ion_atom->MolecularDynamicAtom::SetMass(ion_mass);
                        ion_atom->MolecularDynamicAtom::SetRadius(ion_radius);
                        std::stringstream atom_id;
                        atom_id << ion_atom->GetName() << "_" << gmml::MAX_PDB_ATOM - i << "_" << residue_id.str();
                        ion_atom->SetId(atom_id.str());

                        atoms.push_back(ion_atom);
                        ion->SetAtoms(atoms);

                        this->AddResidue(ion);
                    }
                }
            }
        }
        else
        {
            gmml::log(__LINE__, __FILE__,  gmml::WAR, "The ion has not been found in the library file.");
//            std::cout << "The ion has not been found in the library file." << std::endl;
        }
    }
    else if (ion_count > 0)
    {
        gmml::log(__LINE__, __FILE__,  gmml::INF, "Ionizing .......");
//        std::cout << "Ionizing ......." << std::endl;
        LibraryFileSpace::LibraryFile* lib = new LibraryFileSpace::LibraryFile(lib_file);
        ParameterFileSpace::ParameterFile* param = new ParameterFileSpace::ParameterFile(parameter_file, gmml::IONICMOD);
        double charge = this->GetTotalCharge();
        std::stringstream ss;
        ss << "Total charge of the assembly is " << charge;
        gmml::log(__LINE__, __FILE__,  gmml::INF, ss.str());
//        std::cout << ss.str() << std::endl;
        double ion_charge = 0;
        std::string ion_residue_name = "";
        std::vector<std::string> ion_list = lib->GetAllResidueNames();
        if(find(ion_list.begin(), ion_list.end(), ion_name) != ion_list.end())
        {
            LibraryFileSpace::LibraryFileResidue* lib_ion_residue = lib->GetLibraryResidueByResidueName(ion_name);
            ion_charge = lib_ion_residue->GetLibraryAtomByAtomName(ion_name)->GetCharge();
            ion_residue_name = lib_ion_residue->GetName();

            if(ion_charge == 0)
            {
                gmml::log(__LINE__, __FILE__,  gmml::INF, "The ion has 0 charge");
//                std::cout << "The ion has 0 charge" << std::endl;
                return;
            }
            else
            {
                std::stringstream ss;
                ss << "The assembly will be charged by " << ion_count << " ion(s)" ;
                gmml::log(__LINE__, __FILE__,  gmml::INF, ss.str());
//                std::cout << ss.str() << std::endl;

                ParameterFileSpace::ParameterFile::AtomTypeMap atom_type_map = param->GetAtomTypes();
                double ion_radius = gmml::MINIMUM_RADIUS;
                double ion_mass = gmml::dNotSet;
                if(atom_type_map.find(ion_name) != atom_type_map.end())
                {
                    ion_radius = atom_type_map[ion_name]->GetRadius();
                    ion_mass = atom_type_map[ion_name]->GetMass();
                }
                GeometryTopology::Coordinate* minimum_boundary = new GeometryTopology::Coordinate();
                GeometryTopology::Coordinate* maximum_boundary = new GeometryTopology::Coordinate();
                this->GetBoundary(minimum_boundary, maximum_boundary);

                if(minimum_boundary->GetX() == INFINITY || minimum_boundary->GetY() == INFINITY || minimum_boundary->GetZ() == INFINITY ||
                        maximum_boundary->GetX() == -INFINITY || maximum_boundary->GetY() == -INFINITY || maximum_boundary->GetZ() == -INFINITY)
                    return;
                minimum_boundary->operator +(-gmml::GRID_OFFSET - 2 * ion_radius - gmml::MARGIN);
                maximum_boundary->operator +(gmml::GRID_OFFSET + 2 * ion_radius + gmml::MARGIN);

                for(int i = 0; i < ion_count; i++)
                {
                    GeometryTopology::Grid* grid = new GeometryTopology::Grid(this, minimum_boundary, maximum_boundary, ion_radius, ion_charge);
                    grid->CalculateCellsCharge();
                    grid->CalculateCellsPotentialEnergy(ion_radius);
                    GeometryTopology::CoordinateVector best_positions = grid->GetBestPositions(ion_charge);

                    if(best_positions.size() == 0)
                    {
                        gmml::log(__LINE__, __FILE__,  gmml::ERR, "There is no optimum position to place the ion");
//                        std::cout << "There is no optimum position to place the ion" << std::endl;
                        return;
                    }
                    else
                    {
                        int index = rand() % best_positions.size();
                        GeometryTopology::Coordinate* best_position = new GeometryTopology::Coordinate(best_positions.at(index)->GetX(),
                                                                   best_positions.at(index)->GetY(), best_positions.at(index)->GetZ());
                        GeometryTopology::Grid::CellVector cells = grid->GetCells();
                        for(GeometryTopology::Grid::CellVector::iterator it = cells.begin(); it != cells.end(); it++)
                        {
                            if(best_position->GetX() + gmml::CRITICAL_RADIOUS * ion_radius + gmml::GRID_OFFSET > (*it)->GetCellCenter()->GetX() &&
                                    best_position->GetY() + gmml::CRITICAL_RADIOUS * ion_radius + gmml::GRID_OFFSET > (*it)->GetCellCenter()->GetY() &&
                                    best_position->GetZ() + gmml::CRITICAL_RADIOUS * ion_radius + gmml::GRID_OFFSET > (*it)->GetCellCenter()->GetZ() &&
                                    best_position->GetX() - gmml::CRITICAL_RADIOUS * ion_radius - gmml::GRID_OFFSET < (*it)->GetCellCenter()->GetX() &&
                                    best_position->GetY() - gmml::CRITICAL_RADIOUS * ion_radius - gmml::GRID_OFFSET < (*it)->GetCellCenter()->GetY() &&
                                    best_position->GetZ() - gmml::CRITICAL_RADIOUS * ion_radius - gmml::GRID_OFFSET < (*it)->GetCellCenter()->GetZ())
                            {
                                (*it)->SetCellPotentialEnergy(INFINITY);
                            }
                        }

                        Residue* ion = new Residue(this, ion_residue_name);
                        AtomVector atoms = AtomVector();
                        std::stringstream residue_id;
                        residue_id << ion->GetName() << "_" << gmml::BLANK_SPACE << "_" << (i+1) << "_" << gmml::BLANK_SPACE << "_" << gmml::BLANK_SPACE << "_" << id_;
                        ion->SetId(residue_id.str());

                        GeometryTopology::CoordinateVector atom_coordinates;
                        atom_coordinates.push_back(best_position);
                        Atom* ion_atom = new Atom(ion, lib_ion_residue->GetLibraryAtomByAtomName(ion_name)->GetName(), atom_coordinates);
                        ion_atom->MolecularDynamicAtom::SetAtomType(lib_ion_residue->GetLibraryAtomByAtomName(ion_name)->GetType());
                        ion_atom->MolecularDynamicAtom::SetCharge(lib_ion_residue->GetLibraryAtomByAtomName(ion_name)->GetCharge());
                        ion_atom->MolecularDynamicAtom::SetMass(ion_mass);
                        ion_atom->MolecularDynamicAtom::SetRadius(ion_radius);
                        std::stringstream atom_id;
                        atom_id << ion_atom->GetName() << "_" << gmml::MAX_PDB_ATOM - i << "_" << residue_id.str();
                        ion_atom->SetId(atom_id.str());

                        atoms.push_back(ion_atom);
                        ion->SetAtoms(atoms);

                        this->AddResidue(ion);
                    }
                }
            }
        }
        else
        {
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "The ion has not been found in the library file.");
//            std::cout << "The ion has not been found in the library file." << std::endl;
        }
    }
    else
    {
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Please have a non-negative number as the number of ion(s) want to add");
//        std::cout << "Please have a non-negative number as the number of ion(s) want to add" << std::endl;
    }
}

void Assembly::SplitIons(Assembly *assembly, ResidueVector ions)
{
    for(AssemblyVector::iterator it = this->assemblies_.begin(); it != this->assemblies_.end(); it++)
        (*it)->SplitIons(assembly, ions);
    for(ResidueVector::iterator it = this->residues_.begin(); it != this->residues_.end(); it++)
    {
        Residue* residue = *it;
        if(residue->GetAtoms().size() == 1)
            ions.push_back(residue);
        else
            assembly->AddResidue(residue);
    }
}
