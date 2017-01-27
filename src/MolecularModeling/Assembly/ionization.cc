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
#include "../../../includes/InputSet/CondensedSequenceSpace/condensedsequenceamberprepresidue.hpp"
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
#include "../../../includes/InputSet/PdbFileSpace/pdbtitlecard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbmodelcard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbmodel.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbmodelresidueset.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbatomcard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbheterogenatomcard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbatom.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbconnectcard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdblinkcard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdblink.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdblinkresidue.hpp"
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

using namespace std;
using namespace MolecularModeling;
using namespace TopologyFileSpace;
using namespace CoordinateFileSpace;
using namespace PrepFileSpace;
using namespace PdbFileSpace;
using namespace PdbqtFileSpace;
using namespace ParameterFileSpace;
using namespace GeometryTopology;
using namespace LibraryFileSpace;
using namespace gmml;
using namespace Glycan;
using namespace CondensedSequenceSpace;

//////////////////////////////////////////////////////////
//                       FUNCTIONS                      //
//////////////////////////////////////////////////////////
void Assembly::AddIon(string ion_name, string lib_file, string parameter_file, int ion_count)
{
    if(ion_count == 0)
    {
        gmml::log(__LINE__, __FILE__,  gmml::INF, "Neutralizing .......");
        cout << "Neutralizing ......." << endl;
        LibraryFile* lib = new LibraryFile(lib_file);
        ParameterFile* param = new ParameterFile(parameter_file, gmml::IONICMOD);
        double charge = this->GetTotalCharge();
        if(fabs(charge) < CHARGE_TOLERANCE)
        {
            gmml::log(__LINE__, __FILE__,  gmml::INF, "The assembly has 0 charge and is neutral.");
            cout << "The assembly has 0 charge and is neutral." << endl;
            return;
        }
        else
        {
            stringstream ss;
            ss << "Total charge of the assembly is " << charge;
            gmml::log(__LINE__, __FILE__,  gmml::INF, ss.str());
            cout << ss.str() << endl;
        }
        double ion_charge = 0;
        string ion_residue_name = "";
        vector<string> ion_list = lib->GetAllResidueNames();
        if(find(ion_list.begin(), ion_list.end(), ion_name) != ion_list.end())
        {
            LibraryFileResidue* lib_ion_residue = lib->GetLibraryResidueByResidueName(ion_name);
            ion_charge = lib_ion_residue->GetLibraryAtomByAtomName(ion_name)->GetCharge();
            ion_residue_name = lib_ion_residue->GetName();

            if(ion_charge == 0)
            {
                gmml::log(__LINE__, __FILE__,  gmml::INF, "The ion has 0 charge");
                cout << "The ion has 0 charge" << endl;
                return;
            }
            else if(ion_charge > 0 && charge > 0)
            {
                gmml::log(__LINE__, __FILE__,  gmml::ERR, "The assembly and the given ion have positive charges, neutralizing process is aborted.");
                cout << "The assembly and the given ion have positive charges, neutralizing process is aborted." << endl;
                return;
            }
            else if(ion_charge < 0 && charge < 0)
            {
                gmml::log(__LINE__, __FILE__,  gmml::ERR, "The assembly and the given ion have positive charges, neutralizing process is aborted.");
                cout << "The assembly and the given ion have negative charges, neutralizing process is aborted." << endl;
                return;
            }
            else
            {
                int number_of_neutralizing_ion = (int)(fabs(charge) + gmml::CHARGE_TOLERANCE) / (int)(fabs(ion_charge) + gmml::CHARGE_TOLERANCE);
                stringstream ss;
                ss << "The assembly will be neutralized by " << number_of_neutralizing_ion << " ion(s)";
                gmml::log(__LINE__, __FILE__,  gmml::INF, ss.str());

                ParameterFile::AtomTypeMap atom_type_map = param->GetAtomTypes();
                double ion_radius = MINIMUM_RADIUS;
                double ion_mass = dNotSet;
                if(atom_type_map.find(ion_name) != atom_type_map.end())
                {
                    ion_radius = atom_type_map[ion_name]->GetRadius();
                    ion_mass = atom_type_map[ion_name]->GetMass();
                }
                Coordinate* minimum_boundary = new Coordinate();
                Coordinate* maximum_boundary = new Coordinate();
                this->GetBoundary(minimum_boundary, maximum_boundary);
                if(minimum_boundary->GetX() == INFINITY || minimum_boundary->GetY() == INFINITY || minimum_boundary->GetZ() == INFINITY ||
                        maximum_boundary->GetX() == -INFINITY || maximum_boundary->GetY() == -INFINITY || maximum_boundary->GetZ() == -INFINITY)
                    return;

                minimum_boundary->operator +(-GRID_OFFSET - 2 * ion_radius - MARGIN);
                maximum_boundary->operator +(GRID_OFFSET + 2 * ion_radius + MARGIN);

                for(int i = 0; i < number_of_neutralizing_ion; i++)
                {
                    Grid* grid = new Grid(this, minimum_boundary, maximum_boundary, ion_radius, ion_charge);
                    grid->CalculateCellsCharge();
                    grid->CalculateCellsPotentialEnergy(ion_radius);
                    CoordinateVector best_positions = grid->GetBestPositions(ion_charge);

                    if(best_positions.size() == 0)
                    {
                        gmml::log(__LINE__, __FILE__,  gmml::ERR, "There is no optimum position to place the ion");
                        cout << "There is no optimum position to place the ion" << endl;
                        return;
                    }
                    else
                    {
                        int index = rand() % best_positions.size();
                        Coordinate* best_position = new Coordinate(best_positions.at(index)->GetX(),
                                                                   best_positions.at(index)->GetY(), best_positions.at(index)->GetZ());
                        Grid::CellVector cells = grid->GetCells();
                        for(Grid::CellVector::iterator it = cells.begin(); it != cells.end(); it++)
                        {
                            if(best_position->GetX() + CRITICAL_RADIOUS * ion_radius + GRID_OFFSET > (*it)->GetCellCenter()->GetX() &&
                                    best_position->GetY() + CRITICAL_RADIOUS * ion_radius + GRID_OFFSET > (*it)->GetCellCenter()->GetY() &&
                                    best_position->GetZ() + CRITICAL_RADIOUS * ion_radius + GRID_OFFSET > (*it)->GetCellCenter()->GetZ() &&
                                    best_position->GetX() - CRITICAL_RADIOUS * ion_radius - GRID_OFFSET < (*it)->GetCellCenter()->GetX() &&
                                    best_position->GetY() - CRITICAL_RADIOUS * ion_radius - GRID_OFFSET < (*it)->GetCellCenter()->GetY() &&
                                    best_position->GetZ() - CRITICAL_RADIOUS * ion_radius - GRID_OFFSET < (*it)->GetCellCenter()->GetZ())
                            {
                                (*it)->SetCellPotentialEnergy(INFINITY);
                            }
                        }

                        Residue* ion = new Residue(this, ion_residue_name);
                        AtomVector atoms = AtomVector();
                        stringstream residue_id;
                        residue_id << ion->GetName() << "_" << BLANK_SPACE << "_" << (i+1) << "_" << BLANK_SPACE << "_" << BLANK_SPACE << "_" << id_;
                        ion->SetId(residue_id.str());

                        CoordinateVector atom_coordinates = CoordinateVector();
                        atom_coordinates.push_back(best_position);
                        Atom* ion_atom = new Atom(ion, lib_ion_residue->GetLibraryAtomByAtomName(ion_name)->GetName(), atom_coordinates);
                        ion_atom->MolecularDynamicAtom::SetAtomType(lib_ion_residue->GetLibraryAtomByAtomName(ion_name)->GetType());
                        ion_atom->MolecularDynamicAtom::SetCharge(lib_ion_residue->GetLibraryAtomByAtomName(ion_name)->GetCharge());
                        ion_atom->MolecularDynamicAtom::SetMass(ion_mass);
                        ion_atom->MolecularDynamicAtom::SetRadius(ion_radius);
                        stringstream atom_id;
                        atom_id << ion_atom->GetName() << "_" << MAX_PDB_ATOM - i << "_" << residue_id.str();
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
            cout << "The ion has not been found in the library file." << endl;
        }
    }
    else if (ion_count > 0)
    {
        gmml::log(__LINE__, __FILE__,  gmml::INF, "Ionizing .......");
        cout << "Ionizing ......." << endl;
        LibraryFile* lib = new LibraryFile(lib_file);
        ParameterFile* param = new ParameterFile(parameter_file, gmml::IONICMOD);
        double charge = this->GetTotalCharge();
        stringstream ss;
        ss << "Total charge of the assembly is " << charge;
        gmml::log(__LINE__, __FILE__,  gmml::INF, ss.str());
        cout << ss.str() << endl;
        double ion_charge = 0;
        string ion_residue_name = "";
        vector<string> ion_list = lib->GetAllResidueNames();
        if(find(ion_list.begin(), ion_list.end(), ion_name) != ion_list.end())
        {
            LibraryFileResidue* lib_ion_residue = lib->GetLibraryResidueByResidueName(ion_name);
            ion_charge = lib_ion_residue->GetLibraryAtomByAtomName(ion_name)->GetCharge();
            ion_residue_name = lib_ion_residue->GetName();

            if(ion_charge == 0)
            {
                gmml::log(__LINE__, __FILE__,  gmml::INF, "The ion has 0 charge");
                cout << "The ion has 0 charge" << endl;
                return;
            }
            else
            {
                stringstream ss;
                ss << "The assembly will be charged by " << ion_count << " ion(s)" ;
                gmml::log(__LINE__, __FILE__,  gmml::INF, ss.str());
                cout << ss.str() << endl;

                ParameterFile::AtomTypeMap atom_type_map = param->GetAtomTypes();
                double ion_radius = MINIMUM_RADIUS;
                double ion_mass = dNotSet;
                if(atom_type_map.find(ion_name) != atom_type_map.end())
                {
                    ion_radius = atom_type_map[ion_name]->GetRadius();
                    ion_mass = atom_type_map[ion_name]->GetMass();
                }
                Coordinate* minimum_boundary = new Coordinate();
                Coordinate* maximum_boundary = new Coordinate();
                this->GetBoundary(minimum_boundary, maximum_boundary);

                if(minimum_boundary->GetX() == INFINITY || minimum_boundary->GetY() == INFINITY || minimum_boundary->GetZ() == INFINITY ||
                        maximum_boundary->GetX() == -INFINITY || maximum_boundary->GetY() == -INFINITY || maximum_boundary->GetZ() == -INFINITY)
                    return;
                minimum_boundary->operator +(-GRID_OFFSET - 2 * ion_radius - MARGIN);
                maximum_boundary->operator +(GRID_OFFSET + 2 * ion_radius + MARGIN);

                for(int i = 0; i < ion_count; i++)
                {
                    Grid* grid = new Grid(this, minimum_boundary, maximum_boundary, ion_radius, ion_charge);
                    grid->CalculateCellsCharge();
                    grid->CalculateCellsPotentialEnergy(ion_radius);
                    CoordinateVector best_positions = grid->GetBestPositions(ion_charge);

                    if(best_positions.size() == 0)
                    {
                        gmml::log(__LINE__, __FILE__,  gmml::ERR, "There is no optimum position to place the ion");
                        cout << "There is no optimum position to place the ion" << endl;
                        return;
                    }
                    else
                    {
                        int index = rand() % best_positions.size();
                        Coordinate* best_position = new Coordinate(best_positions.at(index)->GetX(),
                                                                   best_positions.at(index)->GetY(), best_positions.at(index)->GetZ());
                        Grid::CellVector cells = grid->GetCells();
                        for(Grid::CellVector::iterator it = cells.begin(); it != cells.end(); it++)
                        {
                            if(best_position->GetX() + CRITICAL_RADIOUS * ion_radius + GRID_OFFSET > (*it)->GetCellCenter()->GetX() &&
                                    best_position->GetY() + CRITICAL_RADIOUS * ion_radius + GRID_OFFSET > (*it)->GetCellCenter()->GetY() &&
                                    best_position->GetZ() + CRITICAL_RADIOUS * ion_radius + GRID_OFFSET > (*it)->GetCellCenter()->GetZ() &&
                                    best_position->GetX() - CRITICAL_RADIOUS * ion_radius - GRID_OFFSET < (*it)->GetCellCenter()->GetX() &&
                                    best_position->GetY() - CRITICAL_RADIOUS * ion_radius - GRID_OFFSET < (*it)->GetCellCenter()->GetY() &&
                                    best_position->GetZ() - CRITICAL_RADIOUS * ion_radius - GRID_OFFSET < (*it)->GetCellCenter()->GetZ())
                            {
                                (*it)->SetCellPotentialEnergy(INFINITY);
                            }
                        }

                        Residue* ion = new Residue(this, ion_residue_name);
                        AtomVector atoms = AtomVector();
                        stringstream residue_id;
                        residue_id << ion->GetName() << "_" << BLANK_SPACE << "_" << (i+1) << "_" << BLANK_SPACE << "_" << BLANK_SPACE << "_" << id_;
                        ion->SetId(residue_id.str());

                        CoordinateVector atom_coordinates = CoordinateVector();
                        atom_coordinates.push_back(best_position);
                        Atom* ion_atom = new Atom(ion, lib_ion_residue->GetLibraryAtomByAtomName(ion_name)->GetName(), atom_coordinates);
                        ion_atom->MolecularDynamicAtom::SetAtomType(lib_ion_residue->GetLibraryAtomByAtomName(ion_name)->GetType());
                        ion_atom->MolecularDynamicAtom::SetCharge(lib_ion_residue->GetLibraryAtomByAtomName(ion_name)->GetCharge());
                        ion_atom->MolecularDynamicAtom::SetMass(ion_mass);
                        ion_atom->MolecularDynamicAtom::SetRadius(ion_radius);
                        stringstream atom_id;
                        atom_id << ion_atom->GetName() << "_" << MAX_PDB_ATOM - i << "_" << residue_id.str();
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
            cout << "The ion has not been found in the library file." << endl;
        }
    }
    else
    {
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Please have a non-negative number as the number of ion(s) want to add");
        cout << "Please have a non-negative number as the number of ion(s) want to add" << endl;
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

