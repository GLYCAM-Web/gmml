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
void Assembly::AddSolvent(double extension, double closeness, string lib_file)
{
    Coordinate* solvent_component_min_boundary = new Coordinate();
    Coordinate* solvent_component_max_boundary = new Coordinate();
    Assembly* solvent_component = new Assembly();                   //Box of water (TIP3p | TIP5p)
    solvent_component->BuildAssemblyFromLibraryFile(lib_file);      //Reading the box of water from library file
    solvent_component->SetSourceFile(lib_file);
    solvent_component->BuildStructureByLIBFileInformation();        //Building the structure of the box of water

    //Bounding box calculation of the water box (TIP3p | TIP5p)
    solvent_component->GetBoundary(solvent_component_min_boundary, solvent_component_max_boundary);

    //Reading the exact dimension of the water box from library file
    LibraryFile* lib = new LibraryFile(lib_file);
    LibraryFile::ResidueMap lib_residues = lib->GetResidues();
    LibraryFileResidue* lib_residue  = lib_residues.begin()->second;

    double solvent_length = lib_residue->GetBoxLength();
    double solvent_width = lib_residue->GetBoxWidth();
    double solvent_height = lib_residue->GetBoxHeight();

    //Bounding box calculation of the solute
    Coordinate* solute_min_boundary = new Coordinate();
    Coordinate* solute_max_boundary = new Coordinate();
    this->GetBoundary(solute_min_boundary, solute_max_boundary);
    double solute_length = solute_max_boundary->GetX() - solute_min_boundary->GetX();
    double solute_width = solute_max_boundary->GetY() - solute_min_boundary->GetY();
    double solute_height = solute_max_boundary->GetZ() - solute_min_boundary->GetZ();

    //Solvent cube dimension calculation
    double solvent_box_dimension = 0;
    if(solute_length/2 > solvent_box_dimension)
        solvent_box_dimension = solute_length/2;
    if(solute_width/2 > solvent_box_dimension)
        solvent_box_dimension = solute_width/2;
    if(solute_height/2 > solvent_box_dimension)
        solvent_box_dimension = solute_height/2;
    solvent_box_dimension += extension;

    //Center of solute calculation
    Coordinate* center_of_box = new Coordinate(solute_min_boundary->GetX() + solute_length/2, solute_min_boundary->GetY() + solute_width/2,
                                               solute_min_boundary->GetZ() + solute_height/2);
    //Creating the solvent cube around the center of solute
    Coordinate* solvent_box_min_boundary = new Coordinate(center_of_box->GetX(), center_of_box->GetY(), center_of_box->GetZ());
    solvent_box_min_boundary->operator +(-solvent_box_dimension);
    Coordinate* solvent_box_max_boundary = new Coordinate(center_of_box->GetX(), center_of_box->GetY(), center_of_box->GetZ());
    solvent_box_max_boundary->operator +(solvent_box_dimension);
    double solvent_box_length = solvent_box_max_boundary->GetX() - solvent_box_min_boundary->GetX();
    double solvent_box_width = solvent_box_max_boundary->GetY() - solvent_box_min_boundary->GetY();
    double solvent_box_height = solvent_box_max_boundary->GetZ() - solvent_box_min_boundary->GetZ();

    //Number of copies of water box along x, y and z axis inside solvent cube
    int x_copy = solvent_box_length/solvent_length + 1;
    int y_copy = solvent_box_width/solvent_width + 1;
    int z_copy = solvent_box_height/solvent_height + 1;

    //Amount of shifting of the default water box along x, y and z axis to be placed in the min corner of the solvent cube
    double shift_x = solvent_box_min_boundary->GetX() - solvent_component_min_boundary->GetX();// - solvent_length/2;
    double shift_y = solvent_box_min_boundary->GetY() - solvent_component_min_boundary->GetY();// - solvent_width/2;
    double shift_z = solvent_box_min_boundary->GetZ() - solvent_component_min_boundary->GetZ();// - solvent_height/2;

    int sequence_number = this->GetResidues().size() + 1;
    int serial_number = this->GetAllAtomsOfAssembly().size() + 1;

    //Filling the solvent cube with water boxes
    for(int i = 0; i < x_copy; i ++)
    {
        for(int j = 0; j < y_copy; j ++)
        {
            for(int k = 0; k < z_copy; k++)
            {
                //solute_min/max_with_extension = solute_min/max_boundary + closeness
                //check coordinate of each atom to see if it's between solute_min/max_with_extension
                //if so, check distance with all solutes atoms (this)
                //if < closeness remove else add to to_be_added_atoms vector
                Coordinate* solute_min_boundary_with_extension = new Coordinate();
                Coordinate* solute_max_boundary_with_extension = new Coordinate();
                solute_min_boundary_with_extension->SetX(solute_min_boundary->GetX() - closeness);
                solute_min_boundary_with_extension->SetY(solute_min_boundary->GetY() - closeness);
                solute_min_boundary_with_extension->SetZ(solute_min_boundary->GetZ() - closeness);
                solute_max_boundary_with_extension->SetX(solute_max_boundary->GetX() + closeness);
                solute_max_boundary_with_extension->SetY(solute_max_boundary->GetY() + closeness);
                solute_max_boundary_with_extension->SetZ(solute_max_boundary->GetZ() + closeness);
                Assembly* tip_box = new Assembly();
                tip_box->BuildAssemblyFromLibraryFile(lib_file);
                tip_box->SetSourceFile(lib_file);
                tip_box->BuildStructureByLIBFileInformation();
                //                    tip_box->BuildStructureByDistance();
                AtomVector all_atoms_of_tip = tip_box->GetAllAtomsOfAssembly();
                Residue* tip_residue = new Residue();
                vector<string> removed_atom_id_list = vector<string>();
                for(AtomVector::iterator it = all_atoms_of_tip.begin(); it != all_atoms_of_tip.end(); it++)
                {
                    (*it)->GetCoordinates().at(model_index_)->Translate(shift_x + i * solvent_length,
                                                                        shift_y + j * solvent_width, shift_z + k * solvent_height);
                    (*it)->SetDescription("Het;");

                    //Check if the atom of the water box residue is outside the solvent cube and mark it as to be removed
                    if(((*it)->GetCoordinates().at(model_index_)->GetX()) >= solvent_box_max_boundary->GetX())
                        removed_atom_id_list.push_back((*it)->GetId());
                    if(((*it)->GetCoordinates().at(model_index_)->GetY()) >= solvent_box_max_boundary->GetY())
                        removed_atom_id_list.push_back((*it)->GetId());
                    if(((*it)->GetCoordinates().at(model_index_)->GetZ()) >= solvent_box_max_boundary->GetZ())
                        removed_atom_id_list.push_back((*it)->GetId());
                    GeometryTopology::Coordinate* tip_atom_coords = (*it)->GetCoordinates().at(model_index_);

                    //Check if the atom of water box residue is overlaping the solute or is not far enough from the boundary of the solute
                    //and mark it as to be removed
                    if(solute_min_boundary_with_extension->GetX() <= tip_atom_coords->GetX() &&
                            tip_atom_coords->GetX() <= solute_max_boundary_with_extension->GetX() &&
                            solute_min_boundary_with_extension->GetY() <= tip_atom_coords->GetY() &&
                            tip_atom_coords->GetY() <= solute_max_boundary_with_extension->GetY() &&
                            solute_min_boundary_with_extension->GetZ() <= tip_atom_coords->GetZ() &&
                            tip_atom_coords->GetZ() <= solute_max_boundary_with_extension->GetZ() )
                    {
                        AtomVector all_atoms_of_solute = this->GetAllAtomsOfAssembly();
                        bool flag = false;
                        for(AtomVector::iterator it1 = all_atoms_of_solute.begin(); it1 != all_atoms_of_solute.end(); it1++)
                        {
                            Atom* solute_atom = (*it1);
                            if(((*it)->GetCoordinates().at(model_index_)->Distance(*(solute_atom->GetCoordinates().at(model_index_)))) <= closeness)
                            {
                                removed_atom_id_list.push_back((*it)->GetId());
                                flag = true;
                                break;
                            }
                        }
                        if(flag)
                            continue;
                    }
                }
                //Check if one atom of HOH is in to-be-removed list and add the other two atoms belonging to the same HOH
                for(AtomVector::iterator it = all_atoms_of_tip.begin(); it != all_atoms_of_tip.end(); it++)
                {
                    Atom* tip_atom = *it;
                    if(find(removed_atom_id_list.begin(), removed_atom_id_list.end(), tip_atom->GetId()) != removed_atom_id_list.end())
                    {
                        AtomNode* node = tip_atom->GetNode();
                        if(node != NULL)
                        {
                            AtomVector neighbors = node->GetNodeNeighbors();
                            for(AtomVector::iterator it1 = neighbors.begin(); it1 != neighbors.end(); it1++)
                            {
                                Atom* neighbor = *it1;
                                removed_atom_id_list.push_back(neighbor->GetId());

                                AtomNode* neighbor_node = neighbor->GetNode();
                                if(neighbor_node != NULL)
                                {
                                    AtomVector neighbors_of_neighbor = neighbor_node->GetNodeNeighbors();
                                    for(AtomVector::iterator it2 = neighbors_of_neighbor.begin(); it2 != neighbors_of_neighbor.end(); it2++)
                                    {
                                        Atom* neighbor_of_neighbor = *it2;
                                        removed_atom_id_list.push_back(neighbor_of_neighbor->GetId());
                                    }
                                }
                            }
                        }
                    }
                }
                //Add all water molecules of a water box which are not in the to-be-removed list to one residue
                for(AtomVector::iterator it = all_atoms_of_tip.begin(); it != all_atoms_of_tip.end(); it++)
                {
                    Atom* tip_atom = *it;
                    if(find(removed_atom_id_list.begin(), removed_atom_id_list.end(), tip_atom->GetId()) == removed_atom_id_list.end())
                    {
                        tip_residue->AddAtom(tip_atom);
                        string residue_name = tip_atom->GetResidue()->GetName().substr(0,4);
                        tip_residue->SetName("HOH");
                        tip_atom->SetResidue(tip_residue);
                        string id = residue_name + "_" + BLANK_SPACE + "_" + ConvertT<int>(sequence_number) + "_" +
                                BLANK_SPACE + "_" + BLANK_SPACE + "_" + this->GetId();
                        tip_residue->SetId(id);
                        string atom_id = tip_atom->GetName() + "_" + ConvertT<int>(serial_number) + "_" + id;
                        serial_number++;
                        tip_atom->SetId(atom_id);
                    }
                }
                //Add the residue to the assembly
                this->AddResidue(tip_residue);
                sequence_number++;
                //                }
            }
        }
    }
}

void Assembly::SplitSolvent(Assembly* solvent, Assembly* solute)
{
    for(AssemblyVector::iterator it = this->assemblies_.begin(); it != this->assemblies_.end(); it++)
    {
        (*it)->SplitSolvent(solvent, solute);
    }
    for(ResidueVector::iterator it = this->residues_.begin(); it != this->residues_.end(); it++)
    {
        Residue* residue = *it;
        if(residue->GetName().compare("HOH") == 0 || residue->GetName().compare("TP3") == 0 ||
                residue->GetName().compare("TP5") == 0)
            solvent->AddResidue(residue);
        else
            solute->AddResidue(residue);
    }
}

