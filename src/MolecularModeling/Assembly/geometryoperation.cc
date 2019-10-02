#include <math.h>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <set>
#include <queue>
#include <stack>
#include <sstream>

#include "../../../includes/MolecularModeling/assembly.hpp"
#include "../../../includes/MolecularModeling/residue.hpp"
#include "../../../includes/MolecularModeling/atom.hpp"
#include "../../../includes/MolecularModeling/atomnode.hpp"
#include "../../../includes/MolecularModeling/overlaps.hpp"
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
#include "../../../includes/MolecularMetadata/GLYCAM/bondlengthbytypepair.hpp"
#include "../../../includes/MolecularMetadata/GLYCAM/amberatomtypeinfo.hpp"
#include "../../../includes/MolecularMetadata/GLYCAM/glycam06residueinfo.hpp"

#include <unistd.h>
#include <errno.h>
#include <string.h>
#include <cstdlib> //Yao Xiao: call exit() instead of throwing exception, if there's an error.

typedef std::vector<MolecularModeling::Atom*> AtomVector;

using MolecularModeling::Assembly;

//////////////////////////////////////////////////////////
//                       FUNCTIONS                      //
//////////////////////////////////////////////////////////
void Assembly::SetDerivativeAngle(MolecularModeling::Residue *residue, MolecularModeling::Residue *parent_residue, int branch_index)
{
    MolecularModeling::Atom* atom3 = residue->GetHeadAtoms().at(0);
    MolecularModeling::Atom* atom2 = parent_residue->GetTailAtoms().at(branch_index);
    MolecularModeling::Atom* atom1 = NULL;
    if(atom2 != NULL)
    {
        int atom2_index = 1;
        if(atom2->GetName().size() > 1 && isdigit(atom2->GetName().at(1)))
            atom2_index = gmml::ConvertString<int>(gmml::ConvertT<char>(atom2->GetName().at(1)));

        AtomVector atom2_neighbors = atom2->GetNode()->GetNodeNeighbors();
        for(AtomVector::iterator it = atom2_neighbors.begin(); it != atom2_neighbors.end(); it++)
        {
            MolecularModeling::Atom* neighbor = *it;
            if(neighbor->GetId().compare(atom3->GetId()) != 0)
            {
                if(neighbor->GetName().at(0) == 'C' &&
                        (neighbor->GetName().size() > 1 && isdigit(neighbor->GetName().at(1)) &&
                         gmml::ConvertString<int>(gmml::ConvertT<char>(neighbor->GetName().at(1))) == atom2_index))
                {
                    atom1 = neighbor;
                    break;
                }
            }
        }
        if(atom1 != NULL && atom2 != NULL && atom3 != NULL)
            this->SetAngle(atom1, atom2, atom3, 124.0);

    }
}

void Assembly::SetAttachedResidueBond(MolecularModeling::Residue *residue, MolecularModeling::Residue *parent_residue, int branch_index, std::string parameter_file)
{
    ParameterFileSpace::ParameterFile* parameter = new ParameterFileSpace::ParameterFile(parameter_file);
    ParameterFileSpace::ParameterFile::BondMap parameter_bonds = parameter->GetBonds();
    MolecularModeling::Atom* residue_head_atom = residue->GetHeadAtoms().at(0);
    MolecularModeling::Atom* parent_target_atom = parent_residue->GetTailAtoms().at(branch_index);
    AtomVector residue_head_atom_adjacent_atoms = AtomVector();
    AtomVector parent_target_atom_adjacent_atoms = AtomVector();
    double bond_length = gmml::BOND_LENGTH;
    std::vector<std::string> bond = std::vector<std::string>();
    bond.push_back(residue_head_atom->GetAtomType());
    bond.push_back(parent_target_atom->GetAtomType());
    std::vector<std::string> reverse_bond = std::vector<std::string>();
    reverse_bond.push_back(parent_target_atom->GetAtomType());
    reverse_bond.push_back(residue_head_atom->GetAtomType());
    if(parameter_bonds.find(bond) != parameter_bonds.end())
        bond_length = parameter_bonds[bond]->GetLength();
    else if(parameter_bonds.find(reverse_bond) != parameter_bonds.end())
        bond_length = parameter_bonds[reverse_bond]->GetLength();
    residue->GetHeadAtoms().at(0)->GetNode()->AddNodeNeighbor(parent_target_atom);
    parent_residue->GetTailAtoms().at(branch_index)->GetNode()->AddNodeNeighbor(residue_head_atom);
    residue_head_atom_adjacent_atoms = residue_head_atom->GetNode()->GetNodeNeighbors();
    parent_target_atom_adjacent_atoms = parent_target_atom->GetNode()->GetNodeNeighbors();

    GeometryTopology::Coordinate* residue_direction = new GeometryTopology::Coordinate();
    for(AtomVector::iterator it = residue_head_atom_adjacent_atoms.begin(); it != residue_head_atom_adjacent_atoms.end(); it++)
    {
        MolecularModeling::Atom* atom = *it;
        if(atom->GetId().compare(parent_target_atom->GetId()) != 0)
        {
            GeometryTopology::Coordinate* dist = new GeometryTopology::Coordinate(*residue_head_atom->GetCoordinates().at(model_index_));
            dist->operator -(*atom->GetCoordinates().at(model_index_));
            dist->Normalize();
            residue_direction->operator +(*dist);
            residue_direction->Normalize();
        }
    }

    residue_direction->Normalize();
    residue_direction->operator *(bond_length);

    residue_direction->operator +(*residue_head_atom->GetCoordinates().at(model_index_));

    GeometryTopology::Coordinate* oxygen_position = new GeometryTopology::Coordinate(residue_direction->GetX(), residue_direction->GetY(), residue_direction->GetZ());
    GeometryTopology::Coordinate* offset = new GeometryTopology::Coordinate(*parent_target_atom->GetCoordinates().at(model_index_));
    offset->operator -(*oxygen_position);

    AtomVector atomsToTranslate = AtomVector();
    atomsToTranslate.push_back(parent_target_atom);
    residue_head_atom->FindConnectedAtoms(atomsToTranslate);

    for(AtomVector::iterator it = atomsToTranslate.begin() + 1; it != atomsToTranslate.end(); it++)
    {
        (*it)->GetCoordinates().at(model_index_)->Translate(offset->GetX(), offset->GetY(), offset->GetZ());
    }
}

void Assembly::SetResidueResidueBondDistance(MolecularModeling::Atom* tail_atom, MolecularModeling::Atom* head_atom_of_child_residue)
{
    //Attention: SetAtomType()function is overloaded as MolecularModeling::Atom::SetAtomType() and MolecularModeling::MolecularDynamicAtom::SetAtomType(). You don't really know which one to use.
    //Likewise, GetAtomType() is also overloaded.
    //In my situation, I called MolecularModeling::MolecularModelingAtom::SetAtomType(), but later called MolecularModeling::Atom::GetAtomType(). The result is empty.
    //We need to talk about this later
    /*
    std::string tail_atom_type = tail_atom->MolecularDynamicAtom::GetAtomType();
    std::string head_atom_type = head_atom_of_child_residue->MolecularDynamicAtom::GetAtomType();
    std::string tail_atom_hybridization = "";
    std::string head_atom_hybridization = "";
    double head_tail_bond_length = gmml::BOND_LENGTH;
    //look up hybridization state based on atom type
    int GLYCAM06J1ATOMTYPESSIZE = sizeof(gmml::MolecularMetadata::GLYCAM::Glycam06j1AtomTypes)/sizeof(gmml::MolecularMetadata::GLYCAM::Glycam06j1AtomTypes[0]);
    for (int j = 0; j < GLYCAM06J1ATOMTYPESSIZE; j++){
        gmml::MolecularMetadata::GLYCAM::AmberAtomTypeInfo entry = gmml::MolecularMetadata::GLYCAM::Glycam06j1AtomTypes[j];
        if (tail_atom_type == entry.type_){
            tail_atom_hybridization = entry.hybridization_;
        }
        if (head_atom_type == entry.type_){
            head_atom_hybridization = entry.hybridization_;
        }
    }
    
    //look up bond length based on hybridization of head and tail.
    int GLYCAM06J1BONDLENGTHSSIZE = sizeof(gmml::MolecularMetadata::GLYCAM::Glycam06j1BondLengths)/sizeof(gmml::MolecularMetadata::GLYCAM::Glycam06j1BondLengths[0]);
    for (int j = 0; j < GLYCAM06J1BONDLENGTHSSIZE; j++){
        gmml::MolecularMetadata::GLYCAM::BondLengthByTypePair entry = gmml::MolecularMetadata::GLYCAM::Glycam06j1BondLengths[j];
        //Search bidirectionally e.g Cg-Os, Os-Cg
        if ( (tail_atom_type == entry.type1_ && head_atom_type == entry.type2_) || (tail_atom_type == entry.type2_ && head_atom_type == entry.type1_) ){
            head_tail_bond_length = entry.length_;
        }
    }
   */   //Above code replaced by Oliver 2018-10-26 with the code below that uses the new metadata classes:
    std::string tail_atom_type = tail_atom->MolecularDynamicAtom::GetAtomType();
    std::string head_atom_type = head_atom_of_child_residue->MolecularDynamicAtom::GetAtomType();

    gmml::MolecularMetadata::GLYCAM::AmberAtomTypeInfoContainer amberAtomTypeInfoContainer;
    gmml::MolecularMetadata::GLYCAM::AmberAtomTypeInfo entry = amberAtomTypeInfoContainer.GetEntryWithAtomType(head_atom_type);
    std::string head_atom_hybridization = entry.hybridization_;
    //entry = amberAtomTypeInfoContainer.GetEntryWithAtomType(tail_atom_type);
    //std::string tail_atom_hybridization = entry.hybridization_; // OG: This is never used...

    gmml::MolecularMetadata::GLYCAM::BondLengthByTypePairContainer bondLengthByTypePairContainer;
    double head_tail_bond_length = bondLengthByTypePairContainer.GetBondLengthForAtomTypes(tail_atom_type, head_atom_type);
    // End Oliver 2018-10-26 replacement code

    //for sp3 and sp2 hybridzation, sum up all bond vectors between tail atom and neighbors excluding head. Then tail-head bond vector is the opposite of the sum.
    //Normallize length of tail-head bond vector to bond length
    //if (tail_atom_hybridization == "sp3" || tail_atom_hybridization == "sp2"){
    if (head_atom_hybridization == "sp3" || head_atom_hybridization == "sp2"){
	this->Grafting (tail_atom ,head_atom_of_child_residue, head_tail_bond_length);
    }
}

void Assembly::Grafting (MolecularModeling::Atom* tail_atom, MolecularModeling::Atom* head_atom_of_child_residue, double head_tail_bond_length)
{

    GeometryTopology::Coordinate* head_atom_coordinate = head_atom_of_child_residue->GetCoordinates().at(0);
    GeometryTopology::Coordinate* tail_head_atom_bond_vector = new GeometryTopology::Coordinate();
    MolecularModeling::AtomVector head_atom_neighbors = head_atom_of_child_residue->GetNode()->GetNodeNeighbors();
    std::vector<GeometryTopology::Coordinate*> head_atom_bonds = std::vector<GeometryTopology::Coordinate*>();
    for(unsigned int i = 0; i < head_atom_neighbors.size(); i++){
        MolecularModeling::Atom* neighbor = head_atom_neighbors[i];
        if (neighbor != tail_atom){
            GeometryTopology::Coordinate* neighbor_coordinate = neighbor->GetCoordinates().at(0);
            GeometryTopology::Coordinate* head_neighbor_bond_vector = new GeometryTopology::Coordinate();
            head_neighbor_bond_vector->SetX(neighbor_coordinate->GetX() - head_atom_coordinate->GetX());
            head_neighbor_bond_vector->SetY(neighbor_coordinate->GetY() - head_atom_coordinate->GetY());
            head_neighbor_bond_vector->SetZ(neighbor_coordinate->GetZ() - head_atom_coordinate->GetZ());
            head_neighbor_bond_vector->Normalize();
            head_atom_bonds.push_back(head_neighbor_bond_vector);
        }
    }
    double total_x=0.0;
    double total_y=0.0;
    double total_z=0.0;
    for (unsigned int j = 0; j < head_atom_bonds.size(); j++){
        total_x += head_atom_bonds[j]->GetX();
        total_y += head_atom_bonds[j]->GetY();
        total_z += head_atom_bonds[j]->GetZ();
    }
    tail_head_atom_bond_vector->SetX(-1 * total_x);
    tail_head_atom_bond_vector->SetY(-1 * total_y);
    tail_head_atom_bond_vector->SetZ(-1 * total_z);
    double vector_length = tail_head_atom_bond_vector-> length();
    tail_head_atom_bond_vector->SetX(head_tail_bond_length * tail_head_atom_bond_vector->GetX() / vector_length);
    tail_head_atom_bond_vector->SetY(head_tail_bond_length * tail_head_atom_bond_vector->GetY() / vector_length);
    tail_head_atom_bond_vector->SetZ(head_tail_bond_length * tail_head_atom_bond_vector->GetZ() / vector_length);

    //Obtain relative tail position
    GeometryTopology::Coordinate* relative_tail_atom_position = new GeometryTopology::Coordinate();
    relative_tail_atom_position->SetX(head_atom_coordinate->GetX() + tail_head_atom_bond_vector->GetX());
    relative_tail_atom_position->SetY(head_atom_coordinate->GetY() + tail_head_atom_bond_vector->GetY());
    relative_tail_atom_position->SetZ(head_atom_coordinate->GetZ() + tail_head_atom_bond_vector->GetZ());
    //Obtain translation vector by current tail positon - relative tail position
    GeometryTopology::Coordinate* current_tail_atom_position = tail_atom->GetCoordinates().at(0);
    GeometryTopology::Coordinate* translation_vector = new GeometryTopology::Coordinate();
    translation_vector->SetX(current_tail_atom_position->GetX() - relative_tail_atom_position->GetX());
    translation_vector->SetY(current_tail_atom_position->GetY() - relative_tail_atom_position->GetY());
    translation_vector->SetZ(current_tail_atom_position->GetZ() - relative_tail_atom_position->GetZ());
    //Translate coordinates of all atoms in child residue by translation vector.
    MolecularModeling::AtomVector child_residue_atoms = head_atom_of_child_residue->GetResidue()->GetAtoms();
    for (unsigned int i = 0; i < child_residue_atoms.size(); i++){
        GeometryTopology::Coordinate* current_position = child_residue_atoms[i]->GetCoordinates().at(0);
        current_position->Translate(translation_vector->GetX(), translation_vector->GetY(), translation_vector->GetZ());
    }
}

void Assembly::SetAttachedResidueAngle(MolecularModeling::Residue *residue, MolecularModeling::Residue *parent_residue, int branch_index, std::string parameter_file)
{
    ParameterFileSpace::ParameterFile* parameter = new ParameterFileSpace::ParameterFile(parameter_file);
    ParameterFileSpace::ParameterFile::BondMap parameter_bonds = parameter->GetBonds();
    MolecularModeling::Atom* residue_head_atom = residue->GetHeadAtoms().at(0);
    MolecularModeling::Atom* parent_target_atom = parent_residue->GetTailAtoms().at(branch_index);
    AtomVector residue_head_atom_adjacent_atoms = AtomVector();
    AtomVector parent_target_atom_adjacent_atoms = AtomVector();

    double bond_length = gmml::BOND_LENGTH;
    std::vector<std::string> bond = std::vector<std::string>();
    bond.push_back(residue_head_atom->GetAtomType());
    bond.push_back(parent_target_atom->GetAtomType());
    std::vector<std::string> reverse_bond = std::vector<std::string>();
    reverse_bond.push_back(parent_target_atom->GetAtomType());
    reverse_bond.push_back(residue_head_atom->GetAtomType());
    if(parameter_bonds.find(bond) != parameter_bonds.end())
        bond_length = parameter_bonds[bond]->GetLength();
    else if(parameter_bonds.find(reverse_bond) != parameter_bonds.end())
        bond_length = parameter_bonds[reverse_bond]->GetLength();
    residue_head_atom_adjacent_atoms = residue_head_atom->GetNode()->GetNodeNeighbors();
    parent_target_atom_adjacent_atoms = parent_target_atom->GetNode()->GetNodeNeighbors();

    GeometryTopology::Coordinate* carbon_direction = new GeometryTopology::Coordinate();
    for(AtomVector::iterator it = parent_target_atom_adjacent_atoms.begin(); it != parent_target_atom_adjacent_atoms.end(); it++)
    {
        MolecularModeling::Atom* atom = *it;
        if(atom->GetId().compare(residue_head_atom->GetId()) != 0)
        {
            GeometryTopology::Coordinate* dist = new GeometryTopology::Coordinate(*parent_target_atom->GetCoordinates().at(model_index_));
            dist->operator -(*atom->GetCoordinates().at(model_index_));
            dist->Normalize();
            carbon_direction->operator +(*dist);
            carbon_direction->Normalize();
        }
    }

    carbon_direction->Normalize();
    carbon_direction->operator *(bond_length);
    carbon_direction->operator +(*parent_target_atom->GetCoordinates().at(model_index_));

    GeometryTopology::Coordinate* carbon_position = new GeometryTopology::Coordinate(*carbon_direction);

    GeometryTopology::Coordinate* carbon_target = new GeometryTopology::Coordinate(*carbon_position);
    carbon_target->operator -(*parent_target_atom->GetCoordinates().at(model_index_));

    GeometryTopology::Coordinate* head_target = new GeometryTopology::Coordinate(*residue_head_atom->GetCoordinates().at(model_index_));
    head_target->operator -(*parent_target_atom->GetCoordinates().at(model_index_));

    double angle = acos((carbon_target->DotProduct(*head_target)) / (carbon_target->length() * head_target->length() + gmml::DIST_EPSILON));
    double rotation_angle = gmml::ConvertDegree2Radian(gmml::PI_DEGREE - gmml::ROTATION_ANGLE) - angle;

    GeometryTopology::Coordinate* direction = new GeometryTopology::Coordinate(*carbon_target);
    direction->CrossProduct(*head_target);
    direction->Normalize();
    double** rotation_matrix = gmml::GenerateRotationMatrix(direction, parent_target_atom->GetCoordinates().at(model_index_), rotation_angle);

    AtomVector atomsToRotate = AtomVector();
    atomsToRotate.push_back(parent_target_atom);
    residue_head_atom->FindConnectedAtoms(atomsToRotate);

    for(AtomVector::iterator it = atomsToRotate.begin() + 1; it != atomsToRotate.end(); it++)
    {
        GeometryTopology::Coordinate* atom_coordinate = (*it)->GetCoordinates().at(model_index_);
        GeometryTopology::Coordinate* result = new GeometryTopology::Coordinate();
        result->SetX(rotation_matrix[0][0] * atom_coordinate->GetX() + rotation_matrix[0][1] * atom_coordinate->GetY() +
                rotation_matrix[0][2] * atom_coordinate->GetZ() + rotation_matrix[0][3]);
        result->SetY(rotation_matrix[1][0] * atom_coordinate->GetX() + rotation_matrix[1][1] * atom_coordinate->GetY() +
                rotation_matrix[1][2] * atom_coordinate->GetZ() + rotation_matrix[1][3]);
        result->SetZ(rotation_matrix[2][0] * atom_coordinate->GetX() + rotation_matrix[2][1] * atom_coordinate->GetY() +
                rotation_matrix[2][2] * atom_coordinate->GetZ() + rotation_matrix[2][3]);

        (*it)->GetCoordinates().at(model_index_)->SetX(result->GetX());
        (*it)->GetCoordinates().at(model_index_)->SetY(result->GetY());
        (*it)->GetCoordinates().at(model_index_)->SetZ(result->GetZ());
    }
}

void Assembly::SetAttachedResidueTorsion(MolecularModeling::Residue *residue, MolecularModeling::Residue *parent_residue, int branch_index)
{
    MolecularModeling::Atom* residue_head_atom = residue->GetHeadAtoms().at(0);
    MolecularModeling::Atom* parent_target_atom = parent_residue->GetTailAtoms().at(branch_index);

    MolecularModeling::Atom* carbon = residue_head_atom; ///The carbon atom of the new residue that is attached to the parent residue
    if(carbon != NULL)
    {
        MolecularModeling::Atom* oxygen = parent_target_atom; ///The oxygen atom of the parent residue that is attached to the new residue
        if(oxygen != NULL)
        {
            ///i: parent residue oxygen atom index from which the new residue is attached to the parent residue
            ///j: attached residue carbon atom index from which the residue is attached to the parent residue
            ///if i == 5 || i == 6
            ///Set C(i-1)-C(i)-O(i)-C(j) dihedral to 180.0
            ///else
            ///Set H(i)-C(i)-O(i)-C(j) dihedral to 0.0

            ///if parent residue is ROH
            ///Set H(i)-O(i)-C(j)-C(j+1) dihedral to 180.0
            ///else
            ///Set C(i)-O(i)-C(j)-C(j+1) dihedral to 180.0

            ///if i == 6
            ///Set O(i-1)-C(i-1)-C(i)-O(i) dihedral to 60.0

            int oxygen_index = 1;
            if(oxygen->GetName().size() > 1 && isdigit(oxygen->GetName().at(1)))
                oxygen_index = gmml::ConvertString<int>(gmml::ConvertT<char>(oxygen->GetName().at(1)));

            if(oxygen_index == 5 || oxygen_index == 6)
            {
                MolecularModeling::Atom* atom1 = NULL;
                MolecularModeling::Atom* atom2 = NULL;
                MolecularModeling::Atom* atom3 = oxygen;
                MolecularModeling::Atom* atom4 = carbon;

                AtomVector oxygen_neighbors = oxygen->GetNode()->GetNodeNeighbors();
                for(AtomVector::iterator it = oxygen_neighbors.begin(); it != oxygen_neighbors.end(); it++)
                {
                    MolecularModeling::Atom* neighbor = *it;
                    if(neighbor->GetId().compare(carbon->GetId()) != 0)
                    {
                        if(neighbor->GetName().at(0) == 'C')
                        {
                            atom2 = neighbor;
                            break;
                        }
                    }
                }
                if(atom2 != NULL)
                {
                    AtomVector atom2_neighbors = atom2->GetNode()->GetNodeNeighbors();
                    for(AtomVector::iterator it = atom2_neighbors.begin(); it != atom2_neighbors.end(); it++)
                    {
                        MolecularModeling::Atom* neighbor = *it;
                        if(neighbor->GetId().compare(oxygen->GetId()) != 0)
                        {
                            if(neighbor->GetName().at(0) == 'C')
                            {
                                atom1 = neighbor;
                                break;
                            }
                        }
                    }
                }
                if(atom1 != NULL && atom2 != NULL && atom3 != NULL && atom4 != NULL)
                    this->SetDihedral(atom1, atom2, atom3, atom4, 180.0);
            }
            else
            {
                MolecularModeling::Atom* atom1 = NULL;
                MolecularModeling::Atom* atom2 = NULL;
                MolecularModeling::Atom* atom3 = oxygen;
                MolecularModeling::Atom* atom4 = carbon;

                AtomVector oxygen_neighbors = oxygen->GetNode()->GetNodeNeighbors();
                for(AtomVector::iterator it = oxygen_neighbors.begin(); it != oxygen_neighbors.end(); it++)
                {
                    MolecularModeling::Atom* neighbor = *it;
                    if(neighbor->GetId().compare(carbon->GetId()) != 0)
                    {
                        if(neighbor->GetName().at(0) == 'C')
                        {
                            atom2 = neighbor;
                            break;
                        }
                    }
                }
                if(atom2 != NULL)
                {
                    AtomVector atom2_neighbors = atom2->GetNode()->GetNodeNeighbors();
                    for(AtomVector::iterator it = atom2_neighbors.begin(); it != atom2_neighbors.end(); it++)
                    {
                        MolecularModeling::Atom* neighbor = *it;
                        if(neighbor->GetId().compare(oxygen->GetId()) != 0)
                        {
                            if(neighbor->GetName().at(0) == 'H')
                            {
                                atom1 = neighbor;
                                break;
                            }
                        }
                    }
                }

                if(atom1 != NULL && atom2 != NULL && atom3 != NULL && atom4 != NULL)
                    this->SetDihedral(atom1, atom2, atom3, atom4, 0.0);
            }

            if(oxygen->GetResidue()->GetName().compare("ROH") == 0)
            {
                MolecularModeling::Atom* atom1 = NULL;
                MolecularModeling::Atom* atom2 = oxygen;
                MolecularModeling::Atom* atom3 = carbon;
                MolecularModeling::Atom* atom4 = NULL;

                AtomVector oxygen_neighbors = oxygen->GetNode()->GetNodeNeighbors();
                for(AtomVector::iterator it = oxygen_neighbors.begin(); it != oxygen_neighbors.end(); it++)
                {
                    MolecularModeling::Atom* neighbor = *it;
                    if(neighbor->GetId().compare(carbon->GetId()) != 0)
                    {
                        if(neighbor->GetName().at(0) == 'H')
                        {
                            atom1 = neighbor;
                            break;
                        }
                    }
                }
                AtomVector carbon_neighbors = carbon->GetNode()->GetNodeNeighbors();
                for(AtomVector::iterator it = carbon_neighbors.begin(); it != carbon_neighbors.end(); it++)
                {
                    MolecularModeling::Atom* neighbor = *it;
                    if(neighbor->GetId().compare(oxygen->GetId()) != 0)
                    {
                        if(neighbor->GetName().at(0) == 'C')
                        {
                            atom4 = neighbor;
                            break;
                        }
                    }
                }

                if(atom1 != NULL && atom2 != NULL && atom3 != NULL && atom4 != NULL)
                    this->SetDihedral(atom1, atom2, atom3, atom4, 180.0);
            }
            else
            {
                MolecularModeling::Atom* atom1 = NULL;
                MolecularModeling::Atom* atom2 = oxygen;
                MolecularModeling::Atom* atom3 = carbon;
                MolecularModeling::Atom* atom4 = NULL;

                AtomVector oxygen_neighbors = oxygen->GetNode()->GetNodeNeighbors();
                for(AtomVector::iterator it = oxygen_neighbors.begin(); it != oxygen_neighbors.end(); it++)
                {
                    MolecularModeling::Atom* neighbor = *it;
                    if(neighbor->GetId().compare(carbon->GetId()) != 0)
                    {
                        if(neighbor->GetName().at(0) == 'C')
                        {
                            atom1 = neighbor;
                            break;
                        }
                    }
                }
                AtomVector carbon_neighbors = carbon->GetNode()->GetNodeNeighbors();
                for(AtomVector::iterator it = carbon_neighbors.begin(); it != carbon_neighbors.end(); it++)
                {
                    MolecularModeling::Atom* neighbor = *it;
                    if(neighbor->GetId().compare(oxygen->GetId()) != 0)
                    {
                        if(neighbor->GetName().at(0) == 'C')
                        {
                            atom4 = neighbor;
                            break;
                        }
                    }
                }

                if(atom1 != NULL && atom2 != NULL && atom3 != NULL && atom4 != NULL)
                    this->SetDihedral(atom1, atom2, atom3, atom4, 180.0);
            }

            if(oxygen_index == 6)
            {
                MolecularModeling::Atom* atom1 = NULL;
                MolecularModeling::Atom* atom2 = NULL;
                MolecularModeling::Atom* atom3 = NULL;
                MolecularModeling::Atom* atom4 = oxygen;

                AtomVector oxygen_neighbors = oxygen->GetNode()->GetNodeNeighbors();
                for(AtomVector::iterator it = oxygen_neighbors.begin(); it != oxygen_neighbors.end(); it++)
                {
                    MolecularModeling::Atom* neighbor = *it;
                    if(neighbor->GetId().compare(carbon->GetId()) != 0)
                    {
                        if(neighbor->GetName().at(0) == 'C')
                        {
                            atom3 = neighbor;
                            break;
                        }
                    }
                }
                if(atom3 != NULL)
                {
                    AtomVector atom3_neighbors = atom3->GetNode()->GetNodeNeighbors();
                    for(AtomVector::iterator it = atom3_neighbors.begin(); it != atom3_neighbors.end(); it++)
                    {
                        MolecularModeling::Atom* neighbor = *it;
                        if(neighbor->GetId().compare(oxygen->GetId()) != 0)
                        {
                            if(neighbor->GetName().at(0) == 'C')
                            {
                                atom2 = neighbor;
                                break;
                            }
                        }
                    }

                    if(atom2 != NULL)
                    {
                        AtomVector atom2_neighbors = atom2->GetNode()->GetNodeNeighbors();
                        for(AtomVector::iterator it = atom2_neighbors.begin(); it != atom2_neighbors.end(); it++)
                        {
                            MolecularModeling::Atom* neighbor = *it;
                            if(neighbor->GetId().compare(atom3->GetId()) != 0)
                            {
                                if(neighbor->GetName().at(0) == 'O')
                                {
                                    atom1 = neighbor;
                                    break;
                                }
                            }
                        }
                    }
                }

                if(atom1 != NULL && atom2 != NULL && atom3 != NULL && atom4 != NULL)
                    this->SetDihedral(atom1, atom2, atom3, atom4, 60.0);
            }
        }
    }
}

void Assembly::SetPhiTorsion(MolecularModeling::Residue *residue, MolecularModeling::Residue *parent_residue, int branch_index, double torsion)
{
    MolecularModeling::Atom* residue_head_atom = residue->GetHeadAtoms().at(0);
    MolecularModeling::Atom* parent_target_atom = parent_residue->GetTailAtoms().at(branch_index);

    MolecularModeling::Atom* carbon = residue_head_atom; ///The carbon atom of the new residue that is attached to the parent residue
    if(carbon != NULL)
    {
        MolecularModeling::Atom* oxygen = parent_target_atom; ///The oxygen atom of the parent residue that is attached to the new residue
        if(oxygen != NULL)
        {
            int oxygen_index = 1;
            if(oxygen->GetName().size() > 1 && isdigit(oxygen->GetName().at(1)))
                oxygen_index = gmml::ConvertString<int>(gmml::ConvertT<char>(oxygen->GetName().at(1)));

            int carbon_index = 1;
            if(carbon->GetName().size() > 1 && isdigit(carbon->GetName().at(1)))
                carbon_index = gmml::ConvertString<int>(gmml::ConvertT<char>(carbon->GetName().at(1)));

            MolecularModeling::Atom* atom1 = NULL;
            MolecularModeling::Atom* atom2 = carbon;
            MolecularModeling::Atom* atom3 = oxygen;
            MolecularModeling::Atom* atom4 = NULL;

            AtomVector oxygen_neighbors = oxygen->GetNode()->GetNodeNeighbors();
            for(AtomVector::iterator it = oxygen_neighbors.begin(); it != oxygen_neighbors.end(); it++)
            {
                MolecularModeling::Atom* neighbor = *it;
                if(neighbor->GetId().compare(carbon->GetId()) != 0)
                {
                    if(neighbor->GetName().at(0) == 'C' &&
                            (neighbor->GetName().size() > 1 && isdigit(neighbor->GetName().at(1)) &&
                             gmml::ConvertString<int>(gmml::ConvertT<char>(neighbor->GetName().at(1))) == oxygen_index))
                    {
                        atom4 = neighbor;
                        break;
                    }
                }
            }
            AtomVector atom2_neighbors = atom2->GetNode()->GetNodeNeighbors();
            for(AtomVector::iterator it = atom2_neighbors.begin(); it != atom2_neighbors.end(); it++)
            {
                MolecularModeling::Atom* neighbor = *it;
                if(neighbor->GetId().compare(oxygen->GetId()) != 0)
                {
                    if(neighbor->GetName().at(0) == 'C' &&
                            (neighbor->GetName().size() > 1 && isdigit(neighbor->GetName().at(1)) &&
                             gmml::ConvertString<int>(gmml::ConvertT<char>(neighbor->GetName().at(1))) == carbon_index - 1))
                    {
                        atom1 = neighbor;
                        break;
                    }
                }
            }

            if(atom1 != NULL && atom2 != NULL && atom3 != NULL && atom4 != NULL)
                this->SetDihedral(atom1, atom2, atom3, atom4, torsion);

        }
    }
}

void Assembly::SetPsiTorsion(MolecularModeling::Residue *residue, MolecularModeling::Residue *parent_residue, int branch_index, double torsion, bool crystallographic_definition)
{
    MolecularModeling::Atom* residue_head_atom = residue->GetHeadAtoms().at(0);
    MolecularModeling::Atom* parent_target_atom = parent_residue->GetTailAtoms().at(branch_index);

    MolecularModeling::Atom* carbon = residue_head_atom; ///The carbon atom of the new residue that is attached to the parent residue
    if(carbon != NULL)
    {
        MolecularModeling::Atom* oxygen = parent_target_atom; ///The oxygen atom of the parent residue that is attached to the new residue
        if(oxygen != NULL)
        {
            int oxygen_index = 1;
            if(oxygen->GetName().size() > 1 && isdigit(oxygen->GetName().at(1)))
                oxygen_index = gmml::ConvertString<int>(gmml::ConvertT<char>(oxygen->GetName().at(1)));

            MolecularModeling::Atom* atom1 = carbon;
            MolecularModeling::Atom* atom2 = oxygen;
            MolecularModeling::Atom* atom3 = NULL;
            MolecularModeling::Atom* atom4 = NULL;

            AtomVector oxygen_neighbors = oxygen->GetNode()->GetNodeNeighbors();
            for(AtomVector::iterator it = oxygen_neighbors.begin(); it != oxygen_neighbors.end(); it++)
            {
                MolecularModeling::Atom* neighbor = *it;
                if(neighbor->GetId().compare(carbon->GetId()) != 0)
                {
                    if(neighbor->GetName().at(0) == 'C' &&
                            (neighbor->GetName().size() > 1 && isdigit(neighbor->GetName().at(1)) &&
                             gmml::ConvertString<int>(gmml::ConvertT<char>(neighbor->GetName().at(1))) == oxygen_index))
                    {
                        atom3 = neighbor;
                        break;
                    }
                }
            }
            if(atom3 != NULL)
            {
                AtomVector atom3_neighbors = atom3->GetNode()->GetNodeNeighbors();
                for(AtomVector::iterator it = atom3_neighbors.begin(); it != atom3_neighbors.end(); it++)
                {
                    MolecularModeling::Atom* neighbor = *it;
                    if(neighbor->GetId().compare(oxygen->GetId()) != 0)
                    {
                        if(crystallographic_definition)
                        {
                            if(neighbor->GetName().at(0) == 'C' &&
                                    (neighbor->GetName().size() > 1 && isdigit(neighbor->GetName().at(1)) &&
                                     gmml::ConvertString<int>(gmml::ConvertT<char>(neighbor->GetName().at(1))) == oxygen_index + 1))
                            {
                                atom4 = neighbor;
                                break;
                            }
                        }
                        else
                        {
                            if(neighbor->GetName().at(0) == 'H' &&
                                    (neighbor->GetName().size() > 1 && isdigit(neighbor->GetName().at(1)) &&
                                     gmml::ConvertString<int>(gmml::ConvertT<char>(neighbor->GetName().at(1))) == oxygen_index))
                            {
                                atom4 = neighbor;
                                break;
                            }
                        }
                    }
                }
            }


            if(atom1 != NULL && atom2 != NULL && atom3 != NULL && atom4 != NULL)
                this->SetDihedral(atom1, atom2, atom3, atom4, torsion);

        }
    }
}

void Assembly::SetOmegaTorsion(MolecularModeling::Residue *residue, MolecularModeling::Residue *parent_residue, int branch_index, double torsion, int type)
{
    MolecularModeling::Atom* residue_head_atom = residue->GetHeadAtoms().at(0);
    MolecularModeling::Atom* parent_target_atom = parent_residue->GetTailAtoms().at(branch_index);

    MolecularModeling::Atom* carbon = residue_head_atom; ///The carbon atom of the new residue that is attached to the parent residue
    if(carbon != NULL)
    {
        MolecularModeling::Atom* oxygen = parent_target_atom; ///The oxygen atom of the parent residue that is attached to the new residue
        if(oxygen != NULL)
        {
            int oxygen_index = 1;
            if(oxygen->GetName().size() > 1 && isdigit(oxygen->GetName().at(1)))
                oxygen_index = gmml::ConvertString<int>(gmml::ConvertT<char>(oxygen->GetName().at(1)));
            if(oxygen_index == 6 && type == 6)
            {
                MolecularModeling::Atom* atom1 = NULL;
                MolecularModeling::Atom* atom2 = NULL;
                MolecularModeling::Atom* atom3 = NULL;
                MolecularModeling::Atom* atom4 = oxygen;

                AtomVector oxygen_neighbors = oxygen->GetNode()->GetNodeNeighbors();
                for(AtomVector::iterator it = oxygen_neighbors.begin(); it != oxygen_neighbors.end(); it++)
                {
                    MolecularModeling::Atom* neighbor = *it;
                    if(neighbor->GetId().compare(carbon->GetId()) != 0)
                    {
                        if(neighbor->GetName().at(0) == 'C' &&
                                (neighbor->GetName().size() > 1 && isdigit(neighbor->GetName().at(1)) &&
                                 gmml::ConvertString<int>(gmml::ConvertT<char>(neighbor->GetName().at(1))) == oxygen_index))
                        {
                            atom3 = neighbor;
                            break;
                        }
                    }
                }
                if(atom3 != NULL)
                {
                    AtomVector atom3_neighbors = atom3->GetNode()->GetNodeNeighbors();
                    for(AtomVector::iterator it = atom3_neighbors.begin(); it != atom3_neighbors.end(); it++)
                    {
                        MolecularModeling::Atom* neighbor = *it;
                        if(neighbor->GetId().compare(oxygen->GetId()) != 0)
                        {
                            if(neighbor->GetName().at(0) == 'C' &&
                                    (neighbor->GetName().size() > 1 && isdigit(neighbor->GetName().at(1)) &&
                                     gmml::ConvertString<int>(gmml::ConvertT<char>(neighbor->GetName().at(1))) == oxygen_index - 1))
                            {
                                atom2 = neighbor;
                                break;
                            }
                        }
                    }
                }
                if(atom2 != NULL)
                {
                    AtomVector atom2_neighbors = atom2->GetNode()->GetNodeNeighbors();
                    for(AtomVector::iterator it = atom2_neighbors.begin(); it != atom2_neighbors.end(); it++)
                    {
                        MolecularModeling::Atom* neighbor = *it;
                        if(neighbor->GetId().compare(atom3->GetId()) != 0)
                        {
                            if(neighbor->GetName().at(0) == 'O' &&
                                    (neighbor->GetName().size() > 1 && isdigit(neighbor->GetName().at(1)) &&
                                     gmml::ConvertString<int>(gmml::ConvertT<char>(neighbor->GetName().at(1))) == oxygen_index - 1))
                            {
                                atom1 = neighbor;
                                break;
                            }
                        }
                    }
                }
                if(atom1 != NULL && atom2 != NULL && atom3 != NULL && atom4 != NULL)
                {
                    this->SetDihedral(atom1, atom2, atom3, atom4, torsion);
                }
            }
            if(oxygen_index == 8 && type == 7)
            {
                MolecularModeling::Atom* atom1 = NULL;
                MolecularModeling::Atom* atom2 = NULL;
                MolecularModeling::Atom* atom3 = NULL;
                MolecularModeling::Atom* atom4 = NULL;
                MolecularModeling::Atom* atom5 = NULL;

                AtomVector oxygen_neighbors = oxygen->GetNode()->GetNodeNeighbors();
                for(AtomVector::iterator it = oxygen_neighbors.begin(); it != oxygen_neighbors.end(); it++)
                {
                    MolecularModeling::Atom* neighbor = *it;
                    if(neighbor->GetId().compare(carbon->GetId()) != 0)
                    {
                        if(neighbor->GetName().at(0) == 'C' &&
                                (neighbor->GetName().size() > 1 && isdigit(neighbor->GetName().at(1)) &&
                                 gmml::ConvertString<int>(gmml::ConvertT<char>(neighbor->GetName().at(1))) == oxygen_index))
                        {
                            atom5 = neighbor;
                            break;
                        }
                    }
                }
                if(atom5 != NULL)
                {
                    AtomVector atom5_neighbors = atom5->GetNode()->GetNodeNeighbors();
                    for(AtomVector::iterator it = atom5_neighbors.begin(); it != atom5_neighbors.end(); it++)
                    {
                        MolecularModeling::Atom* neighbor = *it;
                        if(neighbor->GetId().compare(oxygen->GetId()) != 0)
                        {
                            if(neighbor->GetName().at(0) == 'C' &&
                                    (neighbor->GetName().size() > 1 && isdigit(neighbor->GetName().at(1)) &&
                                     gmml::ConvertString<int>(gmml::ConvertT<char>(neighbor->GetName().at(1))) == oxygen_index - 1))
                            {
                                atom2 = neighbor;
                                break;
                            }
                        }
                    }
                }
                if(atom2 != NULL)
                {
                    AtomVector atom2_neighbors = atom2->GetNode()->GetNodeNeighbors();
                    for(AtomVector::iterator it = atom2_neighbors.begin(); it != atom2_neighbors.end(); it++)
                    {
                        MolecularModeling::Atom* neighbor = *it;
                        if(neighbor->GetId().compare(atom5->GetId()) != 0)
                        {
                            if(neighbor->GetName().at(0) == 'H' &&
                                    (neighbor->GetName().size() > 1 && isdigit(neighbor->GetName().at(1)) &&
                                     gmml::ConvertString<int>(gmml::ConvertT<char>(neighbor->GetName().at(1))) == oxygen_index - 1))
                            {
                                atom1 = neighbor;
                                break;
                            }
                        }
                    }
                    for(AtomVector::iterator it = atom2_neighbors.begin(); it != atom2_neighbors.end(); it++)
                    {
                        MolecularModeling::Atom* neighbor = *it;
                        if(neighbor->GetId().compare(atom5->GetId()) != 0)
                        {
                            if(neighbor->GetName().at(0) == 'C' &&
                                    (neighbor->GetName().size() > 1 && isdigit(neighbor->GetName().at(1)) &&
                                     gmml::ConvertString<int>(gmml::ConvertT<char>(neighbor->GetName().at(1))) == oxygen_index - 2))
                            {
                                atom3 = neighbor;
                                break;
                            }
                        }
                    }
                }
                if(atom3 != NULL)
                {
                    AtomVector atom3_neighbors = atom3->GetNode()->GetNodeNeighbors();
                    for(AtomVector::iterator it = atom3_neighbors.begin(); it != atom3_neighbors.end(); it++)
                    {
                        MolecularModeling::Atom* neighbor = *it;
                        if(neighbor->GetId().compare(atom2->GetId()) != 0)
                        {
                            if(neighbor->GetName().at(0) == 'H' &&
                                    (neighbor->GetName().size() > 1 && isdigit(neighbor->GetName().at(1)) &&
                                     gmml::ConvertString<int>(gmml::ConvertT<char>(neighbor->GetName().at(1))) == oxygen_index - 2))
                            {
                                atom4 = neighbor;
                                break;
                            }
                        }
                    }
                }
                if(atom1 != NULL && atom2 != NULL && atom3 != NULL && atom4 != NULL)
                {
                    this->SetDihedral(atom1, atom2, atom3, atom4, torsion);
                }
            }
            if(oxygen_index == 8 && type == 8)
            {
                MolecularModeling::Atom* atom1 = NULL;
                MolecularModeling::Atom* atom2 = NULL;
                MolecularModeling::Atom* atom3 = NULL;
                MolecularModeling::Atom* atom4 = NULL;

                AtomVector oxygen_neighbors = oxygen->GetNode()->GetNodeNeighbors();
                for(AtomVector::iterator it = oxygen_neighbors.begin(); it != oxygen_neighbors.end(); it++)
                {
                    MolecularModeling::Atom* neighbor = *it;
                    if(neighbor->GetId().compare(carbon->GetId()) != 0)
                    {
                        if(neighbor->GetName().at(0) == 'C' &&
                                (neighbor->GetName().size() > 1 && isdigit(neighbor->GetName().at(1)) &&
                                 gmml::ConvertString<int>(gmml::ConvertT<char>(neighbor->GetName().at(1))) == oxygen_index))
                        {
                            atom2 = neighbor;
                            break;
                        }
                    }
                }
                if(atom2 != NULL)
                {
                    AtomVector atom2_neighbors = atom2->GetNode()->GetNodeNeighbors();
                    for(AtomVector::iterator it = atom2_neighbors.begin(); it != atom2_neighbors.end(); it++)
                    {
                        MolecularModeling::Atom* neighbor = *it;
                        if(neighbor->GetId().compare(oxygen->GetId()) != 0)
                        {
                            if(neighbor->GetName().at(0) == 'H' &&
                                    (neighbor->GetName().size() > 1 && isdigit(neighbor->GetName().at(1)) &&
                                     gmml::ConvertString<int>(gmml::ConvertT<char>(neighbor->GetName().at(1))) == oxygen_index))
                            {
                                atom1 = neighbor;
                                break;
                            }
                        }
                    }
                    for(AtomVector::iterator it = atom2_neighbors.begin(); it != atom2_neighbors.end(); it++)
                    {
                        MolecularModeling::Atom* neighbor = *it;
                        if(neighbor->GetId().compare(oxygen->GetId()) != 0)
                        {
                            if(neighbor->GetName().at(0) == 'C' &&
                                    (neighbor->GetName().size() > 1 && isdigit(neighbor->GetName().at(1)) &&
                                     gmml::ConvertString<int>(gmml::ConvertT<char>(neighbor->GetName().at(1))) == oxygen_index - 1))
                            {
                                atom3 = neighbor;
                                break;
                            }
                        }
                    }
                }
                if(atom3 != NULL)
                {
                    AtomVector atom3_neighbors = atom3->GetNode()->GetNodeNeighbors();
                    for(AtomVector::iterator it = atom3_neighbors.begin(); it != atom3_neighbors.end(); it++)
                    {
                        MolecularModeling::Atom* neighbor = *it;
                        if(neighbor->GetId().compare(atom2->GetId()) != 0)
                        {
                            if(neighbor->GetName().at(0) == 'H' &&
                                    (neighbor->GetName().size() > 1 && isdigit(neighbor->GetName().at(1)) &&
                                     gmml::ConvertString<int>(gmml::ConvertT<char>(neighbor->GetName().at(1))) == oxygen_index - 1))
                            {
                                atom4 = neighbor;
                                break;
                            }
                        }
                    }
                }
                if(atom1 != NULL && atom2 != NULL && atom3 != NULL && atom4 != NULL)
                {
                    this->SetDihedral(atom1, atom2, atom3, atom4, torsion);
                }
            }
            if(oxygen_index == 8 && type == 9)
            {
                MolecularModeling::Atom* atom1 = NULL;
                MolecularModeling::Atom* atom2 = NULL;
                MolecularModeling::Atom* atom3 = NULL;
                MolecularModeling::Atom* atom4 = NULL;

                AtomVector oxygen_neighbors = oxygen->GetNode()->GetNodeNeighbors();
                for(AtomVector::iterator it = oxygen_neighbors.begin(); it != oxygen_neighbors.end(); it++)
                {
                    MolecularModeling::Atom* neighbor = *it;
                    if(neighbor->GetId().compare(carbon->GetId()) != 0)
                    {
                        if(neighbor->GetName().at(0) == 'C' &&
                                (neighbor->GetName().size() > 1 && isdigit(neighbor->GetName().at(1)) &&
                                 gmml::ConvertString<int>(gmml::ConvertT<char>(neighbor->GetName().at(1))) == oxygen_index))
                        {
                            atom3 = neighbor;
                            break;
                        }
                    }
                }
                if(atom3 != NULL)
                {
                    AtomVector atom3_neighbors = atom3->GetNode()->GetNodeNeighbors();
                    for(AtomVector::iterator it = atom3_neighbors.begin(); it != atom3_neighbors.end(); it++)
                    {
                        MolecularModeling::Atom* neighbor = *it;
                        if(neighbor->GetId().compare(oxygen->GetId()) != 0)
                        {
                            if(neighbor->GetName().at(0) == 'C' &&
                                    (neighbor->GetName().size() > 1 && isdigit(neighbor->GetName().at(1)) &&
                                     gmml::ConvertString<int>(gmml::ConvertT<char>(neighbor->GetName().at(1))) == oxygen_index - 1))
                            {
                                atom4 = neighbor;
                                break;
                            }
                        }
                    }
                    for(AtomVector::iterator it = atom3_neighbors.begin(); it != atom3_neighbors.end(); it++)
                    {
                        MolecularModeling::Atom* neighbor = *it;
                        if(neighbor->GetId().compare(oxygen->GetId()) != 0)
                        {
                            if(neighbor->GetName().at(0) == 'C' &&
                                    (neighbor->GetName().size() > 1 && isdigit(neighbor->GetName().at(1)) &&
                                     gmml::ConvertString<int>(gmml::ConvertT<char>(neighbor->GetName().at(1))) == oxygen_index + 1))
                            {
                                atom2 = neighbor;
                                break;
                            }
                        }
                    }
                }
                if(atom2 != NULL)
                {
                    AtomVector atom2_neighbors = atom2->GetNode()->GetNodeNeighbors();
                    for(AtomVector::iterator it = atom2_neighbors.begin(); it != atom2_neighbors.end(); it++)
                    {
                        MolecularModeling::Atom* neighbor = *it;
                        if(neighbor->GetId().compare(atom3->GetId()) != 0)
                        {
                            if(neighbor->GetName().at(0) == 'O' &&
                                    (neighbor->GetName().size() > 1 && isdigit(neighbor->GetName().at(1)) &&
                                     gmml::ConvertString<int>(gmml::ConvertT<char>(neighbor->GetName().at(1))) == oxygen_index + 1))
                            {
                                atom1 = neighbor;
                                break;
                            }
                        }
                    }
                }
                if(atom1 != NULL && atom2 != NULL && atom3 != NULL && atom4 != NULL)
                {
                    this->SetDihedral(atom1, atom2, atom3, atom4, torsion);
                }
            }
        }
    }
}

void Assembly::SetOmegaDerivativeTorsion(MolecularModeling::Residue *residue, MolecularModeling::Residue *parent_residue, int branch_index, double torsion)
{
    std::string residue_name = residue->GetName();
    MolecularModeling::Atom* parent_target_atom = parent_residue->GetTailAtoms().at(branch_index);

    MolecularModeling::Atom* oxygen = parent_target_atom; ///The oxygen atom of the parent residue that is attached to the new residue
    if(oxygen != NULL)
    {
        int oxygen_index = 1;
        if(oxygen->GetName().size() > 1 && isdigit(oxygen->GetName().at(1)))
            oxygen_index = gmml::ConvertString<int>(gmml::ConvertT<char>(oxygen->GetName().at(1)));
        if(residue_name.compare("SO3") == 0 && oxygen_index == 6)
        {
            MolecularModeling::Atom* atom1 = NULL;
            MolecularModeling::Atom* atom2 = NULL;
            MolecularModeling::Atom* atom3 = NULL;
            MolecularModeling::Atom* atom4 = oxygen;

            AtomVector oxygen_neighbors = oxygen->GetNode()->GetNodeNeighbors();
            for(AtomVector::iterator it = oxygen_neighbors.begin(); it != oxygen_neighbors.end(); it++)
            {
                MolecularModeling::Atom* neighbor = *it;
                if(neighbor->GetName().at(0) == 'C' &&
                        (neighbor->GetName().size() > 1 && isdigit(neighbor->GetName().at(1)) &&
                         gmml::ConvertString<int>(gmml::ConvertT<char>(neighbor->GetName().at(1))) == oxygen_index))
                {
                    atom3 = neighbor;
                    break;
                }

            }
            if(atom3 != NULL)
            {
                AtomVector atom2_neighbors = atom3->GetNode()->GetNodeNeighbors();
                for(AtomVector::iterator it = atom2_neighbors.begin(); it != atom2_neighbors.end(); it++)
                {
                    MolecularModeling::Atom* neighbor = *it;
                    if(neighbor->GetId().compare(oxygen->GetId()) != 0)
                    {
                        if(neighbor->GetName().at(0) == 'C' &&
                                (neighbor->GetName().size() > 1 && isdigit(neighbor->GetName().at(1)) &&
                                 gmml::ConvertString<int>(gmml::ConvertT<char>(neighbor->GetName().at(1))) == oxygen_index - 1))
                        {
                            atom2 = neighbor;
                            break;
                        }
                    }
                }
            }
            if(atom2 != NULL)
            {
                AtomVector atom3_neighbors = atom2->GetNode()->GetNodeNeighbors();
                for(AtomVector::iterator it = atom3_neighbors.begin(); it != atom3_neighbors.end(); it++)
                {
                    MolecularModeling::Atom* neighbor = *it;
                    if(neighbor->GetId().compare(atom3->GetId()) != 0)
                    {
                        if(neighbor->GetName().at(0) == 'O' &&
                                (neighbor->GetName().size() > 1 && isdigit(neighbor->GetName().at(1)) &&
                                 gmml::ConvertString<int>(gmml::ConvertT<char>(neighbor->GetName().at(1))) == oxygen_index - 1))
                        {
                            atom1 = neighbor;
                            break;
                        }
                    }
                }
            }
            if(atom1 != NULL && atom2 != NULL && atom3 != NULL && atom4 != NULL)
            {
                this->SetDihedral(atom1, atom2, atom3, atom4, torsion);
            }
        }
    }
}

// OG July 3rd 2018. This function leaks memory and I call it a lot during gp building. I've created a version that doesn't leak below.
//void Assembly::SetDihedral(MolecularModeling::Atom *atom1, MolecularModeling::Atom *atom2, MolecularModeling::Atom *atom3, MolecularModeling::Atom *atom4, double torsion)
//{
//    double current_dihedral = 0.0;
//    GeometryTopology::Coordinate* a1 = atom1->GetCoordinates().at(model_index_);
//    GeometryTopology::Coordinate* a2 = atom2->GetCoordinates().at(model_index_);
//    GeometryTopology::Coordinate* a3 = atom3->GetCoordinates().at(model_index_);
//    GeometryTopology::Coordinate* a4 = atom4->GetCoordinates().at(model_index_);

//    GeometryTopology::Coordinate* b1 = new GeometryTopology::Coordinate(*a2);
//    b1->operator -(*a1);
//    GeometryTopology::Coordinate* b2 = new GeometryTopology::Coordinate(*a3);
//    b2->operator -(*a2);
//    GeometryTopology::Coordinate* b3 = new GeometryTopology::Coordinate(*a4);
//    b3->operator -(*a3);
//    GeometryTopology::Coordinate* b4 = new GeometryTopology::Coordinate(*b2);
//    b4->operator *(-1);

//    GeometryTopology::Coordinate* b2xb3 = new GeometryTopology::Coordinate(*b2);
//    b2xb3->CrossProduct(*b3);

//    GeometryTopology::Coordinate* b1_m_b2n = new GeometryTopology::Coordinate(*b1);
//    b1_m_b2n->operator *(b2->length());

//    GeometryTopology::Coordinate* b1xb2 = new GeometryTopology::Coordinate(*b1);
//    b1xb2->CrossProduct(*b2);

//    current_dihedral = atan2(b1_m_b2n->DotProduct(*b2xb3), b1xb2->DotProduct(*b2xb3));

//    double** torsion_matrix = gmml::GenerateRotationMatrix(b4, a2, current_dihedral - gmml::ConvertDegree2Radian(torsion));

//    AtomVector atomsToRotate = AtomVector();
//    atomsToRotate.push_back(atom2);
//    atom3->FindConnectedAtoms(atomsToRotate);
//   // std::cout << "Moving ";
//    for(AtomVector::iterator it = atomsToRotate.begin(); it != atomsToRotate.end(); it++)
//    {
//   //     std::cout << (*it)->GetName() << ", ";
//        GeometryTopology::Coordinate* atom_coordinate = (*it)->GetCoordinates().at(model_index_);
//        GeometryTopology::Coordinate* result = new GeometryTopology::Coordinate();
//        result->SetX(torsion_matrix[0][0] * atom_coordinate->GetX() + torsion_matrix[0][1] * atom_coordinate->GetY() +
//                torsion_matrix[0][2] * atom_coordinate->GetZ() + torsion_matrix[0][3]);
//        result->SetY(torsion_matrix[1][0] * atom_coordinate->GetX() + torsion_matrix[1][1] * atom_coordinate->GetY() +
//                torsion_matrix[1][2] * atom_coordinate->GetZ() + torsion_matrix[1][3]);
//        result->SetZ(torsion_matrix[2][0] * atom_coordinate->GetX() + torsion_matrix[2][1] * atom_coordinate->GetY() +
//                torsion_matrix[2][2] * atom_coordinate->GetZ() + torsion_matrix[2][3]);

//        (*it)->GetCoordinates().at(model_index_)->SetX(result->GetX());
//        (*it)->GetCoordinates().at(model_index_)->SetY(result->GetY());
//        (*it)->GetCoordinates().at(model_index_)->SetZ(result->GetZ());
//    }
// //   std::cout << "\n";
//}

// Changed to not leak memory by Oliver on July 3rd, 2018. See above for old code.
double Assembly::GetDihedral(MolecularModeling::Atom *atom1, MolecularModeling::Atom *atom2, MolecularModeling::Atom *atom3, MolecularModeling::Atom *atom4)
{
    double current_dihedral = 0.0;
    GeometryTopology::Coordinate* a1 = atom1->GetCoordinates().at(0);
    GeometryTopology::Coordinate* a2 = atom2->GetCoordinates().at(0);
    GeometryTopology::Coordinate* a3 = atom3->GetCoordinates().at(0);
    GeometryTopology::Coordinate* a4 = atom4->GetCoordinates().at(0);

    GeometryTopology::Coordinate b1 = a2;
    b1.operator -(*a1);
    GeometryTopology::Coordinate b2 = a3;
    b2.operator -(*a2);
    GeometryTopology::Coordinate b3 = a4;
    b3.operator -(*a3);
    GeometryTopology::Coordinate b4 = b2;
    b4.operator *(-1);

    GeometryTopology::Coordinate b2xb3 = b2;
    b2xb3.CrossProduct(b3);

    GeometryTopology::Coordinate b1_m_b2n = b1;
    b1_m_b2n.operator *(b2.length());

    GeometryTopology::Coordinate b1xb2 = b1;
    b1xb2.CrossProduct(b2);

    current_dihedral = atan2(b1_m_b2n.DotProduct(b2xb3), b1xb2.DotProduct(b2xb3));
    return current_dihedral /3.1415 * 180;
}

void Assembly::SetDihedral(MolecularModeling::Atom *atom1, MolecularModeling::Atom *atom2, MolecularModeling::Atom *atom3, MolecularModeling::Atom *atom4, double torsion)
{
    double current_dihedral = 0.0;
    GeometryTopology::Coordinate* a1 = atom1->GetCoordinates().at(model_index_);
    GeometryTopology::Coordinate* a2 = atom2->GetCoordinates().at(model_index_);
    GeometryTopology::Coordinate* a3 = atom3->GetCoordinates().at(model_index_);
    GeometryTopology::Coordinate* a4 = atom4->GetCoordinates().at(model_index_);

    GeometryTopology::Coordinate b1 = a2;
    b1.operator -(*a1);
    GeometryTopology::Coordinate b2 = a3;
    b2.operator -(*a2);
    GeometryTopology::Coordinate b3 = a4;
    b3.operator -(*a3);
    GeometryTopology::Coordinate b4 = b2;
    b4.operator *(-1);

    GeometryTopology::Coordinate b2xb3 = b2;
    b2xb3.CrossProduct(b3);

    GeometryTopology::Coordinate b1_m_b2n = b1;
    b1_m_b2n.operator *(b2.length());

    GeometryTopology::Coordinate b1xb2 = b1;
    b1xb2.CrossProduct(b2);

    current_dihedral = atan2(b1_m_b2n.DotProduct(b2xb3), b1xb2.DotProduct(b2xb3));

    double** torsion_matrix = gmml::GenerateRotationMatrix(&b4, a2, current_dihedral - gmml::ConvertDegree2Radian(torsion));

    AtomVector atomsToRotate = AtomVector();
    atomsToRotate.push_back(atom2);
    atom3->FindConnectedAtoms(atomsToRotate);
    for(AtomVector::iterator it = atomsToRotate.begin(); it != atomsToRotate.end(); it++)
    {
   //     std::cout << (*it)->GetName() << ", ";
        GeometryTopology::Coordinate* atom_coordinate = (*it)->GetCoordinates().at(model_index_);
        GeometryTopology::Coordinate result;
        result.SetX(torsion_matrix[0][0] * atom_coordinate->GetX() + torsion_matrix[0][1] * atom_coordinate->GetY() +
                torsion_matrix[0][2] * atom_coordinate->GetZ() + torsion_matrix[0][3]);
        result.SetY(torsion_matrix[1][0] * atom_coordinate->GetX() + torsion_matrix[1][1] * atom_coordinate->GetY() +
                torsion_matrix[1][2] * atom_coordinate->GetZ() + torsion_matrix[1][3]);
        result.SetZ(torsion_matrix[2][0] * atom_coordinate->GetX() + torsion_matrix[2][1] * atom_coordinate->GetY() +
                torsion_matrix[2][2] * atom_coordinate->GetZ() + torsion_matrix[2][3]);

        (*it)->GetCoordinates().at(model_index_)->SetX(result.GetX());
        (*it)->GetCoordinates().at(model_index_)->SetY(result.GetY());
        (*it)->GetCoordinates().at(model_index_)->SetZ(result.GetZ());
    }
    return;
}

void Assembly::SetAngle(MolecularModeling::Atom* atom1, MolecularModeling::Atom* atom2, MolecularModeling::Atom* atom3, double angle)
{
    double current_angle = 0.0;
    GeometryTopology::Coordinate* a1 = atom1->GetCoordinates().at(model_index_);
    GeometryTopology::Coordinate* a2 = atom2->GetCoordinates().at(model_index_);
    GeometryTopology::Coordinate* a3 = atom3->GetCoordinates().at(model_index_);

    GeometryTopology::Coordinate* b1 = new GeometryTopology::Coordinate(*a1);
    b1->operator -(*a2);
    GeometryTopology::Coordinate* b2 = new GeometryTopology::Coordinate(*a3);
    b2->operator -(*a2);

    current_angle = acos((b1->DotProduct(*b2)) / (b1->length() * b2->length() + gmml::DIST_EPSILON));
    double rotation_angle = gmml::ConvertDegree2Radian(angle) - current_angle;

    GeometryTopology::Coordinate* direction = new GeometryTopology::Coordinate(*b1);
    direction->CrossProduct(*b2);
    direction->Normalize();
    double** rotation_matrix = gmml::GenerateRotationMatrix(direction, a2, rotation_angle);

    AtomVector atomsToRotate = AtomVector();
    atomsToRotate.push_back(atom2);
    atom3->FindConnectedAtoms(atomsToRotate);

    for(AtomVector::iterator it = atomsToRotate.begin() + 1; it != atomsToRotate.end(); it++)
    {
        GeometryTopology::Coordinate* atom_coordinate = (*it)->GetCoordinates().at(model_index_);
        GeometryTopology::Coordinate* result = new GeometryTopology::Coordinate();
        result->SetX(rotation_matrix[0][0] * atom_coordinate->GetX() + rotation_matrix[0][1] * atom_coordinate->GetY() +
                rotation_matrix[0][2] * atom_coordinate->GetZ() + rotation_matrix[0][3]);
        result->SetY(rotation_matrix[1][0] * atom_coordinate->GetX() + rotation_matrix[1][1] * atom_coordinate->GetY() +
                rotation_matrix[1][2] * atom_coordinate->GetZ() + rotation_matrix[1][3]);
        result->SetZ(rotation_matrix[2][0] * atom_coordinate->GetX() + rotation_matrix[2][1] * atom_coordinate->GetY() +
                rotation_matrix[2][2] * atom_coordinate->GetZ() + rotation_matrix[2][3]);

        (*it)->GetCoordinates().at(model_index_)->SetX(result->GetX());
        (*it)->GetCoordinates().at(model_index_)->SetY(result->GetY());
        (*it)->GetCoordinates().at(model_index_)->SetZ(result->GetZ());
    }
}

double Assembly::CalculateBondAngleByCoordinates(GeometryTopology::Coordinate* atom1_crd, GeometryTopology::Coordinate* atom2_crd, GeometryTopology::Coordinate* atom3_crd)
{
    GeometryTopology::Coordinate* b1 = new GeometryTopology::Coordinate(*atom1_crd);
    b1->operator -(*atom2_crd);
    GeometryTopology::Coordinate* b2 = new GeometryTopology::Coordinate(*atom3_crd);
    b2->operator -(*atom2_crd);

    return acos(b1->DotProduct((*b2)) / b1->length() / b2->length());

}

double Assembly::CalculateBondAngleByAtoms(MolecularModeling::Atom *atom1, MolecularModeling::Atom *atom2, MolecularModeling::Atom *atom3)
{
    GeometryTopology::Coordinate* a1 = atom1->GetCoordinates().at(model_index_);
    GeometryTopology::Coordinate* a2 = atom2->GetCoordinates().at(model_index_);
    GeometryTopology::Coordinate* a3 = atom3->GetCoordinates().at(model_index_);

    GeometryTopology::Coordinate* b1 = new GeometryTopology::Coordinate(*a1);
    b1->operator -(*a2);
    GeometryTopology::Coordinate* b2 = new GeometryTopology::Coordinate(*a3);
    b2->operator -(*a2);

    return acos(b1->DotProduct((*b2)) / b1->length() / b2->length());
}

double Assembly::CalculateTorsionAngleByCoordinates(GeometryTopology::Coordinate* atom1_crd, GeometryTopology::Coordinate* atom2_crd, GeometryTopology::Coordinate* atom3_crd, GeometryTopology::Coordinate* atom4_crd) //This defines angles in radians!!!!!!!!
{
    double current_dihedral = 0.0;

    GeometryTopology::Coordinate* b1 = new GeometryTopology::Coordinate(*atom2_crd);
    b1->operator -(*atom1_crd);
    GeometryTopology::Coordinate* b2 = new GeometryTopology::Coordinate(*atom3_crd);
    b2->operator -(*atom2_crd);
    GeometryTopology::Coordinate* b3 = new GeometryTopology::Coordinate(*atom4_crd);
    b3->operator -(*atom3_crd);
    GeometryTopology::Coordinate* b4 = new GeometryTopology::Coordinate(*b2);
    b4->operator *(-1);

    GeometryTopology::Coordinate* b2xb3 = new GeometryTopology::Coordinate(*b2);
    b2xb3->CrossProduct(*b3);

    GeometryTopology::Coordinate* b1_m_b2n = new GeometryTopology::Coordinate(*b1);
    b1_m_b2n->operator *(b2->length());

    GeometryTopology::Coordinate* b1xb2 = new GeometryTopology::Coordinate(*b1);
    b1xb2->CrossProduct(*b2);

    current_dihedral = atan2(b1_m_b2n->DotProduct(*b2xb3), b1xb2->DotProduct(*b2xb3));
    return current_dihedral;
}


// atom1Coordinates << atom1->GetId() << "X: " << a1->GetX() << " Y: " << a1->GetY() << " Z: " << a1->GetZ();
// gmml::log(__LINE__, __FILE__, gmml::INF, atom1Coordinates.str());
double Assembly::CalculateTorsionAngleByAtoms(MolecularModeling::Atom *atom1, MolecularModeling::Atom *atom2, MolecularModeling::Atom *atom3, MolecularModeling::Atom *atom4)
{
  
  double current_dihedral = 0.0;
  GeometryTopology::Coordinate* a1 = atom1->GetCoordinates().at(model_index_);
  GeometryTopology::Coordinate* a2 = atom2->GetCoordinates().at(model_index_);
  GeometryTopology::Coordinate* a3 = atom3->GetCoordinates().at(model_index_);
  GeometryTopology::Coordinate* a4 = atom4->GetCoordinates().at(model_index_);

  GeometryTopology::Coordinate b1 = a2;
  b1.operator -(*a1);
  GeometryTopology::Coordinate b2 = a3;
  b2.operator -(*a2);
  GeometryTopology::Coordinate b3 = a4;
  b3.operator -(*a3);
  GeometryTopology::Coordinate b4 = b2;
  b4.operator *(-1);

  GeometryTopology::Coordinate b2xb3 = b2;
  b2xb3.CrossProduct(b3);

  GeometryTopology::Coordinate b1_m_b2n = b1;
  b1_m_b2n.operator *(b2.length());

  GeometryTopology::Coordinate b1xb2 = b1;
  b1xb2.CrossProduct(b2);

  current_dihedral = atan2(b1_m_b2n.DotProduct(b2xb3), b1xb2.DotProduct(b2xb3));
    // 
    // 
    // double current_dihedral = 0.0;
    // GeometryTopology::Coordinate* a1 = atom1->GetCoordinates().at(model_index_);
    // std::stringstream atom1Coordinates;
    // GeometryTopology::Coordinate* a2 = atom2->GetCoordinates().at(model_index_);
    // GeometryTopology::Coordinate* a3 = atom3->GetCoordinates().at(model_index_);
    // GeometryTopology::Coordinate* a4 = atom4->GetCoordinates().at(model_index_);
    // 
    // GeometryTopology::Coordinate* b1 = new GeometryTopology::Coordinate(*a2);
    // b1->operator -(*a1);
    // GeometryTopology::Coordinate* b2 = new GeometryTopology::Coordinate(*a3);
    // b2->operator -(*a2);
    // GeometryTopology::Coordinate* b3 = new GeometryTopology::Coordinate(*a4);
    // b3->operator -(*a3);
    // GeometryTopology::Coordinate* b4 = new GeometryTopology::Coordinate(*b2);
    // b4->operator *(-1);
    // 
    // GeometryTopology::Coordinate* b2xb3 = new GeometryTopology::Coordinate(*b2);
    // b2xb3->CrossProduct(*b3);
    // 
    // GeometryTopology::Coordinate* b1_m_b2n = new GeometryTopology::Coordinate(*b1);
    // b1_m_b2n->operator *(b2->length());
    // 
    // GeometryTopology::Coordinate* b1xb2 = new GeometryTopology::Coordinate(*b1);
    // b1xb2->CrossProduct(*b2);
    // 
    // current_dihedral = atan2(b1_m_b2n->DotProduct(*b2xb3), b1xb2->DotProduct(*b2xb3));
    return current_dihedral;
}

void Assembly::GetCenterOfMass(GeometryTopology::Coordinate *center_of_mass)
{
    //    center_of_mass = new GeometryTopology::Coordinate();
    AtomVector all_atoms_of_assembly = this->GetAllAtomsOfAssembly();
    for(AtomVector::iterator it = all_atoms_of_assembly.begin(); it != all_atoms_of_assembly.end(); it++)
    {
        MolecularModeling::Atom* atom = *it;
        center_of_mass->Translate(atom->MolecularDynamicAtom::GetMass() * atom->GetCoordinates().at(model_index_)->GetX(),
                                  atom->MolecularDynamicAtom::GetMass() * atom->GetCoordinates().at(model_index_)->GetY(),
                                  atom->MolecularDynamicAtom::GetMass() * atom->GetCoordinates().at(model_index_)->GetZ());
    }
    double total_mass = this->GetTotalMass();
    center_of_mass->operator /(GeometryTopology::Coordinate(center_of_mass->GetX() / total_mass,
                                          center_of_mass->GetY() / total_mass,
                                          center_of_mass->GetZ() / total_mass));
}

void Assembly::GetCenterOfGeometry(GeometryTopology::Coordinate *center_of_geometry)
{
    double sumX = 0.0;
    double sumY = 0.0;
    double sumZ = 0.0;
    GeometryTopology::CoordinateVector all_coords = this->GetAllCoordinates();
    for(GeometryTopology::CoordinateVector::iterator it = all_coords.begin(); it != all_coords.end(); it++)
    {
        GeometryTopology::Coordinate coord = *it;
        sumX += coord.GetX();
        sumY += coord.GetY();
        sumZ += coord.GetZ();
    }
    center_of_geometry->SetX( (sumX / all_coords.size()) );
    center_of_geometry->SetY( (sumY / all_coords.size()) );
    center_of_geometry->SetZ( (sumZ / all_coords.size()) );
}

void Assembly::GetBoundary(GeometryTopology::Coordinate* lower_left_back_corner, GeometryTopology::Coordinate* upper_right_front_corner)
{
    //    lower_left_back_corner = new GeometryTopology::Coordinate(-INFINITY, -INFINITY, -INFINITY);
    //    upper_right_front_corner = new GeometryTopology::Coordinate(INFINITY, INFINITY, INFINITY);
    lower_left_back_corner->SetX(INFINITY);
    lower_left_back_corner->SetY(INFINITY);
    lower_left_back_corner->SetZ(INFINITY);
    upper_right_front_corner->SetX(-INFINITY);
    upper_right_front_corner->SetY(-INFINITY);
    upper_right_front_corner->SetZ(-INFINITY);
    AtomVector all_atoms_of_assembly = this->GetAllAtomsOfAssembly();
    for(AtomVector::iterator it = all_atoms_of_assembly.begin(); it != all_atoms_of_assembly.end(); it++)
    {
        MolecularModeling::Atom* atom = *it;
        if(atom->MolecularDynamicAtom::GetRadius() == gmml::dNotSet)
        {
            //            gmml::log(__LINE__, __FILE__,  gmml::ERR, "There is no information of the atom type/radius/charge of the atoms in the given library/parameter file");
            //            cout << "There is no information of the atom type/radius/charge of the atoms in the given library/parameter file" << endl;
            atom->MolecularDynamicAtom::SetRadius(gmml::DEFAULT_RADIUS);
            //            std::stringstream ss;
            //            ss << "The default value has been set for " << atom->GetId();
            //            gmml::log(__LINE__, __FILE__,  gmml::ERR, ss.str());
            //            cout << ss.str() << endl;
            //            return;
        }
        double upper_right_front_x = atom->GetCoordinates().at(model_index_)->GetX() + atom->MolecularDynamicAtom::GetRadius();
        double lower_left_back_x = atom->GetCoordinates().at(model_index_)->GetX() - atom->MolecularDynamicAtom::GetRadius();
        if(upper_right_front_x > upper_right_front_corner->GetX())
            upper_right_front_corner->SetX(upper_right_front_x);
        if(lower_left_back_x < lower_left_back_corner->GetX())
            lower_left_back_corner->SetX(lower_left_back_x);

        double upper_right_front_y = atom->GetCoordinates().at(model_index_)->GetY() + atom->MolecularDynamicAtom::GetRadius();
        double lower_left_back_y = atom->GetCoordinates().at(model_index_)->GetY() - atom->MolecularDynamicAtom::GetRadius();
        if(upper_right_front_y > upper_right_front_corner->GetY())
            upper_right_front_corner->SetY(upper_right_front_y);
        if(lower_left_back_y < lower_left_back_corner->GetY())
            lower_left_back_corner->SetY(lower_left_back_y);

        double upper_right_front_z = atom->GetCoordinates().at(model_index_)->GetZ() + atom->MolecularDynamicAtom::GetRadius();
        double lower_left_back_z = atom->GetCoordinates().at(model_index_)->GetZ() - atom->MolecularDynamicAtom::GetRadius();
        if(upper_right_front_z > upper_right_front_corner->GetZ())
            upper_right_front_corner->SetZ(upper_right_front_z);
        if(lower_left_back_z < lower_left_back_corner->GetZ())
            lower_left_back_corner->SetZ(lower_left_back_z);
    }
    if(all_atoms_of_assembly.size() == 0)
    {
        lower_left_back_corner->SetX(0.0);
        lower_left_back_corner->SetY(0.0);
        lower_left_back_corner->SetZ(0.0);
        upper_right_front_corner->SetX(0.0);
        upper_right_front_corner->SetY(0.0);
        upper_right_front_corner->SetZ(0.0);
    }
}

std::string Assembly::CalculateRSOrientations(MolecularModeling::Atom *prev_atom, MolecularModeling::Atom *target, MolecularModeling::Atom *next_atom)
{
    std::string orientation = "";
    ///Calculating the plane based on the two ring neighbors of the current atom
    GeometryTopology::Coordinate prev_atom_coord = GeometryTopology::Coordinate(*prev_atom->GetCoordinates().at(model_index_));
    GeometryTopology::Coordinate current_atom_coord = GeometryTopology::Coordinate(*target->GetCoordinates().at(model_index_));
    GeometryTopology::Coordinate next_atom_coord = GeometryTopology::Coordinate(*next_atom->GetCoordinates().at(model_index_));
    prev_atom_coord.operator -(current_atom_coord) ;
    next_atom_coord.operator -(current_atom_coord) ;
    GeometryTopology::Plane plane = GeometryTopology::Plane();
    plane.SetV1(prev_atom_coord);
    plane.SetV2(next_atom_coord);
    GeometryTopology::Coordinate normal_v = plane.GetUnitNormalVector();

    MolecularModeling::AtomNode* node = target->GetNode();
    AtomVector neighbors = node->GetNodeNeighbors();
    for(AtomVector::iterator it1 = neighbors.begin(); it1 != neighbors.end(); it1++)
    {
        MolecularModeling::Atom* neighbor = (*it1);
        if(neighbor->GetId().at(0) == 'O')
        {
            GeometryTopology::Coordinate side_atom_coord = GeometryTopology::Coordinate(*neighbor->GetCoordinates().at(model_index_));
            side_atom_coord.operator -(current_atom_coord);
            side_atom_coord.Normalize();
            GeometryTopology::Coordinate normal_v_x_side_atom = normal_v;
            normal_v_x_side_atom.CrossProduct(side_atom_coord); ///cross product (perpendicular vector to plan's normal vector and the normal vector of side atom oxygen)
            double sin_theta = normal_v_x_side_atom.length()/(normal_v.length()*side_atom_coord.length());

            if(sin_theta >= 0) ///theta between plan's normal vector and the vector of side oxygen is between 0 to 180 degree
                orientation = "R";
            else ///theta between plan's normal vector and the vector of side oxygen is between grater than 180 degree
                orientation = "S";
            break;
        }
    }
    return orientation;
}

double Assembly::CalculateAtomicOverlaps(Assembly *assemblyB)
{
    AtomVector assemblyBAtoms = assemblyB->GetAllAtomsOfAssembly();
    AtomVector assemblyAAtoms = this->GetAllAtomsOfAssembly();
    return gmml::CalculateAtomicOverlaps(assemblyAAtoms, assemblyBAtoms);
}
double Assembly::CalculateAtomicOverlaps(AtomVector assemblyBAtoms)
{
    AtomVector assemblyAAtoms = this->GetAllAtomsOfAssembly();
    return gmml::CalculateAtomicOverlaps(assemblyAAtoms, assemblyBAtoms);
}




AtomVector Assembly::GetAllAtomsOfAssemblyWithinXAngstromOf(GeometryTopology::Coordinate *coordinate, double distance)
{
    AtomVector returnAtoms;
    AtomVector assemblyAtoms = this->GetAllAtomsOfAssembly();
    for(AtomVector::iterator it = assemblyAtoms.begin(); it != assemblyAtoms.end(); ++it)
    {
        MolecularModeling::Atom *atom = *it;
        if (atom->GetDistanceToCoordinate(coordinate) <= distance)
        {
            returnAtoms.push_back(atom);
        }
    }
    return returnAtoms;
}
