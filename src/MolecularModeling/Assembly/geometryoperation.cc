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
void Assembly::SetDerivativeAngle(Residue *residue, Residue *parent_residue, int branch_index)
{
    Atom* atom3 = residue->GetHeadAtoms().at(0);
    Atom* atom2 = parent_residue->GetTailAtoms().at(branch_index);
    Atom* atom1 = NULL;
    if(atom2 != NULL)
    {
        int atom2_index = 1;
        if(atom2->GetName().size() > 1 && isdigit(atom2->GetName().at(1)))
            atom2_index = ConvertString<int>(ConvertT<char>(atom2->GetName().at(1)));

        AtomVector atom2_neighbors = atom2->GetNode()->GetNodeNeighbors();
        for(AtomVector::iterator it = atom2_neighbors.begin(); it != atom2_neighbors.end(); it++)
        {
            Atom* neighbor = *it;
            if(neighbor->GetId().compare(atom3->GetId()) != 0)
            {
                if(neighbor->GetName().at(0) == 'C' &&
                        (neighbor->GetName().size() > 1 && isdigit(neighbor->GetName().at(1)) &&
                         ConvertString<int>(ConvertT<char>(neighbor->GetName().at(1))) == atom2_index))
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

void Assembly::SetAttachedResidueBond(Residue *residue, Residue *parent_residue, int branch_index, string parameter_file)
{
    ParameterFile* parameter = new ParameterFile(parameter_file);
    ParameterFile::BondMap parameter_bonds = parameter->GetBonds();
    Atom* residue_head_atom = residue->GetHeadAtoms().at(0);
    Atom* parent_target_atom = parent_residue->GetTailAtoms().at(branch_index);
    AtomVector residue_head_atom_adjacent_atoms = AtomVector();
    AtomVector parent_target_atom_adjacent_atoms = AtomVector();
    double bond_length = BOND_LENGTH;
    vector<string> bond = vector<string>();
    bond.push_back(residue_head_atom->GetAtomType());
    bond.push_back(parent_target_atom->GetAtomType());
    vector<string> reverse_bond = vector<string>();
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

    Coordinate* residue_direction = new Coordinate();
    for(AtomVector::iterator it = residue_head_atom_adjacent_atoms.begin(); it != residue_head_atom_adjacent_atoms.end(); it++)
    {
        Atom* atom = *it;
        if(atom->GetId().compare(parent_target_atom->GetId()) != 0)
        {
            Coordinate* dist = new Coordinate(*residue_head_atom->GetCoordinates().at(model_index_));
            dist->operator -(*atom->GetCoordinates().at(model_index_));
            dist->Normalize();
            residue_direction->operator +(*dist);
            residue_direction->Normalize();
        }
    }

    residue_direction->Normalize();
    residue_direction->operator *(bond_length);

    residue_direction->operator +(*residue_head_atom->GetCoordinates().at(model_index_));

    Coordinate* oxygen_position = new Coordinate(residue_direction->GetX(), residue_direction->GetY(), residue_direction->GetZ());
    Coordinate* offset = new Coordinate(*parent_target_atom->GetCoordinates().at(model_index_));
    offset->operator -(*oxygen_position);

    AtomVector atomsToTranslate = AtomVector();
    atomsToTranslate.push_back(parent_target_atom);
    residue_head_atom->FindConnectedAtoms(atomsToTranslate);

    for(AtomVector::iterator it = atomsToTranslate.begin() + 1; it != atomsToTranslate.end(); it++)
    {
        (*it)->GetCoordinates().at(model_index_)->Translate(offset->GetX(), offset->GetY(), offset->GetZ());
    }
}

void Assembly::SetAttachedResidueAngle(Residue *residue, Residue *parent_residue, int branch_index, string parameter_file)
{
    ParameterFile* parameter = new ParameterFile(parameter_file);
    ParameterFile::BondMap parameter_bonds = parameter->GetBonds();
    Atom* residue_head_atom = residue->GetHeadAtoms().at(0);
    Atom* parent_target_atom = parent_residue->GetTailAtoms().at(branch_index);
    AtomVector residue_head_atom_adjacent_atoms = AtomVector();
    AtomVector parent_target_atom_adjacent_atoms = AtomVector();

    double bond_length = BOND_LENGTH;
    vector<string> bond = vector<string>();
    bond.push_back(residue_head_atom->GetAtomType());
    bond.push_back(parent_target_atom->GetAtomType());
    vector<string> reverse_bond = vector<string>();
    reverse_bond.push_back(parent_target_atom->GetAtomType());
    reverse_bond.push_back(residue_head_atom->GetAtomType());
    if(parameter_bonds.find(bond) != parameter_bonds.end())
        bond_length = parameter_bonds[bond]->GetLength();
    else if(parameter_bonds.find(reverse_bond) != parameter_bonds.end())
        bond_length = parameter_bonds[reverse_bond]->GetLength();
    residue_head_atom_adjacent_atoms = residue_head_atom->GetNode()->GetNodeNeighbors();
    parent_target_atom_adjacent_atoms = parent_target_atom->GetNode()->GetNodeNeighbors();

    Coordinate* carbon_direction = new Coordinate();
    for(AtomVector::iterator it = parent_target_atom_adjacent_atoms.begin(); it != parent_target_atom_adjacent_atoms.end(); it++)
    {
        Atom* atom = *it;
        if(atom->GetId().compare(residue_head_atom->GetId()) != 0)
        {
            Coordinate* dist = new Coordinate(*parent_target_atom->GetCoordinates().at(model_index_));
            dist->operator -(*atom->GetCoordinates().at(model_index_));
            dist->Normalize();
            carbon_direction->operator +(*dist);
            carbon_direction->Normalize();
        }
    }

    carbon_direction->Normalize();
    carbon_direction->operator *(bond_length);
    carbon_direction->operator +(*parent_target_atom->GetCoordinates().at(model_index_));

    Coordinate* carbon_position = new Coordinate(*carbon_direction);

    Coordinate* carbon_target = new Coordinate(*carbon_position);
    carbon_target->operator -(*parent_target_atom->GetCoordinates().at(model_index_));

    Coordinate* head_target = new Coordinate(*residue_head_atom->GetCoordinates().at(model_index_));
    head_target->operator -(*parent_target_atom->GetCoordinates().at(model_index_));

    double angle = acos((carbon_target->DotProduct(*head_target)) / (carbon_target->length() * head_target->length() + DIST_EPSILON));
    double rotation_angle = ConvertDegree2Radian(PI_DEGREE - ROTATION_ANGLE) - angle;

    Coordinate* direction = new Coordinate(*carbon_target);
    direction->CrossProduct(*head_target);
    direction->Normalize();
    double** rotation_matrix = GenerateRotationMatrix(direction, parent_target_atom->GetCoordinates().at(model_index_), rotation_angle);

    AtomVector atomsToRotate = AtomVector();
    atomsToRotate.push_back(parent_target_atom);
    residue_head_atom->FindConnectedAtoms(atomsToRotate);

    for(AtomVector::iterator it = atomsToRotate.begin() + 1; it != atomsToRotate.end(); it++)
    {
        Coordinate* atom_coordinate = (*it)->GetCoordinates().at(model_index_);
        Coordinate* result = new Coordinate();
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

void Assembly::SetAttachedResidueTorsion(Residue *residue, Residue *parent_residue, int branch_index)
{
    Atom* residue_head_atom = residue->GetHeadAtoms().at(0);
    Atom* parent_target_atom = parent_residue->GetTailAtoms().at(branch_index);

    Atom* carbon = residue_head_atom; ///The carbon atom of the new residue that is attached to the parent residue
    if(carbon != NULL)
    {
        Atom* oxygen = parent_target_atom; ///The oxygen atom of the parent residue that is attached to the new residue
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
                oxygen_index = ConvertString<int>(ConvertT<char>(oxygen->GetName().at(1)));

            if(oxygen_index == 5 || oxygen_index == 6)
            {
                Atom* atom1 = NULL;
                Atom* atom2 = NULL;
                Atom* atom3 = oxygen;
                Atom* atom4 = carbon;

                AtomVector oxygen_neighbors = oxygen->GetNode()->GetNodeNeighbors();
                for(AtomVector::iterator it = oxygen_neighbors.begin(); it != oxygen_neighbors.end(); it++)
                {
                    Atom* neighbor = *it;
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
                        Atom* neighbor = *it;
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
                Atom* atom1 = NULL;
                Atom* atom2 = NULL;
                Atom* atom3 = oxygen;
                Atom* atom4 = carbon;

                AtomVector oxygen_neighbors = oxygen->GetNode()->GetNodeNeighbors();
                for(AtomVector::iterator it = oxygen_neighbors.begin(); it != oxygen_neighbors.end(); it++)
                {
                    Atom* neighbor = *it;
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
                        Atom* neighbor = *it;
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
                Atom* atom1 = NULL;
                Atom* atom2 = oxygen;
                Atom* atom3 = carbon;
                Atom* atom4 = NULL;

                AtomVector oxygen_neighbors = oxygen->GetNode()->GetNodeNeighbors();
                for(AtomVector::iterator it = oxygen_neighbors.begin(); it != oxygen_neighbors.end(); it++)
                {
                    Atom* neighbor = *it;
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
                    Atom* neighbor = *it;
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
                Atom* atom1 = NULL;
                Atom* atom2 = oxygen;
                Atom* atom3 = carbon;
                Atom* atom4 = NULL;

                AtomVector oxygen_neighbors = oxygen->GetNode()->GetNodeNeighbors();
                for(AtomVector::iterator it = oxygen_neighbors.begin(); it != oxygen_neighbors.end(); it++)
                {
                    Atom* neighbor = *it;
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
                    Atom* neighbor = *it;
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
                Atom* atom1 = NULL;
                Atom* atom2 = NULL;
                Atom* atom3 = NULL;
                Atom* atom4 = oxygen;

                AtomVector oxygen_neighbors = oxygen->GetNode()->GetNodeNeighbors();
                for(AtomVector::iterator it = oxygen_neighbors.begin(); it != oxygen_neighbors.end(); it++)
                {
                    Atom* neighbor = *it;
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
                        Atom* neighbor = *it;
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
                            Atom* neighbor = *it;
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

void Assembly::SetPhiTorsion(Residue *residue, Residue *parent_residue, int branch_index, double torsion)
{
    Atom* residue_head_atom = residue->GetHeadAtoms().at(0);
    Atom* parent_target_atom = parent_residue->GetTailAtoms().at(branch_index);

    Atom* carbon = residue_head_atom; ///The carbon atom of the new residue that is attached to the parent residue
    if(carbon != NULL)
    {
        Atom* oxygen = parent_target_atom; ///The oxygen atom of the parent residue that is attached to the new residue
        if(oxygen != NULL)
        {
            int oxygen_index = 1;
            if(oxygen->GetName().size() > 1 && isdigit(oxygen->GetName().at(1)))
                oxygen_index = ConvertString<int>(ConvertT<char>(oxygen->GetName().at(1)));

            int carbon_index = 1;
            if(carbon->GetName().size() > 1 && isdigit(carbon->GetName().at(1)))
                carbon_index = ConvertString<int>(ConvertT<char>(carbon->GetName().at(1)));

            Atom* atom1 = NULL;
            Atom* atom2 = carbon;
            Atom* atom3 = oxygen;
            Atom* atom4 = NULL;

            AtomVector oxygen_neighbors = oxygen->GetNode()->GetNodeNeighbors();
            for(AtomVector::iterator it = oxygen_neighbors.begin(); it != oxygen_neighbors.end(); it++)
            {
                Atom* neighbor = *it;
                if(neighbor->GetId().compare(carbon->GetId()) != 0)
                {
                    if(neighbor->GetName().at(0) == 'C' &&
                            (neighbor->GetName().size() > 1 && isdigit(neighbor->GetName().at(1)) &&
                             ConvertString<int>(ConvertT<char>(neighbor->GetName().at(1))) == oxygen_index))
                    {
                        atom4 = neighbor;
                        break;
                    }
                }
            }
            AtomVector atom2_neighbors = atom2->GetNode()->GetNodeNeighbors();
            for(AtomVector::iterator it = atom2_neighbors.begin(); it != atom2_neighbors.end(); it++)
            {
                Atom* neighbor = *it;
                if(neighbor->GetId().compare(oxygen->GetId()) != 0)
                {
                    if(neighbor->GetName().at(0) == 'C' &&
                            (neighbor->GetName().size() > 1 && isdigit(neighbor->GetName().at(1)) &&
                             ConvertString<int>(ConvertT<char>(neighbor->GetName().at(1))) == carbon_index - 1))
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

void Assembly::SetPsiTorsion(Residue *residue, Residue *parent_residue, int branch_index, double torsion, bool crystallographic_definition)
{
    Atom* residue_head_atom = residue->GetHeadAtoms().at(0);
    Atom* parent_target_atom = parent_residue->GetTailAtoms().at(branch_index);

    Atom* carbon = residue_head_atom; ///The carbon atom of the new residue that is attached to the parent residue
    if(carbon != NULL)
    {
        Atom* oxygen = parent_target_atom; ///The oxygen atom of the parent residue that is attached to the new residue
        if(oxygen != NULL)
        {
            int oxygen_index = 1;
            if(oxygen->GetName().size() > 1 && isdigit(oxygen->GetName().at(1)))
                oxygen_index = ConvertString<int>(ConvertT<char>(oxygen->GetName().at(1)));

            Atom* atom1 = carbon;
            Atom* atom2 = oxygen;
            Atom* atom3 = NULL;
            Atom* atom4 = NULL;

            AtomVector oxygen_neighbors = oxygen->GetNode()->GetNodeNeighbors();
            for(AtomVector::iterator it = oxygen_neighbors.begin(); it != oxygen_neighbors.end(); it++)
            {
                Atom* neighbor = *it;
                if(neighbor->GetId().compare(carbon->GetId()) != 0)
                {
                    if(neighbor->GetName().at(0) == 'C' &&
                            (neighbor->GetName().size() > 1 && isdigit(neighbor->GetName().at(1)) &&
                             ConvertString<int>(ConvertT<char>(neighbor->GetName().at(1))) == oxygen_index))
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
                    Atom* neighbor = *it;
                    if(neighbor->GetId().compare(oxygen->GetId()) != 0)
                    {
                        if(crystallographic_definition)
                        {
                            if(neighbor->GetName().at(0) == 'C' &&
                                    (neighbor->GetName().size() > 1 && isdigit(neighbor->GetName().at(1)) &&
                                     ConvertString<int>(ConvertT<char>(neighbor->GetName().at(1))) == oxygen_index + 1))
                            {
                                atom4 = neighbor;
                                break;
                            }
                        }
                        else
                        {
                            if(neighbor->GetName().at(0) == 'H' &&
                                    (neighbor->GetName().size() > 1 && isdigit(neighbor->GetName().at(1)) &&
                                     ConvertString<int>(ConvertT<char>(neighbor->GetName().at(1))) == oxygen_index))
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

void Assembly::SetOmegaTorsion(Residue *residue, Residue *parent_residue, int branch_index, double torsion, int type)
{
    Atom* residue_head_atom = residue->GetHeadAtoms().at(0);
    Atom* parent_target_atom = parent_residue->GetTailAtoms().at(branch_index);

    Atom* carbon = residue_head_atom; ///The carbon atom of the new residue that is attached to the parent residue
    if(carbon != NULL)
    {
        Atom* oxygen = parent_target_atom; ///The oxygen atom of the parent residue that is attached to the new residue
        if(oxygen != NULL)
        {
            int oxygen_index = 1;
            if(oxygen->GetName().size() > 1 && isdigit(oxygen->GetName().at(1)))
                oxygen_index = ConvertString<int>(ConvertT<char>(oxygen->GetName().at(1)));
            if(oxygen_index == 6 && type == 6)
            {
                Atom* atom1 = NULL;
                Atom* atom2 = NULL;
                Atom* atom3 = NULL;
                Atom* atom4 = oxygen;

                AtomVector oxygen_neighbors = oxygen->GetNode()->GetNodeNeighbors();
                for(AtomVector::iterator it = oxygen_neighbors.begin(); it != oxygen_neighbors.end(); it++)
                {
                    Atom* neighbor = *it;
                    if(neighbor->GetId().compare(carbon->GetId()) != 0)
                    {
                        if(neighbor->GetName().at(0) == 'C' &&
                                (neighbor->GetName().size() > 1 && isdigit(neighbor->GetName().at(1)) &&
                                 ConvertString<int>(ConvertT<char>(neighbor->GetName().at(1))) == oxygen_index))
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
                        Atom* neighbor = *it;
                        if(neighbor->GetId().compare(oxygen->GetId()) != 0)
                        {
                            if(neighbor->GetName().at(0) == 'C' &&
                                    (neighbor->GetName().size() > 1 && isdigit(neighbor->GetName().at(1)) &&
                                     ConvertString<int>(ConvertT<char>(neighbor->GetName().at(1))) == oxygen_index - 1))
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
                        Atom* neighbor = *it;
                        if(neighbor->GetId().compare(atom3->GetId()) != 0)
                        {
                            if(neighbor->GetName().at(0) == 'O' &&
                                    (neighbor->GetName().size() > 1 && isdigit(neighbor->GetName().at(1)) &&
                                     ConvertString<int>(ConvertT<char>(neighbor->GetName().at(1))) == oxygen_index - 1))
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
                Atom* atom1 = NULL;
                Atom* atom2 = NULL;
                Atom* atom3 = NULL;
                Atom* atom4 = NULL;
                Atom* atom5 = NULL;

                AtomVector oxygen_neighbors = oxygen->GetNode()->GetNodeNeighbors();
                for(AtomVector::iterator it = oxygen_neighbors.begin(); it != oxygen_neighbors.end(); it++)
                {
                    Atom* neighbor = *it;
                    if(neighbor->GetId().compare(carbon->GetId()) != 0)
                    {
                        if(neighbor->GetName().at(0) == 'C' &&
                                (neighbor->GetName().size() > 1 && isdigit(neighbor->GetName().at(1)) &&
                                 ConvertString<int>(ConvertT<char>(neighbor->GetName().at(1))) == oxygen_index))
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
                        Atom* neighbor = *it;
                        if(neighbor->GetId().compare(oxygen->GetId()) != 0)
                        {
                            if(neighbor->GetName().at(0) == 'C' &&
                                    (neighbor->GetName().size() > 1 && isdigit(neighbor->GetName().at(1)) &&
                                     ConvertString<int>(ConvertT<char>(neighbor->GetName().at(1))) == oxygen_index - 1))
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
                        Atom* neighbor = *it;
                        if(neighbor->GetId().compare(atom5->GetId()) != 0)
                        {
                            if(neighbor->GetName().at(0) == 'H' &&
                                    (neighbor->GetName().size() > 1 && isdigit(neighbor->GetName().at(1)) &&
                                     ConvertString<int>(ConvertT<char>(neighbor->GetName().at(1))) == oxygen_index - 1))
                            {
                                atom1 = neighbor;
                                break;
                            }
                        }
                    }
                    for(AtomVector::iterator it = atom2_neighbors.begin(); it != atom2_neighbors.end(); it++)
                    {
                        Atom* neighbor = *it;
                        if(neighbor->GetId().compare(atom5->GetId()) != 0)
                        {
                            if(neighbor->GetName().at(0) == 'C' &&
                                    (neighbor->GetName().size() > 1 && isdigit(neighbor->GetName().at(1)) &&
                                     ConvertString<int>(ConvertT<char>(neighbor->GetName().at(1))) == oxygen_index - 2))
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
                        Atom* neighbor = *it;
                        if(neighbor->GetId().compare(atom2->GetId()) != 0)
                        {
                            if(neighbor->GetName().at(0) == 'H' &&
                                    (neighbor->GetName().size() > 1 && isdigit(neighbor->GetName().at(1)) &&
                                     ConvertString<int>(ConvertT<char>(neighbor->GetName().at(1))) == oxygen_index - 2))
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
                Atom* atom1 = NULL;
                Atom* atom2 = NULL;
                Atom* atom3 = NULL;
                Atom* atom4 = NULL;

                AtomVector oxygen_neighbors = oxygen->GetNode()->GetNodeNeighbors();
                for(AtomVector::iterator it = oxygen_neighbors.begin(); it != oxygen_neighbors.end(); it++)
                {
                    Atom* neighbor = *it;
                    if(neighbor->GetId().compare(carbon->GetId()) != 0)
                    {
                        if(neighbor->GetName().at(0) == 'C' &&
                                (neighbor->GetName().size() > 1 && isdigit(neighbor->GetName().at(1)) &&
                                 ConvertString<int>(ConvertT<char>(neighbor->GetName().at(1))) == oxygen_index))
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
                        Atom* neighbor = *it;
                        if(neighbor->GetId().compare(oxygen->GetId()) != 0)
                        {
                            if(neighbor->GetName().at(0) == 'H' &&
                                    (neighbor->GetName().size() > 1 && isdigit(neighbor->GetName().at(1)) &&
                                     ConvertString<int>(ConvertT<char>(neighbor->GetName().at(1))) == oxygen_index))
                            {
                                atom1 = neighbor;
                                break;
                            }
                        }
                    }
                    for(AtomVector::iterator it = atom2_neighbors.begin(); it != atom2_neighbors.end(); it++)
                    {
                        Atom* neighbor = *it;
                        if(neighbor->GetId().compare(oxygen->GetId()) != 0)
                        {
                            if(neighbor->GetName().at(0) == 'C' &&
                                    (neighbor->GetName().size() > 1 && isdigit(neighbor->GetName().at(1)) &&
                                     ConvertString<int>(ConvertT<char>(neighbor->GetName().at(1))) == oxygen_index - 1))
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
                        Atom* neighbor = *it;
                        if(neighbor->GetId().compare(atom2->GetId()) != 0)
                        {
                            if(neighbor->GetName().at(0) == 'H' &&
                                    (neighbor->GetName().size() > 1 && isdigit(neighbor->GetName().at(1)) &&
                                     ConvertString<int>(ConvertT<char>(neighbor->GetName().at(1))) == oxygen_index - 1))
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
                Atom* atom1 = NULL;
                Atom* atom2 = NULL;
                Atom* atom3 = NULL;
                Atom* atom4 = NULL;

                AtomVector oxygen_neighbors = oxygen->GetNode()->GetNodeNeighbors();
                for(AtomVector::iterator it = oxygen_neighbors.begin(); it != oxygen_neighbors.end(); it++)
                {
                    Atom* neighbor = *it;
                    if(neighbor->GetId().compare(carbon->GetId()) != 0)
                    {
                        if(neighbor->GetName().at(0) == 'C' &&
                                (neighbor->GetName().size() > 1 && isdigit(neighbor->GetName().at(1)) &&
                                 ConvertString<int>(ConvertT<char>(neighbor->GetName().at(1))) == oxygen_index))
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
                        Atom* neighbor = *it;
                        if(neighbor->GetId().compare(oxygen->GetId()) != 0)
                        {
                            if(neighbor->GetName().at(0) == 'C' &&
                                    (neighbor->GetName().size() > 1 && isdigit(neighbor->GetName().at(1)) &&
                                     ConvertString<int>(ConvertT<char>(neighbor->GetName().at(1))) == oxygen_index - 1))
                            {
                                atom4 = neighbor;
                                break;
                            }
                        }
                    }
                    for(AtomVector::iterator it = atom3_neighbors.begin(); it != atom3_neighbors.end(); it++)
                    {
                        Atom* neighbor = *it;
                        if(neighbor->GetId().compare(oxygen->GetId()) != 0)
                        {
                            if(neighbor->GetName().at(0) == 'C' &&
                                    (neighbor->GetName().size() > 1 && isdigit(neighbor->GetName().at(1)) &&
                                     ConvertString<int>(ConvertT<char>(neighbor->GetName().at(1))) == oxygen_index + 1))
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
                        Atom* neighbor = *it;
                        if(neighbor->GetId().compare(atom3->GetId()) != 0)
                        {
                            if(neighbor->GetName().at(0) == 'O' &&
                                    (neighbor->GetName().size() > 1 && isdigit(neighbor->GetName().at(1)) &&
                                     ConvertString<int>(ConvertT<char>(neighbor->GetName().at(1))) == oxygen_index + 1))
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

void Assembly::SetOmegaDerivativeTorsion(Residue *residue, Residue *parent_residue, int branch_index, double torsion)
{
    string residue_name = residue->GetName();
    Atom* parent_target_atom = parent_residue->GetTailAtoms().at(branch_index);

    Atom* oxygen = parent_target_atom; ///The oxygen atom of the parent residue that is attached to the new residue
    if(oxygen != NULL)
    {
        int oxygen_index = 1;
        if(oxygen->GetName().size() > 1 && isdigit(oxygen->GetName().at(1)))
            oxygen_index = ConvertString<int>(ConvertT<char>(oxygen->GetName().at(1)));
        if(residue_name.compare("SO3") == 0 && oxygen_index == 6)
        {
            Atom* atom1 = NULL;
            Atom* atom2 = NULL;
            Atom* atom3 = NULL;
            Atom* atom4 = oxygen;

            AtomVector oxygen_neighbors = oxygen->GetNode()->GetNodeNeighbors();
            for(AtomVector::iterator it = oxygen_neighbors.begin(); it != oxygen_neighbors.end(); it++)
            {
                Atom* neighbor = *it;
                if(neighbor->GetName().at(0) == 'C' &&
                        (neighbor->GetName().size() > 1 && isdigit(neighbor->GetName().at(1)) &&
                         ConvertString<int>(ConvertT<char>(neighbor->GetName().at(1))) == oxygen_index))
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
                    Atom* neighbor = *it;
                    if(neighbor->GetId().compare(oxygen->GetId()) != 0)
                    {
                        if(neighbor->GetName().at(0) == 'C' &&
                                (neighbor->GetName().size() > 1 && isdigit(neighbor->GetName().at(1)) &&
                                 ConvertString<int>(ConvertT<char>(neighbor->GetName().at(1))) == oxygen_index - 1))
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
                    Atom* neighbor = *it;
                    if(neighbor->GetId().compare(atom3->GetId()) != 0)
                    {
                        if(neighbor->GetName().at(0) == 'O' &&
                                (neighbor->GetName().size() > 1 && isdigit(neighbor->GetName().at(1)) &&
                                 ConvertString<int>(ConvertT<char>(neighbor->GetName().at(1))) == oxygen_index - 1))
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

void Assembly::SetDihedral(Atom *atom1, Atom *atom2, Atom *atom3, Atom *atom4, double torsion)
{
    double current_dihedral = 0.0;
    Coordinate* a1 = atom1->GetCoordinates().at(model_index_);
    Coordinate* a2 = atom2->GetCoordinates().at(model_index_);
    Coordinate* a3 = atom3->GetCoordinates().at(model_index_);
    Coordinate* a4 = atom4->GetCoordinates().at(model_index_);

    Coordinate* b1 = new Coordinate(*a2);
    b1->operator -(*a1);
    Coordinate* b2 = new Coordinate(*a3);
    b2->operator -(*a2);
    Coordinate* b3 = new Coordinate(*a4);
    b3->operator -(*a3);
    Coordinate* b4 = new Coordinate(*b2);
    b4->operator *(-1);

    Coordinate* b2xb3 = new Coordinate(*b2);
    b2xb3->CrossProduct(*b3);

    Coordinate* b1_m_b2n = new Coordinate(*b1);
    b1_m_b2n->operator *(b2->length());

    Coordinate* b1xb2 = new Coordinate(*b1);
    b1xb2->CrossProduct(*b2);

    current_dihedral = atan2(b1_m_b2n->DotProduct(*b2xb3), b1xb2->DotProduct(*b2xb3));

    double** torsion_matrix = GenerateRotationMatrix(b4, a2, current_dihedral - ConvertDegree2Radian(torsion));

    AtomVector atomsToRotate = AtomVector();
    atomsToRotate.push_back(atom2);
    atom3->FindConnectedAtoms(atomsToRotate);

    for(AtomVector::iterator it = atomsToRotate.begin(); it != atomsToRotate.end(); it++)
    {
        Coordinate* atom_coordinate = (*it)->GetCoordinates().at(model_index_);
        Coordinate* result = new Coordinate();
        result->SetX(torsion_matrix[0][0] * atom_coordinate->GetX() + torsion_matrix[0][1] * atom_coordinate->GetY() +
                torsion_matrix[0][2] * atom_coordinate->GetZ() + torsion_matrix[0][3]);
        result->SetY(torsion_matrix[1][0] * atom_coordinate->GetX() + torsion_matrix[1][1] * atom_coordinate->GetY() +
                torsion_matrix[1][2] * atom_coordinate->GetZ() + torsion_matrix[1][3]);
        result->SetZ(torsion_matrix[2][0] * atom_coordinate->GetX() + torsion_matrix[2][1] * atom_coordinate->GetY() +
                torsion_matrix[2][2] * atom_coordinate->GetZ() + torsion_matrix[2][3]);

        (*it)->GetCoordinates().at(model_index_)->SetX(result->GetX());
        (*it)->GetCoordinates().at(model_index_)->SetY(result->GetY());
        (*it)->GetCoordinates().at(model_index_)->SetZ(result->GetZ());
    }
}

void Assembly::SetAngle(Atom* atom1, Atom* atom2, Atom* atom3, double angle)
{
    double current_angle = 0.0;
    Coordinate* a1 = atom1->GetCoordinates().at(model_index_);
    Coordinate* a2 = atom2->GetCoordinates().at(model_index_);
    Coordinate* a3 = atom3->GetCoordinates().at(model_index_);

    Coordinate* b1 = new Coordinate(*a1);
    b1->operator -(*a2);
    Coordinate* b2 = new Coordinate(*a3);
    b2->operator -(*a2);

    current_angle = acos((b1->DotProduct(*b2)) / (b1->length() * b2->length() + DIST_EPSILON));
    double rotation_angle = ConvertDegree2Radian(angle) - current_angle;

    Coordinate* direction = new Coordinate(*b1);
    direction->CrossProduct(*b2);
    direction->Normalize();
    double** rotation_matrix = GenerateRotationMatrix(direction, a2, rotation_angle);

    AtomVector atomsToRotate = AtomVector();
    atomsToRotate.push_back(atom2);
    atom3->FindConnectedAtoms(atomsToRotate);

    for(AtomVector::iterator it = atomsToRotate.begin() + 1; it != atomsToRotate.end(); it++)
    {
        Coordinate* atom_coordinate = (*it)->GetCoordinates().at(model_index_);
        Coordinate* result = new Coordinate();
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

double Assembly::CalculateBondAngleByCoordinates(Coordinate* atom1_crd, Coordinate* atom2_crd, Coordinate* atom3_crd)
{
    Coordinate* b1 = new Coordinate(*atom1_crd);
    b1->operator -(*atom2_crd);
    Coordinate* b2 = new Coordinate(*atom3_crd);
    b2->operator -(*atom2_crd);

    return acos(b1->DotProduct((*b2)) / b1->length() / b2->length());

}

double Assembly::CalculateBondAngleByAtoms(Atom *atom1, Atom *atom2, Atom *atom3)
{
    Coordinate* a1 = atom1->GetCoordinates().at(model_index_);
    Coordinate* a2 = atom2->GetCoordinates().at(model_index_);
    Coordinate* a3 = atom3->GetCoordinates().at(model_index_);

    Coordinate* b1 = new Coordinate(*a1);
    b1->operator -(*a2);
    Coordinate* b2 = new Coordinate(*a3);
    b2->operator -(*a2);

    return acos(b1->DotProduct((*b2)) / b1->length() / b2->length());
}

double Assembly::CalculateTorsionAngleByCoordinates(Coordinate* atom1_crd, Coordinate* atom2_crd, Coordinate* atom3_crd, Coordinate* atom4_crd)
{
    double current_dihedral = 0.0;

    Coordinate* b1 = new Coordinate(*atom2_crd);
    b1->operator -(*atom1_crd);
    Coordinate* b2 = new Coordinate(*atom3_crd);
    b2->operator -(*atom2_crd);
    Coordinate* b3 = new Coordinate(*atom4_crd);
    b3->operator -(*atom3_crd);
    Coordinate* b4 = new Coordinate(*b2);
    b4->operator *(-1);

    Coordinate* b2xb3 = new Coordinate(*b2);
    b2xb3->CrossProduct(*b3);

    Coordinate* b1_m_b2n = new Coordinate(*b1);
    b1_m_b2n->operator *(b2->length());

    Coordinate* b1xb2 = new Coordinate(*b1);
    b1xb2->CrossProduct(*b2);

    current_dihedral = atan2(b1_m_b2n->DotProduct(*b2xb3), b1xb2->DotProduct(*b2xb3));
    return current_dihedral;
}

double Assembly::CalculateTorsionAngleByAtoms(Atom *atom1, Atom *atom2, Atom *atom3, Atom *atom4)
{
    double current_dihedral = 0.0;
    Coordinate* a1 = atom1->GetCoordinates().at(model_index_);
    Coordinate* a2 = atom2->GetCoordinates().at(model_index_);
    Coordinate* a3 = atom3->GetCoordinates().at(model_index_);
    Coordinate* a4 = atom4->GetCoordinates().at(model_index_);

    Coordinate* b1 = new Coordinate(*a2);
    b1->operator -(*a1);
    Coordinate* b2 = new Coordinate(*a3);
    b2->operator -(*a2);
    Coordinate* b3 = new Coordinate(*a4);
    b3->operator -(*a3);
    Coordinate* b4 = new Coordinate(*b2);
    b4->operator *(-1);

    Coordinate* b2xb3 = new Coordinate(*b2);
    b2xb3->CrossProduct(*b3);

    Coordinate* b1_m_b2n = new Coordinate(*b1);
    b1_m_b2n->operator *(b2->length());

    Coordinate* b1xb2 = new Coordinate(*b1);
    b1xb2->CrossProduct(*b2);

    current_dihedral = atan2(b1_m_b2n->DotProduct(*b2xb3), b1xb2->DotProduct(*b2xb3));
    return current_dihedral;
}

void Assembly::GetCenterOfMass(Coordinate *center_of_mass)
{
    //    center_of_mass = new Coordinate();
    AtomVector all_atoms_of_assembly = this->GetAllAtomsOfAssembly();
    for(AtomVector::iterator it = all_atoms_of_assembly.begin(); it != all_atoms_of_assembly.end(); it++)
    {
        Atom* atom = *it;
        center_of_mass->Translate(atom->MolecularDynamicAtom::GetMass() * atom->GetCoordinates().at(model_index_)->GetX(),
                                  atom->MolecularDynamicAtom::GetMass() * atom->GetCoordinates().at(model_index_)->GetY(),
                                  atom->MolecularDynamicAtom::GetMass() * atom->GetCoordinates().at(model_index_)->GetZ());
    }
    double total_mass = this->GetTotalMass();
    center_of_mass->operator /(Coordinate(center_of_mass->GetX() / total_mass,
                                          center_of_mass->GetY() / total_mass,
                                          center_of_mass->GetZ() / total_mass));
}

void Assembly::GetCenterOfGeometry(Coordinate *center_of_geometry)
{
    double sumX = 0.0;
    double sumY = 0.0;
    double sumZ = 0.0;
    CoordinateVector all_coords = this->GetAllCoordinates();
    for(CoordinateVector::iterator it = all_coords.begin(); it != all_coords.end(); it++)
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

void Assembly::GetBoundary(Coordinate* lower_left_back_corner, Coordinate* upper_right_front_corner)
{
    //    lower_left_back_corner = new Coordinate(-INFINITY, -INFINITY, -INFINITY);
    //    upper_right_front_corner = new Coordinate(INFINITY, INFINITY, INFINITY);
    lower_left_back_corner->SetX(INFINITY);
    lower_left_back_corner->SetY(INFINITY);
    lower_left_back_corner->SetZ(INFINITY);
    upper_right_front_corner->SetX(-INFINITY);
    upper_right_front_corner->SetY(-INFINITY);
    upper_right_front_corner->SetZ(-INFINITY);
    AtomVector all_atoms_of_assembly = this->GetAllAtomsOfAssembly();
    for(AtomVector::iterator it = all_atoms_of_assembly.begin(); it != all_atoms_of_assembly.end(); it++)
    {
        Atom* atom = *it;
        if(atom->MolecularDynamicAtom::GetRadius() == dNotSet)
        {
            //            gmml::log(__LINE__, __FILE__,  gmml::ERR, "There is no information of the atom type/radius/charge of the atoms in the given library/parameter file");
            //            cout << "There is no information of the atom type/radius/charge of the atoms in the given library/parameter file" << endl;
            atom->MolecularDynamicAtom::SetRadius(DEFAULT_RADIUS);
            //            stringstream ss;
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

string Assembly::CalculateRSOrientations(Atom *prev_atom, Atom *target, Atom *next_atom)
{
    string orientation = "";
    ///Calculating the plane based on the two ring neighbors of the current atom
    Coordinate prev_atom_coord = Coordinate(*prev_atom->GetCoordinates().at(model_index_));
    Coordinate current_atom_coord = Coordinate(*target->GetCoordinates().at(model_index_));
    Coordinate next_atom_coord = Coordinate(*next_atom->GetCoordinates().at(model_index_));
    prev_atom_coord.operator -(current_atom_coord) ;
    next_atom_coord.operator -(current_atom_coord) ;
    Plane plane = Plane();
    plane.SetV1(prev_atom_coord);
    plane.SetV2(next_atom_coord);
    Coordinate normal_v = plane.GetUnitNormalVector();

    AtomNode* node = target->GetNode();
    AtomVector neighbors = node->GetNodeNeighbors();
    for(AtomVector::iterator it1 = neighbors.begin(); it1 != neighbors.end(); it1++)
    {
        Atom* neighbor = (*it1);
        if(neighbor->GetId().at(0) == 'O')
        {
            Coordinate side_atom_coord = Coordinate(*neighbor->GetCoordinates().at(model_index_));
            side_atom_coord.operator -(current_atom_coord);
            side_atom_coord.Normalize();
            Coordinate normal_v_x_side_atom = normal_v;
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
    AtomVector assemblyAAtoms = this->GetAllAtomsOfAssembly();
    AtomVector assemblyBAtoms = assemblyB->GetAllAtomsOfAssembly();

    double rA = 0.0, rB = 0.0, distance = 0.0, totalOverlap = 0.0;

    for(AtomVector::iterator atomA = assemblyAAtoms.begin(); atomA != assemblyAAtoms.end(); atomA++)
    {
        for(AtomVector::iterator atomB = assemblyBAtoms.begin(); atomB != assemblyBAtoms.end(); atomB++)
        {
            distance = (*atomA)->GetDistanceToAtom( (*atomB));
            if ( ( distance < 3.6 ) && ( distance > 0.0 ) ) //Close enough to overlap, but not the same atom
            {
                // element info not set, so I look at first letter of atom name.
                if ((*atomA)->GetName().at(0) == 'C') rA = 1.70; // Rowland and Taylor modification
                if ((*atomA)->GetName().at(0) == 'O') rA = 1.52;
                if ((*atomA)->GetName().at(0) == 'N') rA = 1.55;
                if ((*atomA)->GetName().at(0) == 'S') rA = 1.80;
                if ((*atomA)->GetName().at(0) == 'P') rA = 1.80;
                if ((*atomA)->GetName().at(0) == 'H') rA = 1.09;

                if ((*atomB)->GetName().at(0) == 'C') rB = 1.70;
                if ((*atomB)->GetName().at(0) == 'O') rB = 1.52;
                if ((*atomB)->GetName().at(0) == 'N') rB = 1.55;
                if ((*atomB)->GetName().at(0) == 'S') rB = 1.80;
                if ((*atomA)->GetName().at(0) == 'P') rA = 1.80;
                if ((*atomA)->GetName().at(0) == 'H') rA = 1.09;

       //         std::cout << "Distance=" << distance << " rA=" << rA << " rB=" << rB << std::endl;
                if (rA + rB > distance + 0.6){ // 0.6 overlap is deemed acceptable. (Copying chimera:)
                    // Eqn 1, Rychkov and Petukhov, J. Comput. Chem., 2006, Joint Neighbours...
                    // Each atom against each atom, so overlap can be "double" counted. See paper.
                    totalOverlap += ( 2 * (PI_RADIAN) * rA* ( rA - distance / 2 - ( ( (rA*rA) - (rB*rB) ) / (2 * distance) ) ) );
                    //std::cout << "Overlap=" << totalOverlap << std::endl;
                }
            }
        }
    }
    return (totalOverlap / CARBON_SURFACE_AREA); //Normalise to area of a buried carbon
}

