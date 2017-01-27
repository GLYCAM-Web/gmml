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
bool Assembly::CheckCondensedSequenceSanity(string sequence, CondensedSequence::CondensedSequenceAmberPrepResidueTree& prep_residues)
{
    try
    {
        CondensedSequence* condensed_sequence = new CondensedSequence(sequence);
        prep_residues = condensed_sequence->GetCondensedSequenceAmberPrepResidueTree();
        for(CondensedSequence::CondensedSequenceAmberPrepResidueTree::iterator it = prep_residues.begin(); it != prep_residues.end(); it++)
        {
            CondensedSequenceAmberPrepResidue* amber_prep_residue = *it;
            string amber_prep_residue_name = amber_prep_residue->GetName();
            if(amber_prep_residue_name.compare("UNK") == 0)
            {
                cout << "The input sequence (" << sequence << ") is not valid" << endl;
                return false;
            }
        }
    }
    catch(exception ex)
    {
        cout << "The input sequence (" << sequence << ") is not valid" << endl;
        return false;
    }

    cout << "The input sequence (" << sequence << ") is valid" << endl;
    return true;
}

void Assembly::BuildAssemblyFromCondensedSequence(string sequence, string prep_file, string parameter_file, bool structure)
{

    try
    {
        CondensedSequence* condensed_sequence = new CondensedSequence(sequence);
        CondensedSequence::CondensedSequenceAmberPrepResidueTree amber_prep_residues = condensed_sequence->GetCondensedSequenceAmberPrepResidueTree();
        PrepFile* prep = new PrepFile(prep_file);
        PrepFile::ResidueMap prep_residue_map = prep->GetResidues();
        ParameterFile* parameter = NULL;
        ParameterFile::AtomTypeMap atom_type_map = ParameterFile::AtomTypeMap();
        if(parameter_file.compare("") != 0)
        {
            parameter = new ParameterFile(parameter_file);
            atom_type_map = parameter->GetAtomTypes();
        }
        ResidueVector parent_residues = ResidueVector();
        ResidueVector branch_residues = ResidueVector();
        vector<bool> derivatives = vector<bool>();
        int sequence_number = 0;
        int serial_number = 0;
        stringstream ss;
        for(CondensedSequence::CondensedSequenceAmberPrepResidueTree::iterator it = amber_prep_residues.begin(); it != amber_prep_residues.end(); ++it)
        {
            CondensedSequenceAmberPrepResidue* amber_prep_residue = *it;
            string amber_prep_residue_name = amber_prep_residue->GetName();
            string amber_prep_residue_parent_oxygen = amber_prep_residue->GetParentOxygen();

            if(prep_residue_map.find(amber_prep_residue_name) != prep_residue_map.end())
            {
                PrepFileResidue* prep_residue = prep_residue_map[amber_prep_residue_name];

                // Build residue from prep residue
                sequence_number++;
                CoordinateVector cartesian_coordinate_list = CoordinateVector();

                Residue* assembly_residue = new Residue();
                assembly_residue->SetAssembly(this);
                string prep_residue_name = prep_residue->GetName();
                assembly_residue->SetName(prep_residue_name);
                stringstream id;
                id << prep_residue_name << "_" << gmml::BLANK_SPACE << "_" << sequence_number << "_" << gmml::BLANK_SPACE << "_"
                   << gmml::BLANK_SPACE << "_" << id_;
                assembly_residue->SetId(id.str());
                if(distance(amber_prep_residues.begin(), it) == (int)amber_prep_residues.size()-1)
                    ss << prep_residue_name;
                else
                    ss << prep_residue_name << "-";

                PrepFileResidue::PrepFileAtomVector prep_atoms = prep_residue->GetAtoms();
                for(PrepFileResidue::PrepFileAtomVector::iterator it1 = prep_atoms.begin(); it1 != prep_atoms.end(); it1++)
                {
                    PrepFileAtom* prep_atom = (*it1);
                    string atom_name = prep_atom->GetName();
                    serial_number++;
                    Atom* assembly_atom = new Atom();
                    assembly_atom->SetResidue(assembly_residue);
                    assembly_atom->SetName(atom_name);
                    stringstream atom_id;
                    atom_id << atom_name << "_" << serial_number << "_" << id.str();
                    assembly_atom->SetId(atom_id.str());

                    assembly_atom->MolecularDynamicAtom::SetAtomType(prep_atom->GetType());
                    assembly_atom->MolecularDynamicAtom::SetCharge(prep_atom->GetCharge());
                    if(parameter != NULL)
                    {
                        if(atom_type_map.find(assembly_atom->GetAtomType()) != atom_type_map.end())
                        {
                            ParameterFileAtom* parameter_atom = atom_type_map[assembly_atom->GetAtomType()];
                            assembly_atom->MolecularDynamicAtom::SetMass(parameter_atom->GetMass());
                            assembly_atom->MolecularDynamicAtom::SetRadius(parameter_atom->GetRadius());
                        }
                        else
                        {
                            assembly_atom->MolecularDynamicAtom::SetMass(dNotSet);
                            assembly_atom->MolecularDynamicAtom::SetRadius(dNotSet);
                        }
                    }
                    else
                    {
                        assembly_atom->MolecularDynamicAtom::SetMass(dNotSet);
                        assembly_atom->MolecularDynamicAtom::SetRadius(dNotSet);
                    }

                    if(atom_name.at(0) == 'P' && prep_residue_name.compare("PO3") == 0)
                        assembly_residue->AddHeadAtom(assembly_atom);
                    else if(atom_name.at(0) == 'S' && prep_residue_name.compare("SO3") == 0)
                    {
                        assembly_residue->AddHeadAtom(assembly_atom);
                    }
                    else if(atom_name.at(0) == 'C' && prep_residue_name.compare("MEX") == 0)
                    {
                        assembly_residue->AddHeadAtom(assembly_atom);
                    }
                    else if(atom_name.find("C1") != string::npos && prep_residue_name.compare("ACX") == 0)
                    {
                        assembly_residue->AddHeadAtom(assembly_atom);
                    }

                    if(prep_residue->GetCoordinateType() == PrepFileSpace::kINT)
                    {
                        vector<Coordinate*> coordinate_list = vector<Coordinate*>();
                        int index = distance(prep_atoms.begin(), it1);
                        if(index == 0)
                        {
                        }
                        if(index == 1)
                        {
                            int parent_index = prep_atom->GetBondIndex() - 1;
                            Coordinate* parent_coordinate = cartesian_coordinate_list.at(parent_index);
                            coordinate_list.push_back(parent_coordinate);
                        }
                        if(index == 2)
                        {
                            int grandparent_index = prep_atom->GetAngleIndex() - 1;
                            int parent_index = prep_atom->GetBondIndex() - 1;
                            Coordinate* grandparent_coordinate = cartesian_coordinate_list.at(grandparent_index);
                            Coordinate* parent_coordinate = cartesian_coordinate_list.at(parent_index);
                            coordinate_list.push_back(grandparent_coordinate);
                            coordinate_list.push_back(parent_coordinate);
                        }
                        if(index > 2)
                        {
                            int great_grandparent_index = prep_atom->GetDihedralIndex() - 1;
                            int grandparent_index = prep_atom->GetAngleIndex() - 1;
                            int parent_index = prep_atom->GetBondIndex() - 1;

                            Coordinate* great_grandparent_coordinate = cartesian_coordinate_list.at(great_grandparent_index);
                            Coordinate* grandparent_coordinate = cartesian_coordinate_list.at(grandparent_index);
                            Coordinate* parent_coordinate = cartesian_coordinate_list.at(parent_index);
                            coordinate_list.push_back(great_grandparent_coordinate);
                            coordinate_list.push_back(grandparent_coordinate);
                            coordinate_list.push_back(parent_coordinate);
                        }
                        Coordinate* coordinate = new Coordinate();
                        coordinate = gmml::ConvertInternalCoordinate2CartesianCoordinate(coordinate_list, prep_atom->GetBondLength(),
                                                                                         prep_atom->GetAngle(), prep_atom->GetDihedral());
                        cartesian_coordinate_list.push_back(coordinate);

                        assembly_atom->AddCoordinate(coordinate);
                    }
                    else if(prep_residue->GetCoordinateType() == PrepFileSpace::kXYZ)
                    {
                        assembly_atom->AddCoordinate(new Coordinate(prep_atom->GetBondLength(), prep_atom->GetAngle(), prep_atom->GetDihedral()));
                    }
                    if(assembly_atom->GetAtomType().compare("DU") != 0)
                        assembly_residue->AddAtom(assembly_atom);
                    if(atom_name.compare(amber_prep_residue->GetAnomericCarbon()) == 0)
                        assembly_residue->AddHeadAtom(assembly_atom);
                }

                if(structure)
                {
                    Assembly* temp_assembly = new Assembly();
                    ResidueVector temp_assembly_residues = ResidueVector();
                    temp_assembly_residues.push_back(assembly_residue);
                    temp_assembly->SetResidues(temp_assembly_residues);
                    temp_assembly->SetSourceFile(prep_file);
                    temp_assembly->BuildStructureByPrepFileInformation();
                }
                residues_.push_back(assembly_residue);
                if(amber_prep_residue->GetParentId() != -1)
                {
                    Residue* parent_residue = residues_.at(amber_prep_residue->GetParentId());
                    AtomVector parent_residue_atoms = parent_residue->GetAtoms();
                    for(AtomVector::iterator it3 = parent_residue_atoms.begin(); it3 != parent_residue_atoms.end(); it3++)
                    {
                        Atom* parent_atom = *it3;
                        if(parent_atom->GetName().compare(amber_prep_residue->GetParentOxygen()) == 0)
                            parent_residue->AddTailAtom(parent_atom);
                    }
                    parent_residues.push_back(parent_residue);
                    branch_residues.push_back(assembly_residue);
                    if(amber_prep_residue->GetIsDerivative())
                        derivatives.push_back(true);
                    else
                        derivatives.push_back(false);
                }
            }
            else
            {
                cout << "Residue " << amber_prep_residue_name << " has not been found in the database" << endl;
            }
        }

        name_ = ss.str();
        this->SetSourceFile(prep_file);

        if(structure)
        {
            map<Residue*, int> parent_branch_map = map<Residue*, int>();
            int linkage_index = -1;
            for(ResidueVector::iterator it = parent_residues.begin(); it != parent_residues.end(); it++)
            {
                Residue* parent_residue = (*it);
                int parent_index = distance(parent_residues.begin(), it);
                if(parent_branch_map.find(parent_residue) == parent_branch_map.end())
                    parent_branch_map[parent_residue] = 0;
                else
                    parent_branch_map[parent_residue]++;
                Residue* assembly_residue = branch_residues.at(parent_index);

                int branch_index = parent_branch_map[parent_residue];

                if(!derivatives.at(parent_index))
                    this->AttachResidues(assembly_residue, parent_residue, branch_index, parameter_file);
                else //Attach derivative
                {
                    this->AttachResidues(assembly_residue, parent_residue, branch_index, parameter_file);
                    this->RemoveHydrogenAtAttachedPosition(parent_residue, branch_index);
                    this->AdjustCharge(assembly_residue, parent_residue, branch_index);
                    this->SetDerivativeAngle(assembly_residue, parent_residue, branch_index);
                }
                if(linkage_index >= 0)
                {
                    // 2-8 default rotamer
                    stringstream linkage_name;
                    linkage_name << assembly_residue->GetName() << assembly_residue->GetHeadAtoms().at(0)->GetName().at(1) << "-" <<
                                    parent_residue->GetTailAtoms().at(branch_index)->GetName().at(1) << parent_residue->GetName();
                    if(linkage_name.str().find("0SA2-8") != string::npos ||
                                                linkage_name.str().find("0SB2-8") != string::npos)
                    {
                        this->SetPhiTorsion(assembly_residue, parent_residue, branch_index, EXTERNAL28LINKAGEROTAMERS[0][0]);
                        this->SetPsiTorsion(assembly_residue, parent_residue, branch_index, EXTERNAL28LINKAGEROTAMERS[0][1], false);
                        this->SetOmegaTorsion(assembly_residue, parent_residue, branch_index, EXTERNAL28LINKAGEROTAMERS[0][2], 7);
                        this->SetOmegaTorsion(assembly_residue, parent_residue, branch_index, EXTERNAL28LINKAGEROTAMERS[0][3], 8);
                        this->SetOmegaTorsion(assembly_residue, parent_residue, branch_index, EXTERNAL28LINKAGEROTAMERS[0][4], 9);
                    }
                    if(linkage_name.str().find("0GL2-8") != string::npos)
                    {
                        this->SetPhiTorsion(assembly_residue, parent_residue, branch_index, INTERNAL28LINKAGEROTAMERS[0][0]);
                        this->SetPsiTorsion(assembly_residue, parent_residue, branch_index, INTERNAL28LINKAGEROTAMERS[0][1], false);
                        this->SetOmegaTorsion(assembly_residue, parent_residue, branch_index, INTERNAL28LINKAGEROTAMERS[0][2], 7);
                        this->SetOmegaTorsion(assembly_residue, parent_residue, branch_index, INTERNAL28LINKAGEROTAMERS[0][3], 8);
                        this->SetOmegaTorsion(assembly_residue, parent_residue, branch_index, INTERNAL28LINKAGEROTAMERS[0][4], 9);
                    }
                }
                linkage_index++;
            }
        }
    }
    catch(exception ex)
    {
        cout << "Building assembly from " << sequence << " failed." << endl;
    }
}

Assembly::AssemblyVector Assembly::BuildAllRotamersFromCondensedSequence(string sequence, string prep_file, string parameter_file,
                                                                         CondensedSequence::CondensedSequenceRotamersAndGlycosidicAnglesInfo rotamers_glycosidic_angles_info,
                                                                         CondensedSequence::IndexNameMap& names)
{

    try
    {
        CondensedSequence* condensed_sequence = new CondensedSequence(sequence);
        AssemblyVector structures = AssemblyVector(condensed_sequence->CountAllPossible28LinkagesRotamers(rotamers_glycosidic_angles_info) *
                                                   condensed_sequence->CountAllPossibleSelectedRotamers(rotamers_glycosidic_angles_info));
        CondensedSequence::IndexLinkageConfigurationMap structure_map = condensed_sequence->CreateIndexLinkageConfigurationMap(
                    rotamers_glycosidic_angles_info, names);
        CondensedSequence::CondensedSequenceAmberPrepResidueTree amber_prep_residues = condensed_sequence->GetCondensedSequenceAmberPrepResidueTree();
        PrepFile* prep = new PrepFile(prep_file);
        PrepFile::ResidueMap prep_residue_map = prep->GetResidues();
        ParameterFile* parameter = NULL;
        ParameterFile::AtomTypeMap atom_type_map = ParameterFile::AtomTypeMap();
        if(parameter_file.compare("") != 0)
        {
            parameter = new ParameterFile(parameter_file);
            atom_type_map = parameter->GetAtomTypes();
        }
        for(unsigned int i = 0; i < structures.size(); i++)
        {
            ResidueVector parent_residues = ResidueVector();
            ResidueVector branch_residues = ResidueVector();
            vector<bool> derivatives = vector<bool>();
            structures.at(i) = new Assembly();
            int sequence_number = 0;
            int serial_number = 0;
            stringstream ss;
            for(CondensedSequence::CondensedSequenceAmberPrepResidueTree::iterator it = amber_prep_residues.begin(); it != amber_prep_residues.end(); ++it)
            {
                CondensedSequenceAmberPrepResidue* amber_prep_residue = *it;
                string amber_prep_residue_name = amber_prep_residue->GetName();
                string amber_prep_residue_parent_oxygen = amber_prep_residue->GetParentOxygen();

                if(prep_residue_map.find(amber_prep_residue_name) != prep_residue_map.end())
                {
                    PrepFileResidue* prep_residue = prep_residue_map[amber_prep_residue_name];

                    // Build residue from prep residue
                    sequence_number++;
                    CoordinateVector cartesian_coordinate_list = CoordinateVector();

                    Residue* assembly_residue = new Residue();
                    assembly_residue->SetAssembly(structures.at(i));
                    string prep_residue_name = prep_residue->GetName();
                    assembly_residue->SetName(prep_residue_name);
                    stringstream id;
                    id << prep_residue_name << "_" << gmml::BLANK_SPACE << "_" << sequence_number << "_" << gmml::BLANK_SPACE << "_"
                       << gmml::BLANK_SPACE << "_" << id_;
                    assembly_residue->SetId(id.str());
                    if(distance(amber_prep_residues.begin(), it) == (int)amber_prep_residues.size()-1)
                        ss << prep_residue_name;
                    else
                        ss << prep_residue_name << "-";

                    PrepFileResidue::PrepFileAtomVector prep_atoms = prep_residue->GetAtoms();
                    for(PrepFileResidue::PrepFileAtomVector::iterator it1 = prep_atoms.begin(); it1 != prep_atoms.end(); it1++)
                    {
                        serial_number++;
                        Atom* assembly_atom = new Atom();
                        PrepFileAtom* prep_atom = (*it1);
                        assembly_atom->SetResidue(assembly_residue);
                        string atom_name = prep_atom->GetName();
                        assembly_atom->SetName(atom_name);
                        stringstream atom_id;
                        atom_id << atom_name << "_" << serial_number << "_" << id.str();
                        assembly_atom->SetId(atom_id.str());

                        assembly_atom->MolecularDynamicAtom::SetAtomType(prep_atom->GetType());
                        assembly_atom->MolecularDynamicAtom::SetCharge(prep_atom->GetCharge());
                        if(parameter != NULL)
                        {
                            if(atom_type_map.find(assembly_atom->GetAtomType()) != atom_type_map.end())
                            {
                                ParameterFileAtom* parameter_atom = atom_type_map[assembly_atom->GetAtomType()];
                                assembly_atom->MolecularDynamicAtom::SetMass(parameter_atom->GetMass());
                                assembly_atom->MolecularDynamicAtom::SetRadius(parameter_atom->GetRadius());
                            }
                            else
                            {
                                assembly_atom->MolecularDynamicAtom::SetMass(dNotSet);
                                assembly_atom->MolecularDynamicAtom::SetRadius(dNotSet);
                            }
                        }
                        else
                        {
                            assembly_atom->MolecularDynamicAtom::SetMass(dNotSet);
                            assembly_atom->MolecularDynamicAtom::SetRadius(dNotSet);
                        }

                        if(atom_name.at(0) == 'P' && prep_residue_name.compare("PO3") == 0)
                            assembly_residue->AddHeadAtom(assembly_atom);
                        else if(atom_name.at(0) == 'S' && prep_residue_name.compare("SO3") == 0)
                        {
                            assembly_residue->AddHeadAtom(assembly_atom);
                        }
                        else if(atom_name.at(0) == 'C' && prep_residue_name.compare("MEX") == 0)
                        {
                            assembly_residue->AddHeadAtom(assembly_atom);
                        }
                        else if(atom_name.find("C1") != string::npos && prep_residue_name.compare("ACX") == 0)
                        {
                            assembly_residue->AddHeadAtom(assembly_atom);
                        }

                        if(prep_residue->GetCoordinateType() == PrepFileSpace::kINT)
                        {
                            vector<Coordinate*> coordinate_list = vector<Coordinate*>();
                            int index = distance(prep_atoms.begin(), it1);
                            if(index == 0)
                            {
                            }
                            if(index == 1)
                            {
                                int parent_index = prep_atom->GetBondIndex() - 1;
                                Coordinate* parent_coordinate = cartesian_coordinate_list.at(parent_index);
                                coordinate_list.push_back(parent_coordinate);
                            }
                            if(index == 2)
                            {
                                int grandparent_index = prep_atom->GetAngleIndex() - 1;
                                int parent_index = prep_atom->GetBondIndex() - 1;
                                Coordinate* grandparent_coordinate = cartesian_coordinate_list.at(grandparent_index);
                                Coordinate* parent_coordinate = cartesian_coordinate_list.at(parent_index);
                                coordinate_list.push_back(grandparent_coordinate);
                                coordinate_list.push_back(parent_coordinate);
                            }
                            if(index > 2)
                            {
                                int great_grandparent_index = prep_atom->GetDihedralIndex() - 1;
                                int grandparent_index = prep_atom->GetAngleIndex() - 1;
                                int parent_index = prep_atom->GetBondIndex() - 1;

                                Coordinate* great_grandparent_coordinate = cartesian_coordinate_list.at(great_grandparent_index);
                                Coordinate* grandparent_coordinate = cartesian_coordinate_list.at(grandparent_index);
                                Coordinate* parent_coordinate = cartesian_coordinate_list.at(parent_index);
                                coordinate_list.push_back(great_grandparent_coordinate);
                                coordinate_list.push_back(grandparent_coordinate);
                                coordinate_list.push_back(parent_coordinate);
                            }
                            Coordinate* coordinate = new Coordinate();
                            coordinate = gmml::ConvertInternalCoordinate2CartesianCoordinate(coordinate_list, prep_atom->GetBondLength(),
                                                                                             prep_atom->GetAngle(), prep_atom->GetDihedral());
                            cartesian_coordinate_list.push_back(coordinate);

                            assembly_atom->AddCoordinate(coordinate);
                        }
                        else if(prep_residue->GetCoordinateType() == PrepFileSpace::kXYZ)
                        {
                            assembly_atom->AddCoordinate(new Coordinate(prep_atom->GetBondLength(), prep_atom->GetAngle(), prep_atom->GetDihedral()));
                        }
                        if(assembly_atom->GetAtomType().compare("DU") != 0)
                            assembly_residue->AddAtom(assembly_atom);
                        if(atom_name.compare(amber_prep_residue->GetAnomericCarbon()) == 0)
                            assembly_residue->AddHeadAtom(assembly_atom);
                    }

                    if(true)
                    {
                        Assembly* temp_assembly = new Assembly();
                        ResidueVector temp_assembly_residues = ResidueVector();
                        temp_assembly_residues.push_back(assembly_residue);
                        temp_assembly->SetResidues(temp_assembly_residues);
                        temp_assembly->SetSourceFile(prep_file);
                        temp_assembly->BuildStructureByPrepFileInformation();
                    }
                    structures.at(i)->residues_.push_back(assembly_residue);
                    if(amber_prep_residue->GetParentId() != -1)
                    {
                        Residue* parent_residue = structures.at(i)->residues_.at(amber_prep_residue->GetParentId());
                        AtomVector parent_residue_atoms = parent_residue->GetAtoms();
                        for(AtomVector::iterator it3 = parent_residue_atoms.begin(); it3 != parent_residue_atoms.end(); it3++)
                        {
                            Atom* parent_atom = *it3;
                            if(parent_atom->GetName().compare(amber_prep_residue->GetParentOxygen()) == 0)
                                parent_residue->AddTailAtom(parent_atom);
                        }
                        parent_residues.push_back(parent_residue);
                        branch_residues.push_back(assembly_residue);
                        if(amber_prep_residue->GetIsDerivative())
                            derivatives.push_back(true);
                        else
                            derivatives.push_back(false);
                    }
                }
                else
                {
                    cout << "Residue " << amber_prep_residue_name << " has not been found in the database" << endl;
                }
            }

            structures.at(i)->name_ = ss.str();
            structures.at(i)->SetSourceFile(prep_file);

            if(true)
            {
                map<Residue*, int> parent_branch_map = map<Residue*, int>();
                int linkage_index = -1;
                for(ResidueVector::iterator it = parent_residues.begin(); it != parent_residues.end(); it++)
                {
                    Residue* parent_residue = (*it);
                    int parent_index = distance(parent_residues.begin(), it);
                    if(parent_branch_map.find(parent_residue) == parent_branch_map.end())
                        parent_branch_map[parent_residue] = 0;
                    else
                        parent_branch_map[parent_residue]++;
                    Residue* assembly_residue = branch_residues.at(parent_index);

                    int branch_index = parent_branch_map[parent_residue];


                    if(!derivatives.at(parent_index))
                        structures.at(i)->AttachResidues(assembly_residue, parent_residue, branch_index, parameter_file);
                    else
                    {
                        //Attach derivative
                        structures.at(i)->AttachResidues(assembly_residue, parent_residue, branch_index, parameter_file);
                        structures.at(i)->RemoveHydrogenAtAttachedPosition(parent_residue, branch_index);
                        structures.at(i)->AdjustCharge(assembly_residue, parent_residue, branch_index);
                        structures.at(i)->SetDerivativeAngle(assembly_residue, parent_residue, branch_index);
                    }
                    if(linkage_index >= 0)
                    {
                        if(!derivatives.at(parent_index))
                        {
                            /*cout << rotamers_glycosidic_angles_info.at(linkage_index).first << " "
                             << assembly_residue->GetName() << "(" << assembly_residue->GetHeadAtoms().at(0)->GetName() << "-"
                             << parent_residue->GetTailAtoms().at(branch_index)->GetName() << ")" << parent_residue->GetName() << endl;*/
                            vector<double> phi_psi_omega = structure_map[i].at(linkage_index);
                            if(phi_psi_omega.at(0) != dNotSet)
                            {
                                double phi = phi_psi_omega.at(0);
                                structures.at(i)->SetPhiTorsion(assembly_residue, parent_residue, branch_index, phi);// Set phi angle of assembly_residue-parent_residue to phi
                            }
                            if(phi_psi_omega.at(1) != dNotSet)
                            {
                                double psi = phi_psi_omega.at(1);
                                structures.at(i)->SetPsiTorsion(assembly_residue, parent_residue, branch_index, psi);// Set psi angle of assembly_residue-parent_residue to psi
                            }
                            if(phi_psi_omega.at(2) != dNotSet)
                            {
                                double omega = phi_psi_omega.at(2);
                                structures.at(i)->SetOmegaTorsion(assembly_residue, parent_residue, branch_index, omega);// Set omega angle of assembly_residue-parent_residue to omega
                            }
                            if(phi_psi_omega.size() > 3) //2-8 linkages
                            {
                                structures.at(i)->SetPhiTorsion(assembly_residue, parent_residue, branch_index, phi_psi_omega.at(3));
                                structures.at(i)->SetPsiTorsion(assembly_residue, parent_residue, branch_index, phi_psi_omega.at(4), false);
                                structures.at(i)->SetOmegaTorsion(assembly_residue, parent_residue, branch_index, phi_psi_omega.at(5), 7);
                                structures.at(i)->SetOmegaTorsion(assembly_residue, parent_residue, branch_index, phi_psi_omega.at(6), 8);
                                structures.at(i)->SetOmegaTorsion(assembly_residue, parent_residue, branch_index, phi_psi_omega.at(7), 9);
                            }
                        }
                        else
                        {
                            vector<double> phi_psi_omega = structure_map[i].at(linkage_index);
                            if(phi_psi_omega.at(2) != dNotSet)
                            {
                                double omega = phi_psi_omega.at(2);
                                structures.at(i)->SetOmegaDerivativeTorsion(assembly_residue, parent_residue, branch_index, omega);
                            }
                        }
                    }
                    linkage_index++;
                }
            }
        }
        return structures;
    }
    catch(exception ex)
    {
        cout << "Building assembly from " << sequence << " failed." << endl;
    }
}

void Assembly::BuildAssemblyFromPdbFile(string pdb_file_path, vector<string> amino_lib_files, vector<string> glycam_lib_files,
                                        vector<string> other_lib_files, vector<string> prep_files, string parameter_file)
{
    cout << "Building assembly from pdb file ..." << endl;
    gmml::log(__LINE__, __FILE__, gmml::INF, "Building assembly from pdb file ...");
    try
    {
        this->ClearAssembly();
        PdbFile* pdb_file = new PdbFile(pdb_file_path);
        ParameterFile* parameter = NULL;
        ParameterFile::AtomTypeMap atom_type_map = ParameterFile::AtomTypeMap();
        if(parameter_file.compare("") != 0)
        {
            parameter = new ParameterFile(parameter_file);
            atom_type_map = parameter->GetAtomTypes();
        }

        LibraryFile::ResidueMap lib_residues = LibraryFile::ResidueMap();
        PrepFile::ResidueMap prep_residues = PrepFile::ResidueMap();
        vector<string> lib_files = vector<string>();
        if(!amino_lib_files.empty())
            for(vector<string>::iterator it = amino_lib_files.begin(); it != amino_lib_files.end(); it++)
                lib_files.push_back(*it);
        if(!glycam_lib_files.empty())
            for(vector<string>::iterator it = glycam_lib_files.begin(); it != glycam_lib_files.end(); it++)
                lib_files.push_back(*it);
        if(!other_lib_files.empty())
            for(vector<string>::iterator it = other_lib_files.begin(); it != other_lib_files.end(); it++)
                lib_files.push_back(*it);

        if(lib_files.size() != 0)
            lib_residues = GetAllResiduesFromMultipleLibFilesMap(lib_files);

        if(prep_files.size() != 0)
            prep_residues = GetAllResiduesFromMultiplePrepFilesMap(prep_files);

        vector<string> key_order = vector<string>();
        PdbFile::PdbResidueAtomsMap residue_atoms_map = pdb_file->GetAllAtomsInOrder(key_order);
        for(vector<string>::iterator it = key_order.begin(); it != key_order.end(); it++)
        {
            string residue_key = *it;
            PdbFile::PdbAtomVector* atoms = residue_atoms_map[residue_key];
            Residue* residue = new Residue();
            residue->SetAssembly(this);

            for(PdbFile::PdbAtomVector::iterator it1 = atoms->begin(); it1 != atoms->end(); it1++)
            {
                PdbAtom* atom = (*it1);
                string residue_name = atom->GetAtomResidueName();
                char chain_id = atom->GetAtomChainId();
                int sequence_number = atom->GetAtomResidueSequenceNumber();
                char insertion_code = atom->GetAtomInsertionCode();
                char alternate_location = atom->GetAtomAlternateLocation();
                stringstream ss;
                ss << residue_name << "_" << chain_id << "_" << sequence_number << "_" << insertion_code << "_"
                   << alternate_location << "_" << id_;
                string key = ss.str();
                residue->SetId(key);

                Atom* new_atom = new Atom();
                residue->SetName(residue_name);
                string atom_name = atom->GetAtomName();
                new_atom->SetName(atom_name);
                if(!lib_residues.empty() || !prep_residues.empty())
                {
                    if(lib_residues.find(residue_name) != lib_residues.end())
                    {
                        LibraryFileResidue* lib_residue = lib_residues[residue_name];
                        LibraryFileAtom* lib_atom = lib_residue->GetLibraryAtomByAtomName(atom_name);
                        if(lib_atom != NULL)
                        {
                            new_atom->MolecularDynamicAtom::SetAtomType(lib_atom->GetType());
                            new_atom->MolecularDynamicAtom::SetCharge(lib_atom->GetCharge());

                            if(parameter != NULL)
                            {
                                if(atom_type_map.find(new_atom->MolecularDynamicAtom::GetAtomType()) != atom_type_map.end())
                                {
                                    ParameterFileAtom* parameter_atom = atom_type_map[new_atom->GetAtomType()];
                                    new_atom->MolecularDynamicAtom::SetMass(parameter_atom->GetMass());
                                    new_atom->MolecularDynamicAtom::SetRadius(parameter_atom->GetRadius());
                                }
                                else
                                {
                                    new_atom->MolecularDynamicAtom::SetMass(dNotSet);
                                    new_atom->MolecularDynamicAtom::SetRadius(dNotSet);
                                }
                            }
                            else
                            {
                                new_atom->MolecularDynamicAtom::SetMass(dNotSet);
                                new_atom->MolecularDynamicAtom::SetRadius(dNotSet);
                            }
                        }
                        else
                        {
                            new_atom->MolecularDynamicAtom::SetAtomType("UNK");
                            new_atom->MolecularDynamicAtom::SetCharge(dNotSet);
                            new_atom->MolecularDynamicAtom::SetMass(dNotSet);
                            new_atom->MolecularDynamicAtom::SetRadius(dNotSet);
                        }
                    }
                    else if(prep_residues.find(residue_name) != prep_residues.end())
                    {
                        PrepFileResidue* prep_residue = prep_residues[residue_name];
                        PrepFileAtom* prep_atom = prep_residue->GetPrepAtomByName(atom_name);
                        if(prep_atom != NULL)
                        {
                            new_atom->MolecularDynamicAtom::SetAtomType(prep_atom->GetType());
                            new_atom->MolecularDynamicAtom::SetCharge(prep_atom->GetCharge());

                            if(parameter != NULL)
                            {
                                if(atom_type_map.find(new_atom->MolecularDynamicAtom::GetAtomType()) != atom_type_map.end())
                                {
                                    ParameterFileAtom* parameter_atom = atom_type_map[new_atom->GetAtomType()];
                                    new_atom->MolecularDynamicAtom::SetMass(parameter_atom->GetMass());
                                    new_atom->MolecularDynamicAtom::SetRadius(parameter_atom->GetRadius());
                                }
                                else
                                {
                                    new_atom->MolecularDynamicAtom::SetMass(dNotSet);
                                    new_atom->MolecularDynamicAtom::SetRadius(dNotSet);
                                }
                            }
                            else
                            {
                                new_atom->MolecularDynamicAtom::SetMass(dNotSet);
                                new_atom->MolecularDynamicAtom::SetRadius(dNotSet);
                            }
                        }
                        else
                        {
                            new_atom->MolecularDynamicAtom::SetAtomType("UNK");
                            new_atom->MolecularDynamicAtom::SetCharge(dNotSet);
                            new_atom->MolecularDynamicAtom::SetMass(dNotSet);
                            new_atom->MolecularDynamicAtom::SetRadius(dNotSet);
                        }
                    }
                }
                new_atom->SetResidue(residue);
                stringstream atom_key;
                atom_key << atom_name << "_" << atom->GetAtomSerialNumber() << "_" << key;
                new_atom->SetId(atom_key.str());
                PdbModelCard* models = pdb_file->GetModels();
                PdbModelCard::PdbModelMap model_maps = models->GetModels();
                if(model_maps.size() == 1)
                {
                    new_atom->AddCoordinate(new GeometryTopology::Coordinate(atom->GetAtomOrthogonalCoordinate()));
                    vector<string> card_index = gmml::Split(atom->GetAtomCardIndexInResidueSet(), "_");
                    if(card_index.at(0).compare("ATOM") == 0)
                    {
                        new_atom->SetDescription("Atom;");
                    }
                    else if(card_index.at(0).compare("HETATOM") == 0)
                    {
                        new_atom->SetDescription("Het;");
                    }
                }
                else
                {
                    for(PdbModelCard::PdbModelMap::iterator it2 = model_maps.begin(); it2 != model_maps.end(); it2++)
                    {
                        PdbModel* model = (*it2).second;
                        PdbModelResidueSet* residue_set = model->GetModelResidueSet();
                        PdbModelResidueSet::AtomCardVector atom_cards = residue_set->GetAtoms();
                        vector<string> card_index = gmml::Split(atom->GetAtomCardIndexInResidueSet(), "_");
                        if(card_index.at(0).compare("ATOM") == 0)
                        {
                            PdbAtomCard* atom_card = atom_cards.at(gmml::ConvertString<int>(card_index.at(1)));
                            PdbAtomCard::PdbAtomOrderVector atom_vector = atom_card->GetOrderedAtoms();
                            for(PdbAtomCard::PdbAtomOrderVector::iterator it3 = atom_vector.begin(); it3 != atom_vector.end(); it3++)
                            {
                                PdbAtom* matching_atom = *it3;
                                string matching_residue_name = matching_atom->GetAtomResidueName();
                                char matching_chain_id = matching_atom->GetAtomChainId();
                                int matching_sequence_number = matching_atom->GetAtomResidueSequenceNumber();
                                char matching_insertion_code = matching_atom->GetAtomInsertionCode();
                                char matching_alternate_location = matching_atom->GetAtomAlternateLocation();
                                stringstream sss;
                                sss << matching_residue_name << "_" << matching_chain_id << "_" << matching_sequence_number << "_"
                                    << matching_insertion_code << "_" << matching_alternate_location << "_" << id_;
                                string matching_key = sss.str();

                                if(key.compare(matching_key) == 0)
                                {
                                    GeometryTopology::Coordinate* coordinate = new GeometryTopology::Coordinate(matching_atom->GetAtomOrthogonalCoordinate());
                                    new_atom->AddCoordinate(coordinate);
                                    new_atom->SetDescription("Atom;");
                                }
                            }
                        }
                        else if(card_index.at(0).compare("HETATOM") == 0)
                        {
                            PdbModelResidueSet::HeterogenAtomCardVector heterogen_atom_cards = residue_set->GetHeterogenAtoms();
                            PdbHeterogenAtomCard* heterogen_atom_card = heterogen_atom_cards.at(gmml::ConvertString<int>(card_index.at(1)));
                            PdbHeterogenAtomCard::PdbHeterogenAtomOrderVector heterogen_atom_vector = heterogen_atom_card->GetOrderedHeterogenAtoms();
                            for(PdbHeterogenAtomCard::PdbHeterogenAtomOrderVector::iterator it3 = heterogen_atom_vector.begin(); it3 != heterogen_atom_vector.end(); it3++)
                            {
                                PdbAtom* matching_heterogen_atom = *it3;
                                string matching_heterogen_residue_name = matching_heterogen_atom->GetAtomResidueName();
                                char matching_heterogen_chain_id = matching_heterogen_atom->GetAtomChainId();
                                int matching_heterogen_sequence_number = matching_heterogen_atom->GetAtomResidueSequenceNumber();
                                char matching_heterogen_insertion_code = matching_heterogen_atom->GetAtomInsertionCode();
                                char matching_heterogen_alternate_location = matching_heterogen_atom->GetAtomAlternateLocation();
                                stringstream ssss;
                                ssss << matching_heterogen_residue_name << "_" << matching_heterogen_chain_id << "_" << matching_heterogen_sequence_number << "_"
                                     << matching_heterogen_insertion_code << "_" << matching_heterogen_alternate_location << "_" << id_;
                                string matching_heterogen_key = ssss.str();

                                if(key.compare(matching_heterogen_key) == 0)
                                {
                                    GeometryTopology::Coordinate* coordinate = new GeometryTopology::Coordinate(matching_heterogen_atom->GetAtomOrthogonalCoordinate());
                                    new_atom->AddCoordinate(coordinate);
                                    new_atom->SetDescription("Het;");
                                }
                            }
                        }
                    }
                }
                residue->AddAtom(new_atom);
            }
            this->AddResidue(residue);
        }
    }
    catch(PdbFileSpace::PdbFileProcessingException &ex)
    {}
}

void Assembly::BuildAssemblyFromPdbFile(PdbFile *pdb_file, vector<string> amino_lib_files, vector<string> glycam_lib_files,
                                        vector<string> other_lib_files, vector<string> prep_files, string parameter_file)
{
    cout << "Building assembly from pdb file ..." << endl;
    gmml::log(__LINE__, __FILE__, gmml::INF, "Building assembly from pdb file ...");
    try
    {
        this->ClearAssembly();
        ParameterFile* parameter = NULL;
        ParameterFile::AtomTypeMap atom_type_map = ParameterFile::AtomTypeMap();
        if(parameter_file.compare("") != 0)
        {
            parameter = new ParameterFile(parameter_file);
            atom_type_map = parameter->GetAtomTypes();
        }

        LibraryFile::ResidueMap lib_residues = LibraryFile::ResidueMap();
        PrepFile::ResidueMap prep_residues = PrepFile::ResidueMap();
        vector<string> lib_files = vector<string>();
        if(!amino_lib_files.empty())
            for(vector<string>::iterator it = amino_lib_files.begin(); it != amino_lib_files.end(); it++)
                lib_files.push_back(*it);
        if(!glycam_lib_files.empty())
            for(vector<string>::iterator it = glycam_lib_files.begin(); it != glycam_lib_files.end(); it++)
                lib_files.push_back(*it);
        if(!other_lib_files.empty())
            for(vector<string>::iterator it = other_lib_files.begin(); it != other_lib_files.end(); it++)
                lib_files.push_back(*it);

        if(lib_files.size() != 0)
            lib_residues = GetAllResiduesFromMultipleLibFilesMap(lib_files);

        if(prep_files.size() != 0)
            prep_residues = GetAllResiduesFromMultiplePrepFilesMap(prep_files);

        vector<string> key_order = vector<string>();
        PdbFile::PdbResidueAtomsMap residue_atoms_map = pdb_file->GetAllAtomsInOrder(key_order);
        for(vector<string>::iterator it = key_order.begin(); it != key_order.end(); it++)
        {
            string residue_key = *it;
            PdbFile::PdbAtomVector* atoms = residue_atoms_map[residue_key];
            Residue* residue = new Residue();
            residue->SetAssembly(this);

            for(PdbFile::PdbAtomVector::iterator it1 = atoms->begin(); it1 != atoms->end(); it1++)
            {
                PdbAtom* atom = (*it1);
                string residue_name = atom->GetAtomResidueName();
                char chain_id = atom->GetAtomChainId();
                int sequence_number = atom->GetAtomResidueSequenceNumber();
                char insertion_code = atom->GetAtomInsertionCode();
                char alternate_location = atom->GetAtomAlternateLocation();
                stringstream ss;
                ss << residue_name << "_" << chain_id << "_" << sequence_number << "_" << insertion_code << "_"
                   << alternate_location << "_" << id_;
                string key = ss.str();
                residue->SetId(key);

                Atom* new_atom = new Atom();
                residue->SetName(residue_name);
                string atom_name = atom->GetAtomName();
                new_atom->SetName(atom_name);
                if(!lib_residues.empty() || !prep_residues.empty())
                {
                    if(lib_residues.find(residue_name) != lib_residues.end())
                    {
                        LibraryFileResidue* lib_residue = lib_residues[residue_name];
                        LibraryFileAtom* lib_atom = lib_residue->GetLibraryAtomByAtomName(atom_name);
                        if(lib_atom != NULL)
                        {
                            new_atom->MolecularDynamicAtom::SetAtomType(lib_atom->GetType());
                            new_atom->MolecularDynamicAtom::SetCharge(lib_atom->GetCharge());

                            if(parameter != NULL)
                            {
                                if(atom_type_map.find(new_atom->MolecularDynamicAtom::GetAtomType()) != atom_type_map.end())
                                {
                                    ParameterFileAtom* parameter_atom = atom_type_map[new_atom->GetAtomType()];
                                    new_atom->MolecularDynamicAtom::SetMass(parameter_atom->GetMass());
                                    new_atom->MolecularDynamicAtom::SetRadius(parameter_atom->GetRadius());
                                }
                                else
                                {
                                    new_atom->MolecularDynamicAtom::SetMass(dNotSet);
                                    new_atom->MolecularDynamicAtom::SetRadius(dNotSet);
                                }
                            }
                            else
                            {
                                new_atom->MolecularDynamicAtom::SetMass(dNotSet);
                                new_atom->MolecularDynamicAtom::SetRadius(dNotSet);
                            }
                        }
                        else
                        {
                            new_atom->MolecularDynamicAtom::SetAtomType("UNK");
                            new_atom->MolecularDynamicAtom::SetCharge(dNotSet);
                            new_atom->MolecularDynamicAtom::SetMass(dNotSet);
                            new_atom->MolecularDynamicAtom::SetRadius(dNotSet);
                        }
                    }
                    else if(prep_residues.find(residue_name) != prep_residues.end())
                    {
                        PrepFileResidue* prep_residue = prep_residues[residue_name];
                        PrepFileAtom* prep_atom = prep_residue->GetPrepAtomByName(atom_name);
                        if(prep_atom != NULL)
                        {
                            new_atom->MolecularDynamicAtom::SetAtomType(prep_atom->GetType());
                            new_atom->MolecularDynamicAtom::SetCharge(prep_atom->GetCharge());

                            if(parameter != NULL)
                            {
                                if(atom_type_map.find(new_atom->MolecularDynamicAtom::GetAtomType()) != atom_type_map.end())
                                {
                                    ParameterFileAtom* parameter_atom = atom_type_map[new_atom->GetAtomType()];
                                    new_atom->MolecularDynamicAtom::SetMass(parameter_atom->GetMass());
                                    new_atom->MolecularDynamicAtom::SetRadius(parameter_atom->GetRadius());
                                }
                                else
                                {
                                    new_atom->MolecularDynamicAtom::SetMass(dNotSet);
                                    new_atom->MolecularDynamicAtom::SetRadius(dNotSet);
                                }
                            }
                            else
                            {
                                new_atom->MolecularDynamicAtom::SetMass(dNotSet);
                                new_atom->MolecularDynamicAtom::SetRadius(dNotSet);
                            }
                        }
                        else
                        {
                            new_atom->MolecularDynamicAtom::SetAtomType("UNK");
                            new_atom->MolecularDynamicAtom::SetCharge(dNotSet);
                            new_atom->MolecularDynamicAtom::SetMass(dNotSet);
                            new_atom->MolecularDynamicAtom::SetRadius(dNotSet);
                        }
                    }
                }
                new_atom->SetResidue(residue);
                stringstream atom_key;
                atom_key << atom_name << "_" << atom->GetAtomSerialNumber() << "_" << key;
                new_atom->SetId(atom_key.str());
                PdbModelCard* models = pdb_file->GetModels();
                PdbModelCard::PdbModelMap model_maps = models->GetModels();
                if(model_maps.size() == 1)
                {
                    new_atom->AddCoordinate(new GeometryTopology::Coordinate(atom->GetAtomOrthogonalCoordinate()));
                    vector<string> card_index = gmml::Split(atom->GetAtomCardIndexInResidueSet(), "_");
                    if(card_index.at(0).compare("ATOM") == 0)
                    {
                        new_atom->SetDescription("Atom;");
                    }
                    else if(card_index.at(0).compare("HETATOM") == 0)
                    {
                        new_atom->SetDescription("Het;");
                    }
                }
                else
                {
                    for(PdbModelCard::PdbModelMap::iterator it2 = model_maps.begin(); it2 != model_maps.end(); it2++)
                    {
                        PdbModel* model = (*it2).second;
                        PdbModelResidueSet* residue_set = model->GetModelResidueSet();
                        PdbModelResidueSet::AtomCardVector atom_cards = residue_set->GetAtoms();
                        vector<string> card_index = gmml::Split(atom->GetAtomCardIndexInResidueSet(), "_");
                        if(card_index.at(0).compare("ATOM") == 0)
                        {
                            PdbAtomCard* atom_card = atom_cards.at(gmml::ConvertString<int>(card_index.at(1)));
                            PdbAtomCard::PdbAtomOrderVector atom_vector = atom_card->GetOrderedAtoms();
                            for(PdbAtomCard::PdbAtomOrderVector::iterator it3 = atom_vector.begin(); it3 != atom_vector.end(); it3++)
                            {
                                PdbAtom* matching_atom = *it3;
                                string matching_residue_name = matching_atom->GetAtomResidueName();
                                char matching_chain_id = matching_atom->GetAtomChainId();
                                int matching_sequence_number = matching_atom->GetAtomResidueSequenceNumber();
                                char matching_insertion_code = matching_atom->GetAtomInsertionCode();
                                char matching_alternate_location = matching_atom->GetAtomAlternateLocation();
                                stringstream sss;
                                sss << matching_residue_name << "_" << matching_chain_id << "_" << matching_sequence_number << "_"
                                    << matching_insertion_code << "_" << matching_alternate_location << "_" << id_;
                                string matching_key = sss.str();

                                if(key.compare(matching_key) == 0)
                                {
                                    GeometryTopology::Coordinate* coordinate = new GeometryTopology::Coordinate(matching_atom->GetAtomOrthogonalCoordinate());
                                    new_atom->AddCoordinate(coordinate);
                                    new_atom->SetDescription("Atom;");
                                }
                            }
                        }
                        else if(card_index.at(0).compare("HETATOM") == 0)
                        {
                            PdbModelResidueSet::HeterogenAtomCardVector heterogen_atom_cards = residue_set->GetHeterogenAtoms();
                            PdbHeterogenAtomCard* heterogen_atom_card = heterogen_atom_cards.at(gmml::ConvertString<int>(card_index.at(1)));
                            PdbHeterogenAtomCard::PdbHeterogenAtomOrderVector heterogen_atom_vector = heterogen_atom_card->GetOrderedHeterogenAtoms();
                            for(PdbHeterogenAtomCard::PdbHeterogenAtomOrderVector::iterator it3 = heterogen_atom_vector.begin(); it3 != heterogen_atom_vector.end(); it3++)
                            {
                                PdbAtom* matching_heterogen_atom = *it3;
                                string matching_heterogen_residue_name = matching_heterogen_atom->GetAtomResidueName();
                                char matching_heterogen_chain_id = matching_heterogen_atom->GetAtomChainId();
                                int matching_heterogen_sequence_number = matching_heterogen_atom->GetAtomResidueSequenceNumber();
                                char matching_heterogen_insertion_code = matching_heterogen_atom->GetAtomInsertionCode();
                                char matching_heterogen_alternate_location = matching_heterogen_atom->GetAtomAlternateLocation();
                                stringstream ssss;
                                ssss << matching_heterogen_residue_name << "_" << matching_heterogen_chain_id << "_" << matching_heterogen_sequence_number << "_"
                                     << matching_heterogen_insertion_code << "_" << matching_heterogen_alternate_location << "_" << id_;
                                string matching_heterogen_key = ssss.str();

                                if(key.compare(matching_heterogen_key) == 0)
                                {
                                    GeometryTopology::Coordinate* coordinate = new GeometryTopology::Coordinate(matching_heterogen_atom->GetAtomOrthogonalCoordinate());
                                    new_atom->AddCoordinate(coordinate);
                                    new_atom->SetDescription("Het;");
                                }
                            }
                        }
                    }
                }
                residue->AddAtom(new_atom);
            }
            this->AddResidue(residue);
        }
    }
    catch(PdbFileSpace::PdbFileProcessingException &ex)
    {}
}

void Assembly::BuildAssemblyFromPdbqtFile(string pdbqt_file_path, string parameter_file)
{
    cout << "Building assembly from pdbqt file ..." << endl;
    gmml::log(__LINE__, __FILE__, gmml::INF, "Building assembly from pdbqt file ...");
    try
    {
        this->ClearAssembly();
        PdbqtFile* pdbqt_file = new PdbqtFile(pdbqt_file_path);
        ParameterFile* parameter = NULL;
        ParameterFile::AtomTypeMap atom_type_map = ParameterFile::AtomTypeMap();
        if(parameter_file.compare("") != 0)
        {
            parameter = new ParameterFile(parameter_file);
            atom_type_map = parameter->GetAtomTypes();
        }
        vector<string> key_order = vector<string>();
        PdbqtFile::PdbqtResidueAtomsMap residue_atoms_map = pdbqt_file->GetAllAtomsInOrder(key_order);
        for(vector<string>::iterator it = key_order.begin(); it != key_order.end(); it++)
        {
            string residue_key = *it;
            PdbqtFile::PdbqtAtomVector* atoms = residue_atoms_map[residue_key];
            Residue* residue = new Residue();
            residue->SetAssembly(this);

            for(PdbqtFile::PdbqtAtomVector::iterator it1 = atoms->begin(); it1 != atoms->end(); it1++)
            {
                PdbqtAtom* atom = (*it1);
                string residue_name = atom->GetAtomResidueName();
                char chain_id = atom->GetAtomChainId();
                int sequence_number = atom->GetAtomResidueSequenceNumber();
                char insertion_code = atom->GetAtomInsertionCode();
                char alternate_location = atom->GetAtomAlternateLocation();
                stringstream ss;
                ss << residue_name << "_" << chain_id << "_" << sequence_number << "_" << insertion_code << "_"
                   << alternate_location << "_" << id_;
                string key = ss.str();
                residue->SetId(key);

                Atom* new_atom = new Atom();
                residue->SetName(residue_name);
                string atom_name = atom->GetAtomName();
                new_atom->SetName(atom_name);
                new_atom->MolecularDynamicAtom::SetCharge(atom->GetAtomCharge());
                new_atom->MolecularDynamicAtom::SetAtomType(atom->GetAtomType());
                if(parameter != NULL)
                {
                    if(atom_type_map.find(new_atom->GetAtomType()) != atom_type_map.end())
                    {
                        ParameterFileAtom* parameter_atom = atom_type_map[new_atom->GetAtomType()];
                        new_atom->MolecularDynamicAtom::SetMass(parameter_atom->GetMass());
                        new_atom->MolecularDynamicAtom::SetRadius(parameter_atom->GetRadius());
                    }
                    else
                    {
                        new_atom->MolecularDynamicAtom::SetMass(dNotSet);
                        new_atom->MolecularDynamicAtom::SetRadius(dNotSet);
                    }
                }
                else
                {
                    new_atom->MolecularDynamicAtom::SetMass(dNotSet);
                    new_atom->MolecularDynamicAtom::SetRadius(dNotSet);
                }
                new_atom->SetResidue(residue);
                stringstream atom_key;
                atom_key << atom_name << "_" << atom->GetAtomSerialNumber() << "_" << key;
                new_atom->SetId(atom_key.str());
                PdbqtModelCard* models = pdbqt_file->GetModels();
                PdbqtModelCard::PdbqtModelMap model_maps = models->GetModels();
                if(model_maps.size() == 1)
                {
                    new_atom->AddCoordinate(new GeometryTopology::Coordinate(atom->GetAtomOrthogonalCoordinate()));
                    if(atom->GetType().compare("ATOM") == 0)
                    {
                        new_atom->SetDescription("Atom;");
                    }
                    else if(atom->GetType().compare("HETATOM") == 0)
                    {
                        new_atom->SetDescription("Het;");
                    }
                }
                else
                {
                    for(PdbqtModelCard::PdbqtModelMap::iterator it2 = model_maps.begin(); it2 != model_maps.end(); it2++)
                    {
                        PdbqtModel* model = (*it2).second;
                        PdbqtModelResidueSet* residue_set = model->GetModelResidueSet();
                        PdbqtAtomCard* atom_card = residue_set->GetAtoms();
                        PdbqtAtomCard::PdbqtAtomMap atom_map = atom_card->GetAtoms();
                        PdbqtAtom* matching_atom = atom_map[atom->GetAtomSerialNumber()];
                        string matching_residue_name = matching_atom->GetAtomResidueName();
                        char matching_chain_id = matching_atom->GetAtomChainId();
                        int matching_sequence_number = matching_atom->GetAtomResidueSequenceNumber();
                        char matching_insertion_code = matching_atom->GetAtomInsertionCode();
                        char matching_alternate_location = matching_atom->GetAtomAlternateLocation();
                        stringstream sss;
                        sss << matching_residue_name << "_" << matching_chain_id << "_" << matching_sequence_number << "_"
                            << matching_insertion_code << "_" << matching_alternate_location << "_" << id_;
                        string matching_key = sss.str();

                        if(key.compare(matching_key) == 0)
                        {
                            GeometryTopology::Coordinate* coordinate = new GeometryTopology::Coordinate(matching_atom->GetAtomOrthogonalCoordinate());
                            new_atom->AddCoordinate(coordinate);
                            if(atom->GetType().compare("ATOM") == 0)
                            {
                                new_atom->SetDescription("Atom;");
                            }
                            else if(atom->GetType().compare("HETATOM") == 0)
                            {
                                new_atom->SetDescription("Het;");
                            }
                        }
                    }
                }
                residue->AddAtom(new_atom);
            }
            this->AddResidue(residue);
        }
    }
    catch(PdbFileSpace::PdbFileProcessingException &ex)
    {}
}

void Assembly::BuildAssemblyFromPdbqtFile(PdbqtFile *pdbqt_file, string parameter_file)
{
    cout << "Building assembly from pdbqt file ..." << endl;
    gmml::log(__LINE__, __FILE__, gmml::INF, "Building assembly from pdbqt file ...");
    try
    {
        this->ClearAssembly();
        ParameterFile* parameter = NULL;
        ParameterFile::AtomTypeMap atom_type_map = ParameterFile::AtomTypeMap();
        if(parameter_file.compare("") != 0)
        {
            parameter = new ParameterFile(parameter_file);
            atom_type_map = parameter->GetAtomTypes();
        }
        vector<string> key_order = vector<string>();
        PdbqtFile::PdbqtResidueAtomsMap residue_atoms_map = pdbqt_file->GetAllAtomsInOrder(key_order);
        for(vector<string>::iterator it = key_order.begin(); it != key_order.end(); it++)
        {
            string residue_key = *it;
            PdbqtFile::PdbqtAtomVector* atoms = residue_atoms_map[residue_key];
            Residue* residue = new Residue();
            residue->SetAssembly(this);

            for(PdbqtFile::PdbqtAtomVector::iterator it1 = atoms->begin(); it1 != atoms->end(); it1++)
            {
                PdbqtAtom* atom = (*it1);
                string residue_name = atom->GetAtomResidueName();
                char chain_id = atom->GetAtomChainId();
                int sequence_number = atom->GetAtomResidueSequenceNumber();
                char insertion_code = atom->GetAtomInsertionCode();
                char alternate_location = atom->GetAtomAlternateLocation();
                stringstream ss;
                ss << residue_name << "_" << chain_id << "_" << sequence_number << "_" << insertion_code << "_"
                   << alternate_location << "_" << id_;
                string key = ss.str();
                residue->SetId(key);

                Atom* new_atom = new Atom();
                residue->SetName(residue_name);
                string atom_name = atom->GetAtomName();
                new_atom->SetName(atom_name);
                new_atom->MolecularDynamicAtom::SetCharge(atom->GetAtomCharge());
                new_atom->MolecularDynamicAtom::SetAtomType(atom->GetAtomType());
                if(parameter != NULL)
                {
                    if(atom_type_map.find(new_atom->GetAtomType()) != atom_type_map.end())
                    {
                        ParameterFileAtom* parameter_atom = atom_type_map[new_atom->GetAtomType()];
                        new_atom->MolecularDynamicAtom::SetMass(parameter_atom->GetMass());
                        new_atom->MolecularDynamicAtom::SetRadius(parameter_atom->GetRadius());
                    }
                    else
                    {
                        new_atom->MolecularDynamicAtom::SetMass(dNotSet);
                        new_atom->MolecularDynamicAtom::SetRadius(dNotSet);
                    }
                }
                else
                {
                    new_atom->MolecularDynamicAtom::SetMass(dNotSet);
                    new_atom->MolecularDynamicAtom::SetRadius(dNotSet);
                }
                new_atom->SetResidue(residue);
                stringstream atom_key;
                atom_key << atom_name << "_" << atom->GetAtomSerialNumber() << "_" << key;
                new_atom->SetId(atom_key.str());
                PdbqtModelCard* models = pdbqt_file->GetModels();
                PdbqtModelCard::PdbqtModelMap model_maps = models->GetModels();
                if(model_maps.size() == 1)
                {
                    new_atom->AddCoordinate(new GeometryTopology::Coordinate(atom->GetAtomOrthogonalCoordinate()));
                    if(atom->GetType().compare("ATOM") == 0)
                    {
                        new_atom->SetDescription("Atom;");
                    }
                    else if(atom->GetType().compare("HETATOM") == 0)
                    {
                        new_atom->SetDescription("Het;");
                    }
                }
                else
                {
                    for(PdbqtModelCard::PdbqtModelMap::iterator it2 = model_maps.begin(); it2 != model_maps.end(); it2++)
                    {
                        PdbqtModel* model = (*it2).second;
                        PdbqtModelResidueSet* residue_set = model->GetModelResidueSet();
                        PdbqtAtomCard* atom_card = residue_set->GetAtoms();
                        PdbqtAtomCard::PdbqtAtomMap atom_map = atom_card->GetAtoms();
                        PdbqtAtom* matching_atom = atom_map[atom->GetAtomSerialNumber()];
                        string matching_residue_name = matching_atom->GetAtomResidueName();
                        char matching_chain_id = matching_atom->GetAtomChainId();
                        int matching_sequence_number = matching_atom->GetAtomResidueSequenceNumber();
                        char matching_insertion_code = matching_atom->GetAtomInsertionCode();
                        char matching_alternate_location = matching_atom->GetAtomAlternateLocation();
                        stringstream sss;
                        sss << matching_residue_name << "_" << matching_chain_id << "_" << matching_sequence_number << "_"
                            << matching_insertion_code << "_" << matching_alternate_location << "_" << id_;
                        string matching_key = sss.str();

                        if(key.compare(matching_key) == 0)
                        {
                            GeometryTopology::Coordinate* coordinate = new GeometryTopology::Coordinate(matching_atom->GetAtomOrthogonalCoordinate());
                            new_atom->AddCoordinate(coordinate);
                            if(atom->GetType().compare("ATOM") == 0)
                            {
                                new_atom->SetDescription("Atom;");
                            }
                            else if(atom->GetType().compare("HETATOM") == 0)
                            {
                                new_atom->SetDescription("Het;");
                            }
                        }
                    }
                }
                residue->AddAtom(new_atom);
            }
            this->AddResidue(residue);
        }
    }
    catch(PdbFileSpace::PdbFileProcessingException &ex)
    {}
}

void Assembly::BuildAssemblyFromTopologyFile(string topology_file_path, string parameter_file)
{
    cout << "Building assembly from topology file ..." << endl;
    gmml::log(__LINE__, __FILE__, gmml::INF, "Building assembly from topology file ...");
    this->ClearAssembly();
    TopologyFile* topology_file = new TopologyFile(topology_file_path);
    ParameterFile* parameter = NULL;
    ParameterFile::AtomTypeMap atom_type_map = ParameterFile::AtomTypeMap();
    if(parameter_file.compare("") != 0)
    {
        parameter = new ParameterFile(parameter_file);
        atom_type_map = parameter->GetAtomTypes();
    }
    name_ = topology_file->GetTitle();
    sequence_number_ = 1;
    int serial_number = 0;
    TopologyAssembly::TopologyResidueVector topology_residues = topology_file->GetAssembly()->GetResidues();
    for(TopologyAssembly::TopologyResidueVector::iterator it = topology_residues.begin(); it != topology_residues.end(); it++)
    {
        Residue* assembly_residue = new Residue();
        assembly_residue->SetAssembly(this);
        TopologyResidue* topology_residue = (*it);
        string residue_name = topology_residue->GetResidueName();
        assembly_residue->SetName(residue_name);
        stringstream id;
        id << residue_name << "_" << gmml::BLANK_SPACE << "_" << topology_residue->GetIndex() << "_" << gmml::BLANK_SPACE << "_"
           << gmml::BLANK_SPACE << "_" << id_;
        assembly_residue->SetId(id.str());

        TopologyResidue::TopologyAtomVector topology_atoms = topology_residue->GetAtoms();

        for(TopologyResidue::TopologyAtomVector::iterator it1 = topology_atoms.begin(); it1 != topology_atoms.end(); it1++)
        {
            serial_number++;
            Atom* assembly_atom = new Atom();
            string atom_name = (*it1)->GetAtomName();
            assembly_atom->SetName(atom_name);
            stringstream atom_id;
            atom_id << atom_name << "_" << serial_number << "_" << id.str();
            assembly_atom->SetId(atom_id.str());
            TopologyAtom* topology_atom = (*it1);
            assembly_atom->MolecularDynamicAtom::SetCharge(topology_atom->GetAtomCharge() / CHARGE_DIVIDER);
            assembly_atom->MolecularDynamicAtom::SetMass(topology_atom->GetAtomMass());
            assembly_atom->MolecularDynamicAtom::SetAtomType(topology_atom->GetType());
            if(parameter != NULL)
            {
                if(atom_type_map.find(assembly_atom->GetAtomType()) != atom_type_map.end())
                {
                    ParameterFileAtom* parameter_atom = atom_type_map[assembly_atom->GetAtomType()];
                    assembly_atom->MolecularDynamicAtom::SetRadius(parameter_atom->GetRadius());
                }
                else
                {
                    assembly_atom->MolecularDynamicAtom::SetRadius(dNotSet);
                }
            }
            else
            {
                assembly_atom->MolecularDynamicAtom::SetRadius(dNotSet);
            }
            assembly_atom->SetResidue(assembly_residue);
            assembly_atom->SetName(topology_atom->GetAtomName());

            assembly_residue->AddAtom(assembly_atom);
        }
        this->AddResidue(assembly_residue);

    }
}

void Assembly::BuildAssemblyFromTopologyFile(TopologyFile *topology_file, string parameter_file)
{
    cout << "Building assembly from topology file ..." << endl;
    gmml::log(__LINE__, __FILE__, gmml::INF, "Building assembly from topology file ...");
    this->ClearAssembly();
    ParameterFile* parameter = NULL;
    ParameterFile::AtomTypeMap atom_type_map = ParameterFile::AtomTypeMap();
    if(parameter_file.compare("") != 0)
    {
        parameter = new ParameterFile(parameter_file);
        atom_type_map = parameter->GetAtomTypes();
    }
    name_ = topology_file->GetTitle();
    sequence_number_ = 1;
    int serial_number = 0;
    TopologyAssembly::TopologyResidueVector topology_residues = topology_file->GetAssembly()->GetResidues();
    for(TopologyAssembly::TopologyResidueVector::iterator it = topology_residues.begin(); it != topology_residues.end(); it++)
    {
        Residue* assembly_residue = new Residue();
        assembly_residue->SetAssembly(this);
        TopologyResidue* topology_residue = (*it);
        string residue_name = topology_residue->GetResidueName();
        assembly_residue->SetName(residue_name);
        stringstream id;
        id << residue_name << "_" << gmml::BLANK_SPACE << "_" << topology_residue->GetIndex()
           << "_" << gmml::BLANK_SPACE << "_" << gmml::BLANK_SPACE << "_" << id_;
        assembly_residue->SetId(id.str());

        TopologyResidue::TopologyAtomVector topology_atoms = topology_residue->GetAtoms();

        for(TopologyResidue::TopologyAtomVector::iterator it1 = topology_atoms.begin(); it1 != topology_atoms.end(); it1++)
        {
            serial_number++;
            Atom* assembly_atom = new Atom();
            string atom_name = (*it1)->GetAtomName();
            assembly_atom->SetName(atom_name);
            stringstream atom_id;
            atom_id << atom_name << "_" << serial_number << "_" << id.str();
            assembly_atom->SetId(atom_id.str());
            TopologyAtom* topology_atom = (*it1);
            assembly_atom->MolecularDynamicAtom::SetCharge(topology_atom->GetAtomCharge() / CHARGE_DIVIDER);
            assembly_atom->MolecularDynamicAtom::SetMass(topology_atom->GetAtomMass());
            assembly_atom->MolecularDynamicAtom::SetAtomType(topology_atom->GetType());
            if(parameter != NULL)
            {
                if(atom_type_map.find(assembly_atom->GetAtomType()) != atom_type_map.end())
                {
                    ParameterFileAtom* parameter_atom = atom_type_map[assembly_atom->GetAtomType()];
                    assembly_atom->MolecularDynamicAtom::SetRadius(parameter_atom->GetRadius());
                }
                else
                {
                    assembly_atom->MolecularDynamicAtom::SetRadius(dNotSet);
                }
            }
            else
            {
                assembly_atom->MolecularDynamicAtom::SetRadius(dNotSet);
            }
            assembly_atom->SetResidue(assembly_residue);
            assembly_atom->SetName(topology_atom->GetAtomName());

            assembly_residue->AddAtom(assembly_atom);
        }
        this->AddResidue(assembly_residue);

    }
}

void Assembly::BuildAssemblyFromLibraryFile(string library_file_path, string parameter_file)
{
    cout << "Building assembly from library file ..." << endl;
    gmml::log(__LINE__, __FILE__, gmml::INF, "Building assembly from library file ...");
    this->ClearAssembly();
    LibraryFile* library_file = new LibraryFile(library_file_path);
    ParameterFile* parameter = NULL;
    ParameterFile::AtomTypeMap atom_type_map = ParameterFile::AtomTypeMap();
    if(parameter_file.compare("") != 0)
    {
        parameter = new ParameterFile(parameter_file);
        atom_type_map = parameter->GetAtomTypes();
    }
    sequence_number_ = 1;
    LibraryFile::ResidueMap library_residues = library_file->GetResidues();
    stringstream ss;

    int sequence_number = 0;
    int serial_number = 0;
    for(LibraryFile::ResidueMap::iterator it = library_residues.begin(); it != library_residues.end(); it++)
    {
        sequence_number++;
        Residue* assembly_residue = new Residue();
        assembly_residue->SetAssembly(this);
        string residue_name = (*it).first;
        assembly_residue->SetName(residue_name);
        stringstream id;
        id << residue_name << "_" << gmml::BLANK_SPACE << "_" << sequence_number << "_" << gmml::BLANK_SPACE << "_"
           << gmml::BLANK_SPACE << "_" << id_;
        assembly_residue->SetId(id.str());
        LibraryFileResidue* library_residue = (*it).second;
        int lib_res_tail_atom_index = library_residue->GetTailAtomIndex();
        int lib_res_head_atom_index = library_residue->GetHeadAtomIndex();
        string library_residue_name = library_residue->GetName();
        if(distance(library_residues.begin(), it) == (int)library_residues.size()-1)
            ss << library_residue_name;
        else
            ss << library_residue_name << "-";

        LibraryFileResidue::AtomMap library_atoms = library_residue->GetAtoms();

        for(LibraryFileResidue::AtomMap::iterator it1 = library_atoms.begin(); it1 != library_atoms.end(); it1++)
        {
            serial_number++;
            Atom* assembly_atom = new Atom();
            LibraryFileAtom* library_atom = (*it1).second;
            string atom_name = library_atom->GetName();
            assembly_atom->SetName(atom_name);
            stringstream atom_id;
            atom_id << atom_name << "_" << library_atom->GetAtomOrder() << "_" << id.str();
            assembly_atom->SetId(atom_id.str());

            assembly_atom->SetResidue(assembly_residue);
            assembly_atom->SetName(library_atom->GetName());

            assembly_atom->MolecularDynamicAtom::SetCharge(library_atom->GetCharge());
            assembly_atom->MolecularDynamicAtom::SetAtomType(library_atom->GetType());
            if(parameter != NULL)
            {
                if(atom_type_map.find(assembly_atom->GetAtomType()) != atom_type_map.end())
                {
                    ParameterFileAtom* parameter_atom = atom_type_map[assembly_atom->GetAtomType()];
                    assembly_atom->MolecularDynamicAtom::SetMass(parameter_atom->GetMass());
                    assembly_atom->MolecularDynamicAtom::SetRadius(parameter_atom->GetRadius());
                }
                else
                {
                    assembly_atom->MolecularDynamicAtom::SetMass(dNotSet);
                    assembly_atom->MolecularDynamicAtom::SetRadius(dNotSet);
                }
            }
            else
            {
                assembly_atom->MolecularDynamicAtom::SetMass(dNotSet);
                assembly_atom->MolecularDynamicAtom::SetRadius(dNotSet);
            }

            Coordinate* coordinate = new Coordinate(library_atom->GetCoordinate());
            assembly_atom->AddCoordinate(coordinate);
            assembly_residue->AddAtom(assembly_atom);

            if(library_atom->GetAtomIndex() == lib_res_head_atom_index)
                assembly_residue->AddHeadAtom(assembly_atom);
            if(library_atom->GetAtomIndex() == lib_res_tail_atom_index)
                assembly_residue->AddTailAtom(assembly_atom);
        }
        residues_.push_back(assembly_residue);
    }
    name_ = ss.str();
}

void Assembly::BuildAssemblyFromLibraryFile(LibraryFile *library_file, string parameter_file)
{
    cout << "Building assembly from library file ..." << endl;
    gmml::log(__LINE__, __FILE__, gmml::INF, "Building assembly from library file ...");
    this->ClearAssembly();
    ParameterFile* parameter = NULL;
    ParameterFile::AtomTypeMap atom_type_map = ParameterFile::AtomTypeMap();
    if(parameter_file.compare("") != 0)
    {
        parameter = new ParameterFile(parameter_file);
        atom_type_map = parameter->GetAtomTypes();
    }
    sequence_number_ = 1;
    LibraryFile::ResidueMap library_residues = library_file->GetResidues();
    stringstream ss;

    int sequence_number = 0;
    int serial_number = 0;
    for(LibraryFile::ResidueMap::iterator it = library_residues.begin(); it != library_residues.end(); it++)
    {
        sequence_number++;
        Residue* assembly_residue = new Residue();
        assembly_residue->SetAssembly(this);
        string residue_name = (*it).first;
        assembly_residue->SetName(residue_name);
        stringstream id;
        id << residue_name << "_" << gmml::BLANK_SPACE << "_" << sequence_number << "_" << gmml::BLANK_SPACE << "_"
           << gmml::BLANK_SPACE << "_" << id_;
        assembly_residue->SetId(id.str());
        LibraryFileResidue* library_residue = (*it).second;
        int lib_res_tail_atom_index = library_residue->GetTailAtomIndex();
        int lib_res_head_atom_index = library_residue->GetHeadAtomIndex();
        string library_residue_name = library_residue->GetName();
        if(distance(library_residues.begin(), it) == (int)library_residues.size()-1)
            ss << library_residue_name;
        else
            ss << library_residue_name << "-";

        LibraryFileResidue::AtomMap library_atoms = library_residue->GetAtoms();

        for(LibraryFileResidue::AtomMap::iterator it1 = library_atoms.begin(); it1 != library_atoms.end(); it1++)
        {
            serial_number++;
            Atom* assembly_atom = new Atom();
            LibraryFileAtom* library_atom = (*it1).second;
            string atom_name = library_atom->GetName();
            assembly_atom->SetName(atom_name);
            stringstream atom_id;
            atom_id << atom_name << "_" << library_atom->GetAtomOrder() << "_" << id.str();
            assembly_atom->SetId(atom_id.str());

            assembly_atom->SetResidue(assembly_residue);
            assembly_atom->SetName(library_atom->GetName());

            assembly_atom->MolecularDynamicAtom::SetCharge(library_atom->GetCharge());
            assembly_atom->MolecularDynamicAtom::SetAtomType(library_atom->GetType());
            if(parameter != NULL)
            {
                if(atom_type_map.find(assembly_atom->GetAtomType()) != atom_type_map.end())
                {
                    ParameterFileAtom* parameter_atom = atom_type_map[assembly_atom->GetAtomType()];
                    assembly_atom->MolecularDynamicAtom::SetMass(parameter_atom->GetMass());
                    assembly_atom->MolecularDynamicAtom::SetRadius(parameter_atom->GetRadius());
                }
                else
                {
                    assembly_atom->MolecularDynamicAtom::SetMass(dNotSet);
                    assembly_atom->MolecularDynamicAtom::SetRadius(dNotSet);
                }
            }
            else
            {
                assembly_atom->MolecularDynamicAtom::SetMass(dNotSet);
                assembly_atom->MolecularDynamicAtom::SetRadius(dNotSet);
            }

            Coordinate* coordinate = new Coordinate(library_atom->GetCoordinate());
            assembly_atom->AddCoordinate(coordinate);
            assembly_residue->AddAtom(assembly_atom);

            if(library_atom->GetAtomIndex() == lib_res_head_atom_index)
                assembly_residue->AddHeadAtom(assembly_atom);
            if(library_atom->GetAtomIndex() == lib_res_tail_atom_index)
                assembly_residue->AddTailAtom(assembly_atom);
        }
        residues_.push_back(assembly_residue);
    }
    name_ = ss.str();
}

void Assembly::BuildAssemblyFromTopologyCoordinateFile(string topology_file_path, string coordinate_file_path, string parameter_file)
{
    cout << "Building assembly from topology and coordinate files ..." << endl;
    gmml::log(__LINE__, __FILE__, gmml::INF, "Building assembly from topology and coordinate files ...");
    this->ClearAssembly();
    TopologyFile* topology_file = new TopologyFile(topology_file_path);
    ParameterFile* parameter = NULL;
    ParameterFile::AtomTypeMap atom_type_map = ParameterFile::AtomTypeMap();
    if(parameter_file.compare("") != 0)
    {
        parameter = new ParameterFile(parameter_file);
        atom_type_map = parameter->GetAtomTypes();
    }
    name_ = topology_file->GetTitle();
    sequence_number_ = 1;
    int serial_number = 0;
    TopologyAssembly::TopologyResidueVector topology_residues = topology_file->GetAssembly()->GetResidues();
    for(TopologyAssembly::TopologyResidueVector::iterator it = topology_residues.begin(); it != topology_residues.end(); it++)
    {
        Residue* assembly_residue = new Residue();
        assembly_residue->SetAssembly(this);
        TopologyResidue* topology_residue = (*it);
        string residue_name = topology_residue->GetResidueName();
        assembly_residue->SetName(residue_name);
        stringstream id;
        id << residue_name << "_" << gmml::BLANK_SPACE << "_" << topology_residue->GetIndex() << "_" << gmml::BLANK_SPACE << "_"
           << gmml::BLANK_SPACE << "_" << id_;
        assembly_residue->SetId(id.str());

        TopologyResidue::TopologyAtomVector topology_atoms = topology_residue->GetAtoms();

        for(TopologyResidue::TopologyAtomVector::iterator it1 = topology_atoms.begin(); it1 != topology_atoms.end(); it1++)
        {
            serial_number++;
            Atom* assembly_atom = new Atom();
            string atom_name = (*it1)->GetAtomName();
            assembly_atom->SetName(atom_name);
            stringstream atom_id;
            atom_id << atom_name << "_" << serial_number << "_" << id.str();
            assembly_atom->SetId(atom_id.str());
            TopologyAtom* topology_atom = (*it1);

            assembly_atom->MolecularDynamicAtom::SetCharge(topology_atom->GetAtomCharge() / CHARGE_DIVIDER);
            assembly_atom->MolecularDynamicAtom::SetAtomType(topology_atom->GetType());
            assembly_atom->MolecularDynamicAtom::SetMass(topology_atom->GetAtomMass());
            if(parameter != NULL)
            {
                if(atom_type_map.find(assembly_atom->GetAtomType()) != atom_type_map.end())
                {
                    ParameterFileAtom* parameter_atom = atom_type_map[assembly_atom->GetAtomType()];
                    assembly_atom->MolecularDynamicAtom::SetRadius(parameter_atom->GetRadius());
                }
                else
                {
                    assembly_atom->MolecularDynamicAtom::SetRadius(dNotSet);
                }
            }
            else
            {
                assembly_atom->MolecularDynamicAtom::SetRadius(dNotSet);
            }

            int topology_atom_index = topology_atom->GetIndex();

            assembly_atom->SetResidue(assembly_residue);
            assembly_atom->SetName(topology_atom->GetAtomName());

            CoordinateFile* coordinate_file = new CoordinateFile(coordinate_file_path);
            vector<GeometryTopology::Coordinate*> coord_file_coordinates = coordinate_file->GetCoordinates();
            assembly_atom->AddCoordinate(coord_file_coordinates.at(topology_atom_index-1));
            assembly_residue->AddAtom(assembly_atom);
        }
        residues_.push_back(assembly_residue);
    }
}

void Assembly::BuildAssemblyFromTopologyCoordinateFile(TopologyFile *topology_file, CoordinateFile *coordinate_file, string parameter_file)
{
    cout << "Building assembly from topology and coordinate files ..." << endl;
    gmml::log(__LINE__, __FILE__, gmml::INF, "Building assembly from topology and coordinate files ...");
    this->ClearAssembly();
    ParameterFile* parameter = NULL;
    ParameterFile::AtomTypeMap atom_type_map = ParameterFile::AtomTypeMap();
    if(parameter_file.compare("") != 0)
    {
        parameter = new ParameterFile(parameter_file);
        atom_type_map = parameter->GetAtomTypes();
    }

    name_ = topology_file->GetTitle();
    sequence_number_ = 1;
    int serial_number = 0;
    TopologyAssembly::TopologyResidueVector topology_residues = topology_file->GetAssembly()->GetResidues();
    for(TopologyAssembly::TopologyResidueVector::iterator it = topology_residues.begin(); it != topology_residues.end(); it++)
    {
        Residue* assembly_residue = new Residue();
        assembly_residue->SetAssembly(this);
        TopologyResidue* topology_residue = (*it);
        string residue_name = topology_residue->GetResidueName();
        assembly_residue->SetName(residue_name);
        stringstream id;
        id << residue_name << "_" << gmml::BLANK_SPACE << "_" << topology_residue->GetIndex() << "_" << gmml::BLANK_SPACE << "_"
           << gmml::BLANK_SPACE << "_" << id_;
        assembly_residue->SetId(id.str());

        TopologyResidue::TopologyAtomVector topology_atoms = topology_residue->GetAtoms();

        for(TopologyResidue::TopologyAtomVector::iterator it1 = topology_atoms.begin(); it1 != topology_atoms.end(); it1++)
        {
            serial_number++;
            Atom* assembly_atom = new Atom();
            string atom_name = (*it1)->GetAtomName();
            assembly_atom->SetName(atom_name);
            stringstream atom_id;
            atom_id << atom_name << "_" << serial_number << id.str();
            assembly_atom->SetId(atom_id.str());
            TopologyAtom* topology_atom = (*it1);

            assembly_atom->MolecularDynamicAtom::SetCharge(topology_atom->GetAtomCharge() / CHARGE_DIVIDER);
            assembly_atom->MolecularDynamicAtom::SetAtomType(topology_atom->GetType());
            assembly_atom->MolecularDynamicAtom::SetMass(topology_atom->GetAtomMass());
            if(parameter != NULL)
            {
                if(atom_type_map.find(assembly_atom->GetAtomType()) != atom_type_map.end())
                {
                    ParameterFileAtom* parameter_atom = atom_type_map[assembly_atom->GetAtomType()];
                    assembly_atom->MolecularDynamicAtom::SetRadius(parameter_atom->GetRadius());
                }
                else
                {
                    assembly_atom->MolecularDynamicAtom::SetRadius(dNotSet);
                }
            }
            else
            {
                assembly_atom->MolecularDynamicAtom::SetRadius(dNotSet);
            }

            int topology_atom_index = topology_atom->GetIndex();

            assembly_atom->SetResidue(assembly_residue);
            assembly_atom->SetName(topology_atom->GetAtomName());

            vector<GeometryTopology::Coordinate*> coord_file_coordinates = coordinate_file->GetCoordinates();
            assembly_atom->AddCoordinate(coord_file_coordinates.at(topology_atom_index-1));
            assembly_residue->AddAtom(assembly_atom);
        }
        residues_.push_back(assembly_residue);
    }
}

void Assembly::BuildAssemblyFromPrepFile(string prep_file_path, string parameter_file)
{
    cout << "Building assembly from prep file ..." << endl;
    gmml::log(__LINE__, __FILE__, gmml::INF, "Building assembly from prep file ...");
    this->ClearAssembly();
    PrepFile* prep_file = new PrepFile(prep_file_path);
    ParameterFile* parameter = NULL;
    ParameterFile::AtomTypeMap atom_type_map = ParameterFile::AtomTypeMap();
    if(parameter_file.compare("") != 0)
    {
        parameter = new ParameterFile(parameter_file);
        atom_type_map = parameter->GetAtomTypes();
    }
    sequence_number_ = 1;
    PrepFile::ResidueMap prep_residues = prep_file->GetResidues();
    stringstream ss;

    int sequence_number = 0;
    int serial_number = 0;
    for(PrepFile::ResidueMap::iterator it = prep_residues.begin(); it != prep_residues.end(); it++)
    {
        sequence_number++;
        CoordinateVector cartesian_coordinate_list = CoordinateVector();
        int head_atom_index = INFINITY;
        int tail_atom_index = -INFINITY;
        Atom* head_atom = new Atom();
        Atom* tail_atom = new Atom();

        Residue* assembly_residue = new Residue();
        assembly_residue->SetAssembly(this);
        PrepFileResidue* prep_residue = (*it).second;
        string prep_residue_name = prep_residue->GetName();
        assembly_residue->SetName(prep_residue_name);
        stringstream id;
        id << prep_residue_name << "_" << gmml::BLANK_SPACE << "_" << sequence_number << "_" << gmml::BLANK_SPACE << "_"
           << gmml::BLANK_SPACE << "_" << id_;
        assembly_residue->SetId(id.str());
        if(distance(prep_residues.begin(), it) == (int)prep_residues.size()-1)
            ss << prep_residue_name;
        else
            ss << prep_residue_name << "-";
        PrepFileResidue::PrepFileAtomVector prep_atoms = prep_residue->GetAtoms();

        for(PrepFileResidue::PrepFileAtomVector::iterator it1 = prep_atoms.begin(); it1 != prep_atoms.end(); it1++)
        {
            serial_number++;
            Atom* assembly_atom = new Atom();
            PrepFileAtom* prep_atom = (*it1);
            assembly_atom->SetResidue(assembly_residue);
            string atom_name = prep_atom->GetName();
            assembly_atom->SetName(atom_name);
            stringstream atom_id;
            atom_id << atom_name << "_" << serial_number << "_" << id.str();
            assembly_atom->SetId(atom_id.str());

            assembly_atom->MolecularDynamicAtom::SetAtomType(prep_atom->GetType());
            assembly_atom->MolecularDynamicAtom::SetCharge(prep_atom->GetCharge());
            if(parameter != NULL)
            {
                if(atom_type_map.find(assembly_atom->GetAtomType()) != atom_type_map.end())
                {
                    ParameterFileAtom* parameter_atom = atom_type_map[assembly_atom->GetAtomType()];
                    assembly_atom->MolecularDynamicAtom::SetMass(parameter_atom->GetMass());
                    assembly_atom->MolecularDynamicAtom::SetRadius(parameter_atom->GetRadius());
                }
                else
                {
                    assembly_atom->MolecularDynamicAtom::SetMass(dNotSet);
                    assembly_atom->MolecularDynamicAtom::SetRadius(dNotSet);
                }
            }
            else
            {
                assembly_atom->MolecularDynamicAtom::SetMass(dNotSet);
                assembly_atom->MolecularDynamicAtom::SetRadius(dNotSet);
            }

            if(prep_residue->GetCoordinateType() == PrepFileSpace::kINT)
            {
                vector<Coordinate*> coordinate_list = vector<Coordinate*>();
                int index = distance(prep_atoms.begin(), it1);
                if(index == 0)
                {
                }
                if(index == 1)
                {
                    int parent_index = prep_atom->GetBondIndex() - 1;
                    Coordinate* parent_coordinate = cartesian_coordinate_list.at(parent_index);
                    coordinate_list.push_back(parent_coordinate);
                }
                if(index == 2)
                {
                    int grandparent_index = prep_atom->GetAngleIndex() - 1;
                    int parent_index = prep_atom->GetBondIndex() - 1;
                    Coordinate* grandparent_coordinate = cartesian_coordinate_list.at(grandparent_index);
                    Coordinate* parent_coordinate = cartesian_coordinate_list.at(parent_index);
                    coordinate_list.push_back(grandparent_coordinate);
                    coordinate_list.push_back(parent_coordinate);
                }
                if(index > 2)
                {
                    int great_grandparent_index = prep_atom->GetDihedralIndex() - 1;
                    int grandparent_index = prep_atom->GetAngleIndex() - 1;
                    int parent_index = prep_atom->GetBondIndex() - 1;

                    Coordinate* great_grandparent_coordinate = cartesian_coordinate_list.at(great_grandparent_index);
                    Coordinate* grandparent_coordinate = cartesian_coordinate_list.at(grandparent_index);
                    Coordinate* parent_coordinate = cartesian_coordinate_list.at(parent_index);
                    coordinate_list.push_back(great_grandparent_coordinate);
                    coordinate_list.push_back(grandparent_coordinate);
                    coordinate_list.push_back(parent_coordinate);
                }
                Coordinate* coordinate = new Coordinate();
                coordinate = gmml::ConvertInternalCoordinate2CartesianCoordinate(coordinate_list, prep_atom->GetBondLength(),
                                                                                 prep_atom->GetAngle(), prep_atom->GetDihedral());
                cartesian_coordinate_list.push_back(coordinate);

                assembly_atom->AddCoordinate(coordinate);
            }
            else if(prep_residue->GetCoordinateType() == PrepFileSpace::kXYZ)
            {
                assembly_atom->AddCoordinate(new Coordinate(prep_atom->GetBondLength(), prep_atom->GetAngle(), prep_atom->GetDihedral()));
            }
            if(prep_atom->GetTopologicalType() == kTopTypeM && prep_atom->GetType().compare(prep_residue->GetDummyAtomType()) != 0)
            {
                if(head_atom_index > prep_atom->GetIndex())
                {
                    head_atom_index = prep_atom->GetIndex();
                    head_atom = assembly_atom;
                }
                if(tail_atom_index < prep_atom->GetIndex())
                {
                    tail_atom_index = prep_atom->GetIndex();
                    tail_atom = assembly_atom;
                }
            }
            if(assembly_atom->GetAtomType().compare("DU") != 0)
                assembly_residue->AddAtom(assembly_atom);
        }
        assembly_residue->AddHeadAtom(head_atom);
        assembly_residue->AddTailAtom(tail_atom);
        residues_.push_back(assembly_residue);
    }
    name_ = ss.str();
}

void Assembly::BuildAssemblyFromPrepFile(PrepFile *prep_file, string parameter_file)
{
    cout << "Build assembly from prep file ..." << endl;
    gmml::log(__LINE__, __FILE__, gmml::INF, "Building assembly from prep file ...");
    this->ClearAssembly();
    ParameterFile* parameter = NULL;
    ParameterFile::AtomTypeMap atom_type_map = ParameterFile::AtomTypeMap();
    if(parameter_file.compare("") != 0)
    {
        parameter = new ParameterFile(parameter_file);
        atom_type_map = parameter->GetAtomTypes();
    }
    sequence_number_ = 1;
    PrepFile::ResidueMap prep_residues = prep_file->GetResidues();
    stringstream ss;

    int sequence_number = 0;
    int serial_number = 0;
    for(PrepFile::ResidueMap::iterator it = prep_residues.begin(); it != prep_residues.end(); it++)
    {
        sequence_number++;
        CoordinateVector cartesian_coordinate_list = CoordinateVector();
        int head_atom_index = INFINITY;
        int tail_atom_index = -INFINITY;
        Atom* head_atom = new Atom();
        Atom* tail_atom = new Atom();

        Residue* assembly_residue = new Residue();
        assembly_residue->SetAssembly(this);
        string residue_name = (*it).first;
        assembly_residue->SetName(residue_name);
        stringstream id;
        id << residue_name << "_" << gmml::BLANK_SPACE << "_" << sequence_number << "_" << gmml::BLANK_SPACE << "_"
           << gmml::BLANK_SPACE << "_" << id_;
        assembly_residue->SetId(id.str());
        PrepFileResidue* prep_residue = (*it).second;
        string prep_residue_name = prep_residue->GetName();
        if(distance(prep_residues.begin(), it) == (int)prep_residues.size()-1)
            ss << prep_residue_name;
        else
            ss << prep_residue_name << "-";
        PrepFileResidue::PrepFileAtomVector prep_atoms = prep_residue->GetAtoms();
        PrepFileResidue::PrepFileAtomVector parent_atoms = prep_residue->GetAtomsParentVector();

        for(PrepFileResidue::PrepFileAtomVector::iterator it1 = prep_atoms.begin(); it1 != prep_atoms.end(); it1++)
        {
            serial_number++;
            Atom* assembly_atom = new Atom();
            PrepFileAtom* prep_atom = (*it1);

            assembly_atom->SetResidue(assembly_residue);
            string atom_name = prep_atom->GetName();
            assembly_atom->SetName(atom_name);
            stringstream atom_id;
            atom_id << atom_name << "_" << serial_number << "_" << id.str();
            assembly_atom->SetId(atom_id.str());

            assembly_atom->MolecularDynamicAtom::SetAtomType(prep_atom->GetType());
            assembly_atom->MolecularDynamicAtom::SetCharge(prep_atom->GetCharge());
            if(parameter != NULL)
            {
                if(atom_type_map.find(assembly_atom->GetAtomType()) != atom_type_map.end())
                {
                    ParameterFileAtom* parameter_atom = atom_type_map[assembly_atom->GetAtomType()];
                    assembly_atom->MolecularDynamicAtom::SetMass(parameter_atom->GetMass());
                    assembly_atom->MolecularDynamicAtom::SetRadius(parameter_atom->GetRadius());
                }
                else
                {
                    assembly_atom->MolecularDynamicAtom::SetMass(dNotSet);
                    assembly_atom->MolecularDynamicAtom::SetRadius(dNotSet);
                }
            }
            else
            {
                assembly_atom->MolecularDynamicAtom::SetMass(dNotSet);
                assembly_atom->MolecularDynamicAtom::SetRadius(dNotSet);
            }

            if(prep_residue->GetCoordinateType() == PrepFileSpace::kINT)
            {
                vector<Coordinate*> coordinate_list = vector<Coordinate*>();
                int index = distance(prep_atoms.begin(), it1);
                if(index == 0)
                {

                }
                if(index == 1)
                {
                    int parent_index = parent_atoms.at(index)->GetIndex() - 1;
                    Coordinate* parent_coordinate = cartesian_coordinate_list.at(parent_index);
                    coordinate_list.push_back(parent_coordinate);
                }
                if(index == 2)
                {
                    int parent_index = parent_atoms.at(index)->GetIndex() - 1;
                    int grandparent_index = parent_atoms.at(parent_index)->GetIndex() - 1;
                    Coordinate* grandparent_coordinate = cartesian_coordinate_list.at(grandparent_index);
                    Coordinate* parent_coordinate = cartesian_coordinate_list.at(parent_index);
                    coordinate_list.push_back(grandparent_coordinate);
                    coordinate_list.push_back(parent_coordinate);
                }
                if(index > 2)
                {
                    int parent_index = parent_atoms.at(index)->GetIndex() - 1;
                    int grandparent_index = parent_atoms.at(parent_index)->GetIndex() - 1;
                    int great_grabdparent_index = parent_atoms.at(grandparent_index)->GetIndex() - 1;
                    Coordinate* great_grandparent_coordinate = cartesian_coordinate_list.at(great_grabdparent_index);
                    Coordinate* grandparent_coordinate = cartesian_coordinate_list.at(grandparent_index);
                    Coordinate* parent_coordinate = cartesian_coordinate_list.at(parent_index);
                    coordinate_list.push_back(great_grandparent_coordinate);
                    coordinate_list.push_back(grandparent_coordinate);
                    coordinate_list.push_back(parent_coordinate);
                }
                Coordinate* coordinate = gmml::ConvertInternalCoordinate2CartesianCoordinate(coordinate_list, prep_atom->GetBondLength(),
                                                                                             prep_atom->GetAngle(), prep_atom->GetDihedral());
                cartesian_coordinate_list.push_back(coordinate);

                assembly_atom->AddCoordinate(coordinate);
            }
            else if(prep_residue->GetCoordinateType() == PrepFileSpace::kXYZ)
            {
                assembly_atom->AddCoordinate(new Coordinate(prep_atom->GetBondLength(), prep_atom->GetAngle(), prep_atom->GetDihedral()));
            }
            if(prep_atom->GetTopologicalType() == kTopTypeM && prep_atom->GetType().compare(prep_residue->GetDummyAtomType()) != 0)
            {
                if(head_atom_index > prep_atom->GetIndex())
                {
                    head_atom_index = prep_atom->GetIndex();
                    head_atom = assembly_atom;
                }
                if(tail_atom_index < prep_atom->GetIndex())
                {
                    tail_atom_index = prep_atom->GetIndex();
                    tail_atom = assembly_atom;
                }
            }
            if(assembly_atom->GetAtomType().compare("DU") != 0)
                assembly_residue->AddAtom(assembly_atom);
        }
        assembly_residue->AddHeadAtom(head_atom);
        assembly_residue->AddTailAtom(tail_atom);
        residues_.push_back(assembly_residue);
    }
    name_ = ss.str();
}

