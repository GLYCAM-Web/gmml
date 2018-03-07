#include <math.h>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <set>
#include <queue>
#include <stack>

#include "../../../includes/InputSet/PdbFileSpace/inputfile.hpp"
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
#include "../../../includes/InputSet/TopologyFileSpace/topologyfileprocessingexception.hpp"
#include "../../../includes/InputSet/CoordinateFileSpace/coordinatefile.hpp"
#include "../../../includes/InputSet/CoordinateFileSpace/coordinatefileprocessingexception.hpp"
#include "../../../includes/ParameterSet/PrepFileSpace/prepfile.hpp"
#include "../../../includes/ParameterSet/PrepFileSpace/prepfileresidue.hpp"
#include "../../../includes/ParameterSet/PrepFileSpace/prepfileatom.hpp"
#include "../../../includes/ParameterSet/PrepFileSpace/prepfileprocessingexception.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbfile.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbtitlesection.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbmodelsection.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbmodelcard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbmodelresidueset.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbatomsection.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbheterogenatomsection.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbatomcard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbconnectsection.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdblinksection.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdblinkcard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdblinkcardresidue.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbfileprocessingexception.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbremarksection.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbmastercard.hpp"
#include "../../../includes/InputSet/PdbqtFileSpace/pdbqtfile.hpp"
#include "../../../includes/InputSet/PdbqtFileSpace/pdbqtatom.hpp"
#include "../../../includes/InputSet/PdbqtFileSpace/pdbqtmodel.hpp"
#include "../../../includes/InputSet/PdbqtFileSpace/pdbqtmodelcard.hpp"
#include "../../../includes/InputSet/PdbqtFileSpace/pdbqtatomcard.hpp"
#include "../../../includes/InputSet/PdbqtFileSpace/pdbqtmodelresidueset.hpp"
#include "../../../includes/InputSet/PdbqtFileSpace/pdbqtfileprocessingexception.hpp"
#include "../../../includes/ParameterSet/LibraryFileSpace/libraryfile.hpp"
#include "../../../includes/ParameterSet/LibraryFileSpace/libraryfileatom.hpp"
#include "../../../includes/ParameterSet/LibraryFileSpace/libraryfileresidue.hpp"
#include "../../../includes/ParameterSet/LibraryFileSpace/libraryfileprocessingexception.hpp"
#include "../../../includes/ParameterSet/PrepFileSpace/prepfile.hpp"
#include "../../../includes/ParameterSet/PrepFileSpace/prepfileresidue.hpp"
#include "../../../includes/ParameterSet/PrepFileSpace/prepfileatom.hpp"
#include "../../../includes/ParameterSet/ParameterFileSpace/parameterfile.hpp"
#include "../../../includes/ParameterSet/ParameterFileSpace/parameterfilebond.hpp"
#include "../../../includes/ParameterSet/ParameterFileSpace/parameterfileangle.hpp"
#include "../../../includes/ParameterSet/ParameterFileSpace/parameterfiledihedral.hpp"
#include "../../../includes/ParameterSet/ParameterFileSpace/parameterfiledihedralterm.hpp"
#include "../../../includes/ParameterSet/ParameterFileSpace/parameterfileatom.hpp"
#include "../../../includes/utils.hpp"
#include "../../../includes/common.hpp"
#include "../../../includes/GeometryTopology/grid.hpp"
#include "../../../includes/GeometryTopology/cell.hpp"

#include <unistd.h>
#include <errno.h>
#include <string.h>
#include <iostream>

using MolecularModeling::Assembly;

//////////////////////////////////////////////////////////
//                       FUNCTIONS                      //
//////////////////////////////////////////////////////////
bool Assembly::CheckCondensedSequenceSanity(std::string sequence, CondensedSequenceSpace::CondensedSequence::CondensedSequenceGlycam06ResidueTree& prep_residues)
{
    try
    {
        CondensedSequenceSpace::CondensedSequence* condensed_sequence = new CondensedSequenceSpace::CondensedSequence(sequence);
        prep_residues = condensed_sequence->GetCondensedSequenceGlycam06ResidueTree();
        for(CondensedSequenceSpace::CondensedSequence::CondensedSequenceGlycam06ResidueTree::iterator it = prep_residues.begin(); it != prep_residues.end(); it++)
        {
            CondensedSequenceSpace::CondensedSequenceGlycam06Residue* glycam06_residue = *it;
            std::string glycam06_residue_name = glycam06_residue->GetName();
            if(glycam06_residue_name.compare("UNK") == 0)
            {
                std::cout << "The input sequence (" << sequence << ") is not valid" << std::endl;
                return false;
            }
        }
    }
    catch(std::exception ex)
    {
        std::cout << "The input sequence (" << sequence << ") is not valid" << std::endl;
        return false;
    }

    std::cout << "The input sequence (" << sequence << ") is valid" << std::endl;
    return true;
}

void Assembly::BuildAssemblyFromCondensedSequence(std::string sequence, std::string prep_file, std::string parameter_file, bool structure)
{
    try
    {
        CondensedSequenceSpace::CondensedSequence* condensed_sequence = new CondensedSequenceSpace::CondensedSequence(sequence);
        CondensedSequenceSpace::CondensedSequence::CondensedSequenceGlycam06ResidueTree glycam06_residues = condensed_sequence->GetCondensedSequenceGlycam06ResidueTree();
        PrepFileSpace::PrepFile* prep = new PrepFileSpace::PrepFile(prep_file);
        PrepFileSpace::PrepFile::ResidueMap prep_residue_map = prep->GetResidues();
        ParameterFileSpace::ParameterFile* parameter = NULL;
        ParameterFileSpace::ParameterFile::AtomTypeMap atom_type_map = ParameterFileSpace::ParameterFile::AtomTypeMap();
        if(parameter_file.compare("")!= 0)
        {
            parameter = new ParameterFileSpace::ParameterFile(parameter_file);
            atom_type_map = parameter->GetAtomTypes();
        }
        ResidueVector parent_residues = ResidueVector();
        ResidueVector branch_residues = ResidueVector();
        std::vector<bool> derivatives = std::vector<bool>();
        int sequence_number = 0;
        int serial_number = 0;
        std::stringstream ss;
        for(CondensedSequenceSpace::CondensedSequence::CondensedSequenceGlycam06ResidueTree::iterator it = glycam06_residues.begin(); it != glycam06_residues.end(); ++it)
        {
            CondensedSequenceSpace::CondensedSequenceGlycam06Residue* glycam06_residue = *it;
            std::string glycam06_residue_name = glycam06_residue->GetName();
            std::string glycam06_residue_parent_oxygen = glycam06_residue->GetParentOxygen();

            if(prep_residue_map.find(glycam06_residue_name) != prep_residue_map.end())
            {
                PrepFileSpace::PrepFileResidue* prep_residue = prep_residue_map[glycam06_residue_name];

                // Build residue from prep residue
                sequence_number++;
                CoordinateVector cartesian_coordinate_list = CoordinateVector();
                Residue* assembly_residue = new Residue();
                assembly_residue->SetAssembly(this);
                std::string prep_residue_name = prep_residue->GetName();
                assembly_residue->SetName(prep_residue_name);
                std::stringstream id;
                id << prep_residue_name << "_" << gmml::BLANK_SPACE << "_" << sequence_number << "_" << gmml::BLANK_SPACE << "_"
                   << gmml::BLANK_SPACE << "_" << id_;

                assembly_residue->SetId(id.str());
                if(std::distance(glycam06_residues.begin(), it) == (int)glycam06_residues.size()-1)
                    ss << prep_residue_name;
                else
                    ss << prep_residue_name << "-";
                PrepFileSpace::PrepFileResidue::PrepFileAtomVector prep_atoms = prep_residue->GetAtoms();
                for(PrepFileSpace::PrepFileResidue::PrepFileAtomVector::iterator it1 = prep_atoms.begin(); it1 != prep_atoms.end(); it1++)
                {
                    PrepFileSpace::PrepFileAtom* prep_atom = (*it1);
                    std::string atom_name = prep_atom->GetName();
                    if(prep_atom->GetType() != "DU")
                        serial_number++;

                    Atom* assembly_atom = new Atom();
                    assembly_atom->SetResidue(assembly_residue);
                    assembly_atom->SetName(atom_name);
                    std::stringstream atom_id;
                    atom_id << atom_name << "_" << serial_number << "_" << id.str();
                    assembly_atom->SetId(atom_id.str());

                    assembly_atom->MolecularDynamicAtom::SetAtomType(prep_atom->GetType());
                    assembly_atom->MolecularDynamicAtom::SetCharge(prep_atom->GetCharge());
                    if(parameter != NULL)
                    {
                        if(atom_type_map.find(assembly_atom->GetAtomType()) != atom_type_map.end())
                        {
                            ParameterFileSpace::ParameterFileAtom* parameter_atom = atom_type_map[assembly_atom->GetAtomType()];
                            assembly_atom->MolecularDynamicAtom::SetMass(parameter_atom->GetMass());
                            assembly_atom->MolecularDynamicAtom::SetRadius(parameter_atom->GetRadius());
                        }
                        else
                        {
                            assembly_atom->MolecularDynamicAtom::SetMass(gmml::dNotSet);
                            assembly_atom->MolecularDynamicAtom::SetRadius(gmml::dNotSet);
                        }
                    }
                    else
                    {
                        assembly_atom->MolecularDynamicAtom::SetMass(gmml::dNotSet);
                        assembly_atom->MolecularDynamicAtom::SetRadius(gmml::dNotSet);
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
                    else if(atom_name.find("C1") != std::string::npos && prep_residue_name.compare("ACX") == 0)
                    {
                        assembly_residue->AddHeadAtom(assembly_atom);
                    }

                    if(prep_residue->GetCoordinateType() == PrepFileSpace::kINT)
                    {
                        std::vector<GeometryTopology::Coordinate*> coordinate_list = std::vector<GeometryTopology::Coordinate*>();
                        int index = std::distance(prep_atoms.begin(), it1);
                        if(index == 0)
                        {
                        }
                        if(index == 1)
                        {
                            int parent_index = prep_atom->GetBondIndex() - 1;
                            GeometryTopology::Coordinate* parent_coordinate = cartesian_coordinate_list.at(parent_index);
                            coordinate_list.push_back(parent_coordinate);
                        }
                        if(index == 2)
                        {
                            int grandparent_index = prep_atom->GetAngleIndex() - 1;
                            int parent_index = prep_atom->GetBondIndex() - 1;
                            GeometryTopology::Coordinate* grandparent_coordinate = cartesian_coordinate_list.at(grandparent_index);
                            GeometryTopology::Coordinate* parent_coordinate = cartesian_coordinate_list.at(parent_index);
                            coordinate_list.push_back(grandparent_coordinate);
                            coordinate_list.push_back(parent_coordinate);
                        }
                        if(index > 2)
                        {
                            int great_grandparent_index = prep_atom->GetDihedralIndex() - 1;
                            int grandparent_index = prep_atom->GetAngleIndex() - 1;
                            int parent_index = prep_atom->GetBondIndex() - 1;

                            GeometryTopology::Coordinate* great_grandparent_coordinate = cartesian_coordinate_list.at(great_grandparent_index);
                            GeometryTopology::Coordinate* grandparent_coordinate = cartesian_coordinate_list.at(grandparent_index);
                            GeometryTopology::Coordinate* parent_coordinate = cartesian_coordinate_list.at(parent_index);
                            coordinate_list.push_back(great_grandparent_coordinate);
                            coordinate_list.push_back(grandparent_coordinate);
                            coordinate_list.push_back(parent_coordinate);
                        }
                        GeometryTopology::Coordinate* coordinate = new GeometryTopology::Coordinate();
                        coordinate = gmml::ConvertInternalCoordinate2CartesianCoordinate(coordinate_list, prep_atom->GetBondLength(),
                                                                                         prep_atom->GetAngle(), prep_atom->GetDihedral());
                        cartesian_coordinate_list.push_back(coordinate);

                        assembly_atom->AddCoordinate(coordinate);
                    }
                    else if(prep_residue->GetCoordinateType() == PrepFileSpace::kXYZ)
                    {
                        assembly_atom->AddCoordinate(new GeometryTopology::Coordinate(prep_atom->GetBondLength(), prep_atom->GetAngle(), prep_atom->GetDihedral()));
                    }
                    if(assembly_atom->GetAtomType().compare("DU") != 0)
                        assembly_residue->AddAtom(assembly_atom);
                    if(atom_name.compare(glycam06_residue->GetAnomericCarbon()) == 0)
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
                if(glycam06_residue->GetParentId() != -1)
                {
                    Residue* parent_residue = residues_.at(glycam06_residue->GetParentId());
                    AtomVector parent_residue_atoms = parent_residue->GetAtoms();
                    for(AtomVector::iterator it3 = parent_residue_atoms.begin(); it3 != parent_residue_atoms.end(); it3++)
                    {
                        Atom* parent_atom = *it3;
                        if(parent_atom->GetName().compare(glycam06_residue->GetParentOxygen()) == 0)
                            parent_residue->AddTailAtom(parent_atom);
                    }
                    parent_residues.push_back(parent_residue);
                    branch_residues.push_back(assembly_residue);
                    if(glycam06_residue->GetIsDerivative())
                        derivatives.push_back(true);
                    else
                        derivatives.push_back(false);
                }
            }
            else
            {
                std::cout << "Residue " << glycam06_residue_name << " has not been found in the database" << std::endl;
            }
        }

        name_ = ss.str();
        this->SetSourceFile(prep_file);

        if(structure)
        {
            std::map<Residue*, int> parent_branch_map = std::map<Residue*, int>();
            int linkage_index = -1;
            for(ResidueVector::iterator it = parent_residues.begin(); it != parent_residues.end(); it++)
            {
                Residue* parent_residue = (*it);
                int parent_index = std::distance(parent_residues.begin(), it);
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
                    std::stringstream linkage_name;
                    linkage_name << assembly_residue->GetName() << assembly_residue->GetHeadAtoms().at(0)->GetName().at(1) << "-" <<
                                    parent_residue->GetTailAtoms().at(branch_index)->GetName().at(1) << parent_residue->GetName();
                    if(linkage_name.str().find("0SA2-8") != std::string::npos ||
                                                linkage_name.str().find("0SB2-8") != std::string::npos)
                    {
                        this->SetPhiTorsion(assembly_residue, parent_residue, branch_index, gmml::EXTERNAL28LINKAGEROTAMERS[0][0]);
                        this->SetPsiTorsion(assembly_residue, parent_residue, branch_index, gmml::EXTERNAL28LINKAGEROTAMERS[0][1], false);
                        this->SetOmegaTorsion(assembly_residue, parent_residue, branch_index, gmml::EXTERNAL28LINKAGEROTAMERS[0][2], 7);
                        this->SetOmegaTorsion(assembly_residue, parent_residue, branch_index, gmml::EXTERNAL28LINKAGEROTAMERS[0][3], 8);
                        this->SetOmegaTorsion(assembly_residue, parent_residue, branch_index, gmml::EXTERNAL28LINKAGEROTAMERS[0][4], 9);
                    }
                    if(linkage_name.str().find("0GL2-8") != std::string::npos)
                    {
                        this->SetPhiTorsion(assembly_residue, parent_residue, branch_index, gmml::INTERNAL28LINKAGEROTAMERS[0][0]);
                        this->SetPsiTorsion(assembly_residue, parent_residue, branch_index, gmml::INTERNAL28LINKAGEROTAMERS[0][1], false);
                        this->SetOmegaTorsion(assembly_residue, parent_residue, branch_index, gmml::INTERNAL28LINKAGEROTAMERS[0][2], 7);
                        this->SetOmegaTorsion(assembly_residue, parent_residue, branch_index, gmml::INTERNAL28LINKAGEROTAMERS[0][3], 8);
                        this->SetOmegaTorsion(assembly_residue, parent_residue, branch_index, gmml::INTERNAL28LINKAGEROTAMERS[0][4], 9);
                    }
                }
                linkage_index++;
            }
        }
    }
    catch(std::exception ex)
    {
        std::cout << "Building assembly from " << sequence << " failed." << std::endl;
    }


}

Assembly::AssemblyVector Assembly::BuildAllRotamersFromCondensedSequence(std::string sequence, std::string prep_file, std::string parameter_file,
                                                                         CondensedSequenceSpace::CondensedSequence::CondensedSequenceRotamersAndGlycosidicAnglesInfo rotamers_glycosidic_angles_info,
                                                                         CondensedSequenceSpace::CondensedSequence::IndexNameMap& names)
{

    try
    {
        CondensedSequenceSpace::CondensedSequence* condensed_sequence = new CondensedSequenceSpace::CondensedSequence(sequence);
        AssemblyVector structures = AssemblyVector(condensed_sequence->CountAllPossible28LinkagesRotamers(rotamers_glycosidic_angles_info) *
                                                   condensed_sequence->CountAllPossibleSelectedRotamers(rotamers_glycosidic_angles_info));
        CondensedSequenceSpace::CondensedSequence::IndexLinkageConfigurationMap structure_map = condensed_sequence->CreateIndexLinkageConfigurationMap(
                    rotamers_glycosidic_angles_info, names);
        CondensedSequenceSpace::CondensedSequence::CondensedSequenceGlycam06ResidueTree glycam06_residues = condensed_sequence->GetCondensedSequenceGlycam06ResidueTree();
        PrepFileSpace::PrepFile* prep = new PrepFileSpace::PrepFile(prep_file);
        PrepFileSpace::PrepFile::ResidueMap prep_residue_map = prep->GetResidues();
        ParameterFileSpace::ParameterFile* parameter = NULL;
        ParameterFileSpace::ParameterFile::AtomTypeMap atom_type_map = ParameterFileSpace::ParameterFile::AtomTypeMap();
        if(parameter_file.compare("") != 0)
        {
            parameter = new ParameterFileSpace::ParameterFile(parameter_file);
            atom_type_map = parameter->GetAtomTypes();
        }
        for(unsigned int i = 0; i < structures.size(); i++)
        {
            ResidueVector parent_residues = ResidueVector();
            ResidueVector branch_residues = ResidueVector();
            std::vector<bool> derivatives = std::vector<bool>();
            structures.at(i) = new Assembly();
            int sequence_number = 0;
            int serial_number = 0;
            std::stringstream ss;
            for(CondensedSequenceSpace::CondensedSequence::CondensedSequenceGlycam06ResidueTree::iterator it = glycam06_residues.begin(); it != glycam06_residues.end(); ++it)
            {
                CondensedSequenceSpace::CondensedSequenceGlycam06Residue* glycam06_residue = *it;
                std::string glycam06_residue_name = glycam06_residue->GetName();
                std::string glycam06_residue_parent_oxygen = glycam06_residue->GetParentOxygen();

                if(prep_residue_map.find(glycam06_residue_name) != prep_residue_map.end())
                {
                    PrepFileSpace::PrepFileResidue* prep_residue = prep_residue_map[glycam06_residue_name];

                    // Build residue from prep residue
                    sequence_number++;
                    CoordinateVector cartesian_coordinate_list = CoordinateVector();

                    Residue* assembly_residue = new Residue();
                    assembly_residue->SetAssembly(structures.at(i));
                    std::string prep_residue_name = prep_residue->GetName();
                    assembly_residue->SetName(prep_residue_name);
                    std::stringstream id;
                    id << prep_residue_name << "_" << gmml::BLANK_SPACE << "_" << sequence_number << "_" << gmml::BLANK_SPACE << "_"
                       << gmml::BLANK_SPACE << "_" << id_;
                    assembly_residue->SetId(id.str());
                    if(std::distance(glycam06_residues.begin(), it) == (int)glycam06_residues.size()-1)
                        ss << prep_residue_name;
                    else
                        ss << prep_residue_name << "-";

                    PrepFileSpace::PrepFileResidue::PrepFileAtomVector prep_atoms = prep_residue->GetAtoms();
                    for(PrepFileSpace::PrepFileResidue::PrepFileAtomVector::iterator it1 = prep_atoms.begin(); it1 != prep_atoms.end(); it1++)
                    {
                        Atom* assembly_atom = new Atom();
                        PrepFileSpace::PrepFileAtom* prep_atom = (*it1);
                        if(prep_atom->GetType() != "DU")
                            serial_number++;
                        assembly_atom->SetResidue(assembly_residue);
                        std::string atom_name = prep_atom->GetName();
                        assembly_atom->SetName(atom_name);
                        std::stringstream atom_id;
                        atom_id << atom_name << "_" << serial_number << "_" << id.str();
                        assembly_atom->SetId(atom_id.str());

                        assembly_atom->MolecularDynamicAtom::SetAtomType(prep_atom->GetType());
                        assembly_atom->MolecularDynamicAtom::SetCharge(prep_atom->GetCharge());
                        if(parameter != NULL)
                        {
                            if(atom_type_map.find(assembly_atom->GetAtomType()) != atom_type_map.end())
                            {
                                ParameterFileSpace::ParameterFileAtom* parameter_atom = atom_type_map[assembly_atom->GetAtomType()];
                                assembly_atom->MolecularDynamicAtom::SetMass(parameter_atom->GetMass());
                                assembly_atom->MolecularDynamicAtom::SetRadius(parameter_atom->GetRadius());
                            }
                            else
                            {
                                assembly_atom->MolecularDynamicAtom::SetMass(gmml::dNotSet);
                                assembly_atom->MolecularDynamicAtom::SetRadius(gmml::dNotSet);
                            }
                        }
                        else
                        {
                            assembly_atom->MolecularDynamicAtom::SetMass(gmml::dNotSet);
                            assembly_atom->MolecularDynamicAtom::SetRadius(gmml::dNotSet);
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
                        else if(atom_name.find("C1") != std::string::npos && prep_residue_name.compare("ACX") == 0)
                        {
                            assembly_residue->AddHeadAtom(assembly_atom);
                        }

                        if(prep_residue->GetCoordinateType() == PrepFileSpace::kINT)
                        {
                            std::vector<GeometryTopology::Coordinate*> coordinate_list = std::vector<GeometryTopology::Coordinate*>();
                            int index = std::distance(prep_atoms.begin(), it1);
                            if(index == 0)
                            {
                            }
                            if(index == 1)
                            {
                                int parent_index = prep_atom->GetBondIndex() - 1;
                                GeometryTopology::Coordinate* parent_coordinate = cartesian_coordinate_list.at(parent_index);
                                coordinate_list.push_back(parent_coordinate);
                            }
                            if(index == 2)
                            {
                                int grandparent_index = prep_atom->GetAngleIndex() - 1;
                                int parent_index = prep_atom->GetBondIndex() - 1;
                                GeometryTopology::Coordinate* grandparent_coordinate = cartesian_coordinate_list.at(grandparent_index);
                                GeometryTopology::Coordinate* parent_coordinate = cartesian_coordinate_list.at(parent_index);
                                coordinate_list.push_back(grandparent_coordinate);
                                coordinate_list.push_back(parent_coordinate);
                            }
                            if(index > 2)
                            {
                                int great_grandparent_index = prep_atom->GetDihedralIndex() - 1;
                                int grandparent_index = prep_atom->GetAngleIndex() - 1;
                                int parent_index = prep_atom->GetBondIndex() - 1;

                                GeometryTopology::Coordinate* great_grandparent_coordinate = cartesian_coordinate_list.at(great_grandparent_index);
                                GeometryTopology::Coordinate* grandparent_coordinate = cartesian_coordinate_list.at(grandparent_index);
                                GeometryTopology::Coordinate* parent_coordinate = cartesian_coordinate_list.at(parent_index);
                                coordinate_list.push_back(great_grandparent_coordinate);
                                coordinate_list.push_back(grandparent_coordinate);
                                coordinate_list.push_back(parent_coordinate);
                            }
                            GeometryTopology::Coordinate* coordinate = new GeometryTopology::Coordinate();
                            coordinate = gmml::ConvertInternalCoordinate2CartesianCoordinate(coordinate_list, prep_atom->GetBondLength(),
                                                                                             prep_atom->GetAngle(), prep_atom->GetDihedral());
                            cartesian_coordinate_list.push_back(coordinate);

                            assembly_atom->AddCoordinate(coordinate);
                        }
                        else if(prep_residue->GetCoordinateType() == PrepFileSpace::kXYZ)
                        {
                            assembly_atom->AddCoordinate(new GeometryTopology::Coordinate(prep_atom->GetBondLength(), prep_atom->GetAngle(), prep_atom->GetDihedral()));
                        }
                        if(assembly_atom->GetAtomType().compare("DU") != 0)
                            assembly_residue->AddAtom(assembly_atom);
                        if(atom_name.compare(glycam06_residue->GetAnomericCarbon()) == 0)
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
                    if(glycam06_residue->GetParentId() != -1)
                    {
                        Residue* parent_residue = structures.at(i)->residues_.at(glycam06_residue->GetParentId());
                        AtomVector parent_residue_atoms = parent_residue->GetAtoms();
                        for(AtomVector::iterator it3 = parent_residue_atoms.begin(); it3 != parent_residue_atoms.end(); it3++)
                        {
                            Atom* parent_atom = *it3;
                            if(parent_atom->GetName().compare(glycam06_residue->GetParentOxygen()) == 0)
                                parent_residue->AddTailAtom(parent_atom);
                        }
                        parent_residues.push_back(parent_residue);
                        branch_residues.push_back(assembly_residue);
                        if(glycam06_residue->GetIsDerivative())
                            derivatives.push_back(true);
                        else
                            derivatives.push_back(false);
                    }
                }
                else
                {
                    std::cout << "Residue " << glycam06_residue_name << " has not been found in the database" << std::endl;
                }
            }

            structures.at(i)->name_ = ss.str();
            structures.at(i)->SetSourceFile(prep_file);

            if(true)
            {
                std::map<Residue*, int> parent_branch_map = std::map<Residue*, int>();
                int linkage_index = -1;
                for(ResidueVector::iterator it = parent_residues.begin(); it != parent_residues.end(); it++)
                {
                    Residue* parent_residue = (*it);
                    int parent_index = std::distance(parent_residues.begin(), it);
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
                            /*std::cout << rotamers_glycosidic_angles_info.at(linkage_index).first << " "
                             << assembly_residue->GetName() << "(" << assembly_residue->GetHeadAtoms().at(0)->GetName() << "-"
                             << parent_residue->GetTailAtoms().at(branch_index)->GetName() << ")" << parent_residue->GetName() << std::endl;*/
                            std::vector<double> phi_psi_omega = structure_map[i].at(linkage_index);
                            if(phi_psi_omega.at(0) != gmml::dNotSet)
                            {
                                double phi = phi_psi_omega.at(0);
                                structures.at(i)->SetPhiTorsion(assembly_residue, parent_residue, branch_index, phi);// Set phi angle of assembly_residue-parent_residue to phi
                            }
                            if(phi_psi_omega.at(1) != gmml::dNotSet)
                            {
                                double psi = phi_psi_omega.at(1);
                                structures.at(i)->SetPsiTorsion(assembly_residue, parent_residue, branch_index, psi);// Set psi angle of assembly_residue-parent_residue to psi
                            }
                            if(phi_psi_omega.at(2) != gmml::dNotSet)
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
                            std::vector<double> phi_psi_omega = structure_map[i].at(linkage_index);
                            if(phi_psi_omega.at(2) != gmml::dNotSet)
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
    catch(std::exception ex)
    {
        std::cout << "Building assembly from " << sequence << " failed." << std::endl;
    }
}

/** ***************************************************************************
 *  *          Build from PDB files
 *   * ************************************************************************* **/

void Assembly::BuildAssemblyFromPdbFile(std::string pdb_file_path, std::vector<std::string> amino_lib_files, std::vector<std::string> glycam_lib_files,
                                        std::vector<std::string> other_lib_files, std::vector<std::string> prep_files, std::string parameter_file)
{
    std::cout << "Building assembly from pdb file ..." << std::endl;
//    std::cout << "Reading PDB file into PdbFileSpace::PdbFile structure." << std::endl;
    PdbFileSpace::PdbFile pdb_file;
    try
    {
        gmml::log(__LINE__, __FILE__, gmml::INF, "Reading PDB file into PdbFileSpace::PdbFile structure ...");
        pdb_file = PdbFileSpace::PdbFile(pdb_file_path);
    }
    catch(PdbFileSpace::PdbFileProcessingException &ex)
    {
        std::cout << "Generating PdbFileSpace::PdbFile structure from " << pdb_file_path << "failed." << std::endl;
    }
    this->BuildAssemblyFromPdbFile(&pdb_file, amino_lib_files, glycam_lib_files, other_lib_files, prep_files, parameter_file);
}


void Assembly::BuildAssemblyFromPdbFile(PdbFileSpace::PdbFile *pdb_file, std::vector<std::string> amino_lib_files, std::vector<std::string> glycam_lib_files,
                                        std::vector<std::string> other_lib_files, std::vector<std::string> prep_files, std::string parameter_file)
{
//    std::cout << "Building assembly from pdb file ..." << std::endl;
    gmml::log(__LINE__, __FILE__, gmml::INF, "Building assembly from pdb file ...");
    try
    {
        this->ClearAssembly();
        gmml::log(__LINE__, __FILE__, gmml::INF, "Assembly cleared ...");
        // this->pdb_file_ = pdb_file;
        ParameterFileSpace::ParameterFile* parameter = NULL;
        gmml::log(__LINE__, __FILE__, gmml::INF, "Parameter File Created ...");
        ParameterFileSpace::ParameterFile::AtomTypeMap atom_type_map = ParameterFileSpace::ParameterFile::AtomTypeMap();
        if(parameter_file.compare("") != 0)
        {
            parameter = new ParameterFileSpace::ParameterFile(parameter_file);
            atom_type_map = parameter->GetAtomTypes();
        }

        LibraryFileSpace::LibraryFile::ResidueMap lib_residues = LibraryFileSpace::LibraryFile::ResidueMap();
        PrepFileSpace::PrepFile::ResidueMap prep_residues = PrepFileSpace::PrepFile::ResidueMap();
        std::vector<std::string> lib_files = std::vector<std::string>();
        if(!amino_lib_files.empty())
            for(std::vector<std::string>::iterator it = amino_lib_files.begin(); it != amino_lib_files.end(); it++)
                lib_files.push_back(*it);
        if(!glycam_lib_files.empty())
            for(std::vector<std::string>::iterator it = glycam_lib_files.begin(); it != glycam_lib_files.end(); it++)
                lib_files.push_back(*it);
        if(!other_lib_files.empty())
            for(std::vector<std::string>::iterator it = other_lib_files.begin(); it != other_lib_files.end(); it++)
                lib_files.push_back(*it);

        if(lib_files.size() != 0)
            lib_residues = GetAllResiduesFromMultipleLibFilesMap(lib_files);

        if(prep_files.size() != 0)
            prep_residues = GetAllResiduesFromMultiplePrepFilesMap(prep_files);


        std::vector<std::string> key_order = std::vector<std::string>();
        PdbFileSpace::PdbFile::PdbResidueAtomsMap residue_atoms_map = pdb_file->GetAllAtomsInOrder(key_order);
        
        input_file_ = pdb_file;
        int testPoly = this->input_file_->GetMasterCard()->GetNumRemarks();
        cout << testPoly << endl;


        for(std::vector<std::string>::iterator it = key_order.begin(); it != key_order.end(); it++)
        {
            std::string residue_key = *it;
            PdbFileSpace::PdbFile::PdbAtomCardVector* atoms = residue_atoms_map[residue_key];
            Residue* residue = new Residue();
            residue->SetAssembly(this);

            for(PdbFileSpace::PdbFile::PdbAtomCardVector::iterator it1 = atoms->begin(); it1 != atoms->end(); it1++)
            {
                PdbFileSpace::PdbAtomCard* atom = (*it1);
                std::string residue_name = atom->GetAtomResidueName();
                char chain_id = atom->GetAtomChainId();
                int sequence_number = atom->GetAtomResidueSequenceNumber();
                char insertion_code = atom->GetAtomInsertionCode();
                char alternate_location = atom->GetAtomAlternateLocation();
                std::stringstream ss;
                ss << residue_name << "_" << chain_id << "_" << sequence_number << "_" << insertion_code << "_"
                   << alternate_location << "_" << id_;
                std::string key = ss.str();
                residue->SetId(key);

                Atom* new_atom = new Atom();
                residue->SetName(residue_name);
                std::string atom_name = atom->GetAtomName();
                new_atom->SetName(atom_name);
                if(!lib_residues.empty() || !prep_residues.empty())
                {
                    if(lib_residues.find(residue_name) != lib_residues.end())
                    {
                        LibraryFileSpace::LibraryFileResidue* lib_residue = lib_residues[residue_name];
                        LibraryFileSpace::LibraryFileAtom* lib_atom = lib_residue->GetLibraryAtomByAtomName(atom_name);
                        if(lib_atom != NULL)
                        {
                            new_atom->MolecularDynamicAtom::SetAtomType(lib_atom->GetType());
                            new_atom->MolecularDynamicAtom::SetCharge(lib_atom->GetCharge());

                            if(parameter != NULL)
                            {
                                if(atom_type_map.find(new_atom->MolecularDynamicAtom::GetAtomType()) != atom_type_map.end())
                                {
                                    ParameterFileSpace::ParameterFileAtom* parameter_atom = atom_type_map[new_atom->GetAtomType()];
                                    new_atom->MolecularDynamicAtom::SetMass(parameter_atom->GetMass());
                                    new_atom->MolecularDynamicAtom::SetRadius(parameter_atom->GetRadius());
                                }
                                else
                                {
                                    new_atom->MolecularDynamicAtom::SetMass(gmml::dNotSet);
                                    new_atom->MolecularDynamicAtom::SetRadius(gmml::dNotSet);
                                }
                            }
                            else
                            {
                                new_atom->MolecularDynamicAtom::SetMass(gmml::dNotSet);
                                new_atom->MolecularDynamicAtom::SetRadius(gmml::dNotSet);
                            }
                        }
                        else
                        {
                            new_atom->MolecularDynamicAtom::SetAtomType("UNK");
                            new_atom->MolecularDynamicAtom::SetCharge(gmml::dNotSet);
                            new_atom->MolecularDynamicAtom::SetMass(gmml::dNotSet);
                            new_atom->MolecularDynamicAtom::SetRadius(gmml::dNotSet);
                        }
                    }
                    else if(prep_residues.find(residue_name) != prep_residues.end())
                    {
                        PrepFileSpace::PrepFileResidue* prep_residue = prep_residues[residue_name];
                        PrepFileSpace::PrepFileAtom* prep_atom = prep_residue->GetPrepAtomByName(atom_name);
                        if(prep_atom != NULL)
                        {
                            new_atom->MolecularDynamicAtom::SetAtomType(prep_atom->GetType());
                            new_atom->MolecularDynamicAtom::SetCharge(prep_atom->GetCharge());

                            if(parameter != NULL)
                            {
                                if(atom_type_map.find(new_atom->MolecularDynamicAtom::GetAtomType()) != atom_type_map.end())
                                {
                                    ParameterFileSpace::ParameterFileAtom* parameter_atom = atom_type_map[new_atom->GetAtomType()];
                                    new_atom->MolecularDynamicAtom::SetMass(parameter_atom->GetMass());
                                    new_atom->MolecularDynamicAtom::SetRadius(parameter_atom->GetRadius());
                                }
                                else
                                {
                                    new_atom->MolecularDynamicAtom::SetMass(gmml::dNotSet);
                                    new_atom->MolecularDynamicAtom::SetRadius(gmml::dNotSet);
                                }
                            }
                            else
                            {
                                new_atom->MolecularDynamicAtom::SetMass(gmml::dNotSet);
                                new_atom->MolecularDynamicAtom::SetRadius(gmml::dNotSet);
                            }
                        }
                        else
                        {
                            new_atom->MolecularDynamicAtom::SetAtomType("UNK");
                            new_atom->MolecularDynamicAtom::SetCharge(gmml::dNotSet);
                            new_atom->MolecularDynamicAtom::SetMass(gmml::dNotSet);
                            new_atom->MolecularDynamicAtom::SetRadius(gmml::dNotSet);
                        }
                    }
                }
                new_atom->SetResidue(residue);
                std::stringstream atom_key;
                atom_key << atom_name << "_" << atom->GetAtomSerialNumber() << "_" << key;
                new_atom->SetId(atom_key.str());
                PdbFileSpace::PdbModelSection* models = pdb_file->GetModels();
                PdbFileSpace::PdbModelSection::PdbModelCardMap model_maps = models->GetModels();
                if(model_maps.size() == 1)
                {
                    new_atom->AddCoordinate(new GeometryTopology::Coordinate(atom->GetAtomOrthogonalCoordinate()));
                    std::vector<std::string> card_index = gmml::Split(atom->GetAtomCardIndexInResidueSet(), "_");
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
                    for(PdbFileSpace::PdbModelSection::PdbModelCardMap::iterator it2 = model_maps.begin(); it2 != model_maps.end(); it2++)
                    {
                        PdbFileSpace::PdbModelCard* model = (*it2).second;
                        PdbFileSpace::PdbModelResidueSet* residue_set = model->GetModelResidueSet();
                        PdbFileSpace::PdbModelResidueSet::AtomCardVector atom_cards = residue_set->GetAtomCards();
                        std::vector<std::string> card_index = gmml::Split(atom->GetAtomCardIndexInResidueSet(), "_");
                        if(card_index.at(0).compare("ATOM") == 0)
                        {
                            PdbFileSpace::PdbAtomSection* atom_card = atom_cards.at(gmml::ConvertString<int>(card_index.at(1)));
                            PdbFileSpace::PdbAtomSection::PdbAtomCardOrderVector atom_vector = atom_card->GetOrderedAtomCards();
                            for(PdbFileSpace::PdbAtomSection::PdbAtomCardOrderVector::iterator it3 = atom_vector.begin(); it3 != atom_vector.end(); it3++)
                            {
                                PdbFileSpace::PdbAtomCard* matching_atom = *it3;
                                std::string matching_residue_name = matching_atom->GetAtomResidueName();
                                char matching_chain_id = matching_atom->GetAtomChainId();
                                int matching_sequence_number = matching_atom->GetAtomResidueSequenceNumber();
                                char matching_insertion_code = matching_atom->GetAtomInsertionCode();
                                char matching_alternate_location = matching_atom->GetAtomAlternateLocation();
                                std::stringstream sss;
                                sss << matching_residue_name << "_" << matching_chain_id << "_" << matching_sequence_number << "_"
                                    << matching_insertion_code << "_" << matching_alternate_location << "_" << id_;
                                std::string matching_key = sss.str();

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
                            PdbFileSpace::PdbModelResidueSet::HeterogenAtomCardVector heterogen_atom_cards = residue_set->GetHeterogenAtomCards();
                            PdbFileSpace::PdbHeterogenAtomSection* heterogen_atom_card = heterogen_atom_cards.at(gmml::ConvertString<int>(card_index.at(1)));
                            PdbFileSpace::PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector heterogen_atom_vector = heterogen_atom_card->GetOrderedHeterogenAtomCards();
                            for(PdbFileSpace::PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector::iterator it3 = heterogen_atom_vector.begin(); it3 != heterogen_atom_vector.end(); it3++)
                            {
                                PdbFileSpace::PdbAtomCard* matching_heterogen_atom = *it3;
                                std::string matching_heterogen_residue_name = matching_heterogen_atom->GetAtomResidueName();
                                char matching_heterogen_chain_id = matching_heterogen_atom->GetAtomChainId();
                                int matching_heterogen_sequence_number = matching_heterogen_atom->GetAtomResidueSequenceNumber();
                                char matching_heterogen_insertion_code = matching_heterogen_atom->GetAtomInsertionCode();
                                char matching_heterogen_alternate_location = matching_heterogen_atom->GetAtomAlternateLocation();
                                std::stringstream ssss;
                                ssss << matching_heterogen_residue_name << "_" << matching_heterogen_chain_id << "_" << matching_heterogen_sequence_number << "_"
                                     << matching_heterogen_insertion_code << "_" << matching_heterogen_alternate_location << "_" << id_;
                                std::string matching_heterogen_key = ssss.str();

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

/** ***************************************************************************
 *          Build from PDBQT files
 * ************************************************************************* **/
void Assembly::BuildAssemblyFromPdbqtFile(std::string pdbqt_file_path, std::string parameter_file)
{
    std::cout << "Building assembly from pdbqt file ..." << std::endl;
    std::cout << "Reading PDBQT file into PdbqtFileSpace::PdbqtFile structure." << std::endl;
    PdbqtFileSpace::PdbqtFile pdbqt_file;
    try
    {
        gmml::log(__LINE__, __FILE__, gmml::INF, "Reading PDBQT file into PdbqtFileSpace::PdbqtFile structure ...");
        pdbqt_file = PdbqtFileSpace::PdbqtFile(pdbqt_file_path);
    }
    catch(PdbqtFileSpace::PdbqtFileProcessingException &ex)
    {
        std::cout << "Generating PdbqtFileSpace::PdbqtFile structure from " << pdbqt_file_path << "failed." << std::endl;
    }
    this->BuildAssemblyFromPdbqtFile(&pdbqt_file, parameter_file);
}


void Assembly::BuildAssemblyFromPdbqtFile(PdbqtFileSpace::PdbqtFile *pdbqt_file, std::string parameter_file)
{
    std::cout << "Building assembly from pdbqt file ..." << std::endl;
    gmml::log(__LINE__, __FILE__, gmml::INF, "Building assembly from pdbqt file ...");
    try
    {
        this->ClearAssembly();
        ParameterFileSpace::ParameterFile* parameter = NULL;
        ParameterFileSpace::ParameterFile::AtomTypeMap atom_type_map = ParameterFileSpace::ParameterFile::AtomTypeMap();
        if(parameter_file.compare("") != 0)
        {
            parameter = new ParameterFileSpace::ParameterFile(parameter_file);
            atom_type_map = parameter->GetAtomTypes();
        }
        std::vector<std::string> key_order = std::vector<std::string>();
        PdbqtFileSpace::PdbqtFile::PdbqtResidueAtomsMap residue_atoms_map = pdbqt_file->GetAllAtomsInOrder(key_order);
        for(std::vector<std::string>::iterator it = key_order.begin(); it != key_order.end(); it++)
        {
            std::string residue_key = *it;
            PdbqtFileSpace::PdbqtFile::PdbqtAtomVector* atoms = residue_atoms_map[residue_key];
            Residue* residue = new Residue();
            residue->SetAssembly(this);

            for(PdbqtFileSpace::PdbqtFile::PdbqtAtomVector::iterator it1 = atoms->begin(); it1 != atoms->end(); it1++)
            {
                PdbqtFileSpace::PdbqtAtom* atom = (*it1);
                std::string residue_name = atom->GetAtomResidueName();
                char chain_id = atom->GetAtomChainId();
                int sequence_number = atom->GetAtomResidueSequenceNumber();
                char insertion_code = atom->GetAtomInsertionCode();
                char alternate_location = atom->GetAtomAlternateLocation();
                std::stringstream ss;
                ss << residue_name << "_" << chain_id << "_" << sequence_number << "_" << insertion_code << "_"
                   << alternate_location << "_" << id_;
                std::string key = ss.str();
                residue->SetId(key);

                Atom* new_atom = new Atom();
                residue->SetName(residue_name);
                std::string atom_name = atom->GetAtomName();
                new_atom->SetName(atom_name);
                new_atom->MolecularDynamicAtom::SetCharge(atom->GetAtomCharge());
                new_atom->MolecularDynamicAtom::SetAtomType(atom->GetAtomType());
                if(parameter != NULL)
                {
                    if(atom_type_map.find(new_atom->GetAtomType()) != atom_type_map.end())
                    {
                        ParameterFileSpace::ParameterFileAtom* parameter_atom = atom_type_map[new_atom->GetAtomType()];
                        new_atom->MolecularDynamicAtom::SetMass(parameter_atom->GetMass());
                        new_atom->MolecularDynamicAtom::SetRadius(parameter_atom->GetRadius());
                    }
                    else
                    {
                        new_atom->MolecularDynamicAtom::SetMass(gmml::dNotSet);
                        new_atom->MolecularDynamicAtom::SetRadius(gmml::dNotSet);
                    }
                }
                else
                {
                    new_atom->MolecularDynamicAtom::SetMass(gmml::dNotSet);
                    new_atom->MolecularDynamicAtom::SetRadius(gmml::dNotSet);
                }
                new_atom->SetResidue(residue);
                std::stringstream atom_key;
                atom_key << atom_name << "_" << atom->GetAtomSerialNumber() << "_" << key;
                new_atom->SetId(atom_key.str());
                PdbqtFileSpace::PdbqtModelCard* models = pdbqt_file->GetModels();
                PdbqtFileSpace::PdbqtModelCard::PdbqtModelMap model_maps = models->GetModels();
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
                    for(PdbqtFileSpace::PdbqtModelCard::PdbqtModelMap::iterator it2 = model_maps.begin(); it2 != model_maps.end(); it2++)
                    {
                        PdbqtFileSpace::PdbqtModel* model = (*it2).second;
                        PdbqtFileSpace::PdbqtModelResidueSet* residue_set = model->GetModelResidueSet();
                        PdbqtFileSpace::PdbqtAtomCard* atom_card = residue_set->GetAtoms();
                        PdbqtFileSpace::PdbqtAtomCard::PdbqtAtomMap atom_map = atom_card->GetAtoms();
                        PdbqtFileSpace::PdbqtAtom* matching_atom = atom_map[atom->GetAtomSerialNumber()];
                        std::string matching_residue_name = matching_atom->GetAtomResidueName();
                        char matching_chain_id = matching_atom->GetAtomChainId();
                        int matching_sequence_number = matching_atom->GetAtomResidueSequenceNumber();
                        char matching_insertion_code = matching_atom->GetAtomInsertionCode();
                        char matching_alternate_location = matching_atom->GetAtomAlternateLocation();
                        std::stringstream sss;
                        sss << matching_residue_name << "_" << matching_chain_id << "_" << matching_sequence_number << "_"
                            << matching_insertion_code << "_" << matching_alternate_location << "_" << id_;
                        std::string matching_key = sss.str();

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

/** ***************************************************************************
 *          Build from AMBER Prmtop files
 * ************************************************************************* **/

void Assembly::BuildAssemblyFromTopologyFile(std::string topology_file_path, std::string parameter_file)
{
    std::cout << "Building assembly from an AMBER parameter-topology file ..." << std::endl;
    std::cout << "Reading AMBER parameter-topology file into TopologyFileSpace::TopologyFile structure." << std::endl;
    gmml::log(__LINE__, __FILE__, gmml::INF, "Reading AMBER parameter-topology file into TopologyFileSpace::TopologyFile structure.");
    TopologyFileSpace::TopologyFile topology_file;
    try
    {
        topology_file = TopologyFileSpace::TopologyFile(topology_file_path);
    }
    catch(TopologyFileSpace::TopologyFileProcessingException &ex)
    {
        std::cout << "Generating TopologyFileSpace::TopologyFile structure from " << topology_file_path << "failed." << std::endl;
    }
    this->BuildAssemblyFromTopologyFile(&topology_file, parameter_file);
}

void Assembly::BuildAssemblyFromTopologyFile(TopologyFileSpace::TopologyFile *topology_file, std::string parameter_file)
{
    std::cout << "Building assembly from topology file ..." << std::endl;
    gmml::log(__LINE__, __FILE__, gmml::INF, "Building assembly from topology file ...");
    this->ClearAssembly();
    ParameterFileSpace::ParameterFile* parameter = NULL;
    ParameterFileSpace::ParameterFile::AtomTypeMap atom_type_map = ParameterFileSpace::ParameterFile::AtomTypeMap();
    if(parameter_file.compare("") != 0)
    {
        parameter = new ParameterFileSpace::ParameterFile(parameter_file);
        atom_type_map = parameter->GetAtomTypes();
    }
    name_ = topology_file->GetTitle();
    sequence_number_ = 1;
    int serial_number = 0;
    TopologyFileSpace::TopologyAssembly::TopologyResidueVector topology_residues = topology_file->GetAssembly()->GetResidues();
    for(TopologyFileSpace::TopologyAssembly::TopologyResidueVector::iterator it = topology_residues.begin(); it != topology_residues.end(); it++)
    {
        Residue* assembly_residue = new Residue();
        assembly_residue->SetAssembly(this);
        TopologyFileSpace::TopologyResidue* topology_residue = (*it);
        std::string residue_name = topology_residue->GetResidueName();
        assembly_residue->SetName(residue_name);
        std::stringstream id;
        id << residue_name << "_" << gmml::BLANK_SPACE << "_" << topology_residue->GetIndex()
           << "_" << gmml::BLANK_SPACE << "_" << gmml::BLANK_SPACE << "_" << id_;
        assembly_residue->SetId(id.str());

        TopologyFileSpace::TopologyResidue::TopologyAtomVector topology_atoms = topology_residue->GetAtoms();

        for(TopologyFileSpace::TopologyResidue::TopologyAtomVector::iterator it1 = topology_atoms.begin(); it1 != topology_atoms.end(); it1++)
        {
            serial_number++;
            Atom* assembly_atom = new Atom();
            std::string atom_name = (*it1)->GetAtomName();
            assembly_atom->SetName(atom_name);
            std::stringstream atom_id;
            atom_id << atom_name << "_" << serial_number << "_" << id.str();
            assembly_atom->SetId(atom_id.str());
            TopologyFileSpace::TopologyAtom* topology_atom = (*it1);
            assembly_atom->MolecularDynamicAtom::SetCharge(topology_atom->GetAtomCharge() / gmml::CHARGE_DIVIDER);
            assembly_atom->MolecularDynamicAtom::SetMass(topology_atom->GetAtomMass());
            assembly_atom->MolecularDynamicAtom::SetAtomType(topology_atom->GetType());
            if(parameter != NULL)
            {
                if(atom_type_map.find(assembly_atom->GetAtomType()) != atom_type_map.end())
                {
                    ParameterFileSpace::ParameterFileAtom* parameter_atom = atom_type_map[assembly_atom->GetAtomType()];
                    assembly_atom->MolecularDynamicAtom::SetRadius(parameter_atom->GetRadius());
                }
                else
                {
                    assembly_atom->MolecularDynamicAtom::SetRadius(gmml::dNotSet);
                }
            }
            else
            {
                assembly_atom->MolecularDynamicAtom::SetRadius(gmml::dNotSet);
            }
            assembly_atom->SetResidue(assembly_residue);
            assembly_atom->SetName(topology_atom->GetAtomName());

            assembly_residue->AddAtom(assembly_atom);
        }
        this->AddResidue(assembly_residue);

    }
}

/** ***************************************************************************
 *          Build from AMBER Lib/OFF files
 * ************************************************************************* **/

void Assembly::BuildAssemblyFromLibraryFile(std::string library_file_path, std::string parameter_file)
{
    std::cout << "Building assembly from AMBER library/off file ..." << std::endl;
    std::cout << "Reading AMBER library/off file into LibraryFileSpace::LibraryFile structure." << std::endl;
    gmml::log(__LINE__, __FILE__, gmml::INF, "Reading AMBER library/off file into LibraryFileSpace::LibraryFile structure ...");
    LibraryFileSpace::LibraryFile library_file;
    try
    {
        library_file = LibraryFileSpace::LibraryFile(library_file_path);
    }
    catch(LibraryFileSpace::LibraryFileProcessingException &ex)
    {
        std::cout << "Generating LibraryFileSpace::LibraryFile structure from " << library_file_path << "failed." << std::endl;
    }
    this->BuildAssemblyFromLibraryFile(&library_file, parameter_file);
}



void Assembly::BuildAssemblyFromLibraryFile(LibraryFileSpace::LibraryFile *library_file, std::string parameter_file)
{
    std::cout << "Building assembly from library file ..." << std::endl;
    gmml::log(__LINE__, __FILE__, gmml::INF, "Building assembly from library file ...");
    this->ClearAssembly();
    ParameterFileSpace::ParameterFile* parameter = NULL;
    ParameterFileSpace::ParameterFile::AtomTypeMap atom_type_map = ParameterFileSpace::ParameterFile::AtomTypeMap();
    if(parameter_file.compare("") != 0)
    {
        parameter = new ParameterFileSpace::ParameterFile(parameter_file);
        atom_type_map = parameter->GetAtomTypes();
    }
    sequence_number_ = 1;
    LibraryFileSpace::LibraryFile::ResidueMap library_residues = library_file->GetResidues();
    std::stringstream ss;

    int sequence_number = 0;
    int serial_number = 0;
    for(LibraryFileSpace::LibraryFile::ResidueMap::iterator it = library_residues.begin(); it != library_residues.end(); it++)
    {
        sequence_number++;
        Residue* assembly_residue = new Residue();
        assembly_residue->SetAssembly(this);
        std::string residue_name = (*it).first;
        assembly_residue->SetName(residue_name);
        std::stringstream id;
        id << residue_name << "_" << gmml::BLANK_SPACE << "_" << sequence_number << "_" << gmml::BLANK_SPACE << "_"
           << gmml::BLANK_SPACE << "_" << id_;
        assembly_residue->SetId(id.str());
        LibraryFileSpace::LibraryFileResidue* library_residue = (*it).second;
        int lib_res_tail_atom_index = library_residue->GetTailAtomIndex();
        int lib_res_head_atom_index = library_residue->GetHeadAtomIndex();
        std::string library_residue_name = library_residue->GetName();
        if(std::distance(library_residues.begin(), it) == (int)library_residues.size()-1)
            ss << library_residue_name;
        else
            ss << library_residue_name << "-";

        LibraryFileSpace::LibraryFileResidue::AtomMap library_atoms = library_residue->GetAtoms();

        for(LibraryFileSpace::LibraryFileResidue::AtomMap::iterator it1 = library_atoms.begin(); it1 != library_atoms.end(); it1++)
        {
            serial_number++;
            Atom* assembly_atom = new Atom();
            LibraryFileSpace::LibraryFileAtom* library_atom = (*it1).second;
            std::string atom_name = library_atom->GetName();
            assembly_atom->SetName(atom_name);
            std::stringstream atom_id;
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
                    ParameterFileSpace::ParameterFileAtom* parameter_atom = atom_type_map[assembly_atom->GetAtomType()];
                    assembly_atom->MolecularDynamicAtom::SetMass(parameter_atom->GetMass());
                    assembly_atom->MolecularDynamicAtom::SetRadius(parameter_atom->GetRadius());
                }
                else
                {
                    assembly_atom->MolecularDynamicAtom::SetMass(gmml::dNotSet);
                    assembly_atom->MolecularDynamicAtom::SetRadius(gmml::dNotSet);
                }
            }
            else
            {
                assembly_atom->MolecularDynamicAtom::SetMass(gmml::dNotSet);
                assembly_atom->MolecularDynamicAtom::SetRadius(gmml::dNotSet);
            }

            GeometryTopology::Coordinate* coordinate = new GeometryTopology::Coordinate(library_atom->GetCoordinate());
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

/** ***************************************************************************
 *          Build from AMBER Prmtop & Inpcrd files
 * ************************************************************************* **/

void Assembly::BuildAssemblyFromTopologyCoordinateFile(std::string topology_file_path, std::string coordinate_file_path, std::string parameter_file)
{
    std::cout << "Building assembly from AMBER parameter-topology and unknown-style coordinate files ..." << std::endl;
    std::cout << "Reading AMBER parameter-topology and unknown-style coordinate files into their file structure." << std::endl;
    gmml::log(__LINE__, __FILE__, gmml::INF, "Reading AMBER parameter-topology and unknown-style coordinate files into their file structure.");
    TopologyFileSpace::TopologyFile topology_file;
    CoordinateFileSpace::CoordinateFile coordinate_file;
    try
    {
        topology_file = TopologyFileSpace::TopologyFile(topology_file_path);
    }
    catch(TopologyFileSpace::TopologyFileProcessingException &ex)
    {
        std::cout << "Generating TopologyFileSpace::TopologyFile structure from " << topology_file_path << "failed." << std::endl;
    }
    try
    {
        coordinate_file = CoordinateFileSpace::CoordinateFile(coordinate_file_path);
    }
    catch(CoordinateFileSpace::CoordinateFileProcessingException &ex)
    {
        std::cout << "Generating CoordinateFileSpace::CoordinateFile structure from " << coordinate_file_path << "failed." << std::endl;
    }
    this->BuildAssemblyFromTopologyCoordinateFile(&topology_file, &coordinate_file, parameter_file);
}

void Assembly::BuildAssemblyFromTopologyCoordinateFile(TopologyFileSpace::TopologyFile *topology_file, CoordinateFileSpace::CoordinateFile *coordinate_file, std::string parameter_file)
{
    std::cout << "Building assembly from topology and coordinate files ..." << std::endl;
    gmml::log(__LINE__, __FILE__, gmml::INF, "Building assembly from topology and coordinate files ...");
    this->ClearAssembly();
    ParameterFileSpace::ParameterFile* parameter = NULL;
    ParameterFileSpace::ParameterFile::AtomTypeMap atom_type_map = ParameterFileSpace::ParameterFile::AtomTypeMap();
    if(parameter_file.compare("") != 0)
    {
        parameter = new ParameterFileSpace::ParameterFile(parameter_file);
        atom_type_map = parameter->GetAtomTypes();
    }

    name_ = topology_file->GetTitle();
    sequence_number_ = 1;
    int serial_number = 0;
    TopologyFileSpace::TopologyAssembly::TopologyResidueVector topology_residues = topology_file->GetAssembly()->GetResidues();
    for(TopologyFileSpace::TopologyAssembly::TopologyResidueVector::iterator it = topology_residues.begin(); it != topology_residues.end(); it++)
    {
        Residue* assembly_residue = new Residue();
        assembly_residue->SetAssembly(this);
        TopologyFileSpace::TopologyResidue* topology_residue = (*it);
        std::string residue_name = topology_residue->GetResidueName();
        assembly_residue->SetName(residue_name);
        std::stringstream id;
        id << residue_name << "_" << gmml::BLANK_SPACE << "_" << topology_residue->GetIndex() << "_" << gmml::BLANK_SPACE << "_"
           << gmml::BLANK_SPACE << "_" << id_;
        assembly_residue->SetId(id.str());

        TopologyFileSpace::TopologyResidue::TopologyAtomVector topology_atoms = topology_residue->GetAtoms();

        for(TopologyFileSpace::TopologyResidue::TopologyAtomVector::iterator it1 = topology_atoms.begin(); it1 != topology_atoms.end(); it1++)
        {
            serial_number++;
            Atom* assembly_atom = new Atom();
            std::string atom_name = (*it1)->GetAtomName();
            assembly_atom->SetName(atom_name);
            std::stringstream atom_id;
            atom_id << atom_name << "_" << serial_number << "_" << id.str();
            assembly_atom->SetId(atom_id.str());
            TopologyFileSpace::TopologyAtom* topology_atom = (*it1);

            assembly_atom->MolecularDynamicAtom::SetCharge(topology_atom->GetAtomCharge() / gmml::CHARGE_DIVIDER);
            assembly_atom->MolecularDynamicAtom::SetAtomType(topology_atom->GetType());
            assembly_atom->MolecularDynamicAtom::SetMass(topology_atom->GetAtomMass());
            if(parameter != NULL)
            {
                if(atom_type_map.find(assembly_atom->GetAtomType()) != atom_type_map.end())
                {
                    ParameterFileSpace::ParameterFileAtom* parameter_atom = atom_type_map[assembly_atom->GetAtomType()];
                    assembly_atom->MolecularDynamicAtom::SetRadius(parameter_atom->GetRadius());
                }
                else
                {
                    assembly_atom->MolecularDynamicAtom::SetRadius(gmml::dNotSet);
                }
            }
            else
            {
                assembly_atom->MolecularDynamicAtom::SetRadius(gmml::dNotSet);
            }

            int topology_atom_index = topology_atom->GetIndex();

            assembly_atom->SetResidue(assembly_residue);
            assembly_atom->SetName(topology_atom->GetAtomName());

            std::vector<GeometryTopology::Coordinate*> coord_file_coordinates = coordinate_file->GetCoordinates();
            assembly_atom->AddCoordinate(coord_file_coordinates.at(topology_atom_index-1));
            assembly_residue->AddAtom(assembly_atom);
        }
        residues_.push_back(assembly_residue);
    }
}

/** ***************************************************************************
 *          Build from AMBER Prep files
 * ************************************************************************* **/

void Assembly::BuildAssemblyFromPrepFile(std::string prep_file_path, std::string parameter_file)
{
    std::cout << "Building assembly from prep file ..." << std::endl;
    std::cout << "Reading Prep file into PrepFileSpace::PrepFile structure." << std::endl;
    gmml::log(__LINE__, __FILE__, gmml::INF, "Reading Prep file into PrepFileSpace::PrepFile structure ...");
    PrepFileSpace::PrepFile prep_file;
    try
    {
        prep_file = PrepFileSpace::PrepFile(prep_file_path);
    }
    catch(PrepFileSpace::PrepFileProcessingException &ex)
    {
        std::cout << "Generating PrepFileSpace::PrepFile structure from " << prep_file_path << "failed." << std::endl;
    }
    this->BuildAssemblyFromPrepFile(&prep_file, parameter_file);
}



void Assembly::BuildAssemblyFromPrepFile(PrepFileSpace::PrepFile *prep_file, std::string parameter_file)
{
    std::cout << "Building assembly from prep file ..." << std::endl;
    gmml::log(__LINE__, __FILE__, gmml::INF, "Building assembly from prep file ...");
    this->ClearAssembly();
    ParameterFileSpace::ParameterFile* parameter = NULL;
    ParameterFileSpace::ParameterFile::AtomTypeMap atom_type_map = ParameterFileSpace::ParameterFile::AtomTypeMap();
    if(parameter_file.compare("") != 0)
    {
        parameter = new ParameterFileSpace::ParameterFile(parameter_file);
        atom_type_map = parameter->GetAtomTypes();
    }
    sequence_number_ = 1;
    PrepFileSpace::PrepFile::ResidueMap prep_residues = prep_file->GetResidues();
    std::stringstream ss;

    int sequence_number = 0;
    int serial_number = 0;
    for(PrepFileSpace::PrepFile::ResidueMap::iterator it = prep_residues.begin(); it != prep_residues.end(); it++)
    {
        sequence_number++;
        CoordinateVector cartesian_coordinate_list = CoordinateVector();
        int head_atom_index = INFINITY;
        int tail_atom_index = -INFINITY;
        Atom* head_atom = new Atom();
        Atom* tail_atom = new Atom();

        Residue* assembly_residue = new Residue();
        assembly_residue->SetAssembly(this);
        PrepFileSpace::PrepFileResidue* prep_residue = (*it).second;
	std::string prep_residue_name = prep_residue->GetName();
        assembly_residue->SetName(prep_residue_name);
        std::stringstream id;
        id << prep_residue_name << "_" << gmml::BLANK_SPACE << "_" << sequence_number << "_" << gmml::BLANK_SPACE << "_"
           << gmml::BLANK_SPACE << "_" << id_;
        assembly_residue->SetId(id.str());
        if(std::distance(prep_residues.begin(), it) == (int)prep_residues.size()-1)
            ss << prep_residue_name;
        else
            ss << prep_residue_name << "-";
        PrepFileSpace::PrepFileResidue::PrepFileAtomVector prep_atoms = prep_residue->GetAtoms();
        PrepFileSpace::PrepFileResidue::PrepFileAtomVector parent_atoms = prep_residue->GetAtomsParentVector();

        for(PrepFileSpace::PrepFileResidue::PrepFileAtomVector::iterator it1 = prep_atoms.begin(); it1 != prep_atoms.end(); it1++)
        {
            serial_number++;
            Atom* assembly_atom = new Atom();
            PrepFileSpace::PrepFileAtom* prep_atom = (*it1);

            assembly_atom->SetResidue(assembly_residue);
            std::string atom_name = prep_atom->GetName();
            assembly_atom->SetName(atom_name);
            std::stringstream atom_id;
            atom_id << atom_name << "_" << serial_number << "_" << id.str();
            assembly_atom->SetId(atom_id.str());

            assembly_atom->MolecularDynamicAtom::SetAtomType(prep_atom->GetType());
            assembly_atom->MolecularDynamicAtom::SetCharge(prep_atom->GetCharge());
            if(parameter != NULL)
            {
                if(atom_type_map.find(assembly_atom->GetAtomType()) != atom_type_map.end())
                {
                    ParameterFileSpace::ParameterFileAtom* parameter_atom = atom_type_map[assembly_atom->GetAtomType()];
                    assembly_atom->MolecularDynamicAtom::SetMass(parameter_atom->GetMass());
                    assembly_atom->MolecularDynamicAtom::SetRadius(parameter_atom->GetRadius());
                }
                else
                {
                    assembly_atom->MolecularDynamicAtom::SetMass(gmml::dNotSet);
                    assembly_atom->MolecularDynamicAtom::SetRadius(gmml::dNotSet);
                }
            }
            else
            {
                assembly_atom->MolecularDynamicAtom::SetMass(gmml::dNotSet);
                assembly_atom->MolecularDynamicAtom::SetRadius(gmml::dNotSet);
            }

            if(prep_residue->GetCoordinateType() == PrepFileSpace::kINT)
            {
                std::vector<GeometryTopology::Coordinate*> coordinate_list = std::vector<GeometryTopology::Coordinate*>();
                int index = std::distance(prep_atoms.begin(), it1);
                if(index == 0)
                {

                }
                if(index == 1)
                {
                    int parent_index = parent_atoms.at(index)->GetIndex() - 1;
                    GeometryTopology::Coordinate* parent_coordinate = cartesian_coordinate_list.at(parent_index);
                    coordinate_list.push_back(parent_coordinate);
                }
                if(index == 2)
                {
                    int parent_index = parent_atoms.at(index)->GetIndex() - 1;
                    int grandparent_index = parent_atoms.at(parent_index)->GetIndex() - 1;
                    GeometryTopology::Coordinate* grandparent_coordinate = cartesian_coordinate_list.at(grandparent_index);
                    GeometryTopology::Coordinate* parent_coordinate = cartesian_coordinate_list.at(parent_index);
                    coordinate_list.push_back(grandparent_coordinate);
                    coordinate_list.push_back(parent_coordinate);
                }
                if(index > 2)
                {
                    int parent_index = parent_atoms.at(index)->GetIndex() - 1;
                    int grandparent_index = parent_atoms.at(parent_index)->GetIndex() - 1;
                    int great_grabdparent_index = parent_atoms.at(grandparent_index)->GetIndex() - 1;
                    GeometryTopology::Coordinate* great_grandparent_coordinate = cartesian_coordinate_list.at(great_grabdparent_index);
                    GeometryTopology::Coordinate* grandparent_coordinate = cartesian_coordinate_list.at(grandparent_index);
                    GeometryTopology::Coordinate* parent_coordinate = cartesian_coordinate_list.at(parent_index);
                    coordinate_list.push_back(great_grandparent_coordinate);
                    coordinate_list.push_back(grandparent_coordinate);
                    coordinate_list.push_back(parent_coordinate);
                }
                GeometryTopology::Coordinate* coordinate = gmml::ConvertInternalCoordinate2CartesianCoordinate(coordinate_list, prep_atom->GetBondLength(),
                                                                                             prep_atom->GetAngle(), prep_atom->GetDihedral());
                cartesian_coordinate_list.push_back(coordinate);

                assembly_atom->AddCoordinate(coordinate);
            }
            else if(prep_residue->GetCoordinateType() == PrepFileSpace::kXYZ)
            {
                assembly_atom->AddCoordinate(new GeometryTopology::Coordinate(prep_atom->GetBondLength(), prep_atom->GetAngle(), prep_atom->GetDihedral()));
            }
            if(prep_atom->GetTopologicalType() == gmml::kTopTypeM && prep_atom->GetType().compare(prep_residue->GetDummyAtomType()) != 0)
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
