#include <math.h>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <set>
#include <queue>
#include <stack>

#include "../../includes/MolecularModeling/assembly.hpp"
#include "../../includes/MolecularModeling/residue.hpp"
#include "../../includes/MolecularModeling/atom.hpp"
#include "../../includes/MolecularModeling/atomnode.hpp"
#include "../../includes/InputSet/CondensedSequenceSpace/condensedsequence.hpp"
#include "../../includes/InputSet/CondensedSequenceSpace/condensedsequenceresidue.hpp"
#include "../../includes/InputSet/CondensedSequenceSpace/condensedsequenceamberprepresidue.hpp"
#include "../../includes/InputSet/TopologyFileSpace/topologyfile.hpp"
#include "../../includes/InputSet/TopologyFileSpace/topologyassembly.hpp"
#include "../../includes/InputSet/TopologyFileSpace/topologyresidue.hpp"
#include "../../includes/InputSet/TopologyFileSpace/topologyatom.hpp"
#include "../../includes/InputSet/TopologyFileSpace/topologybond.hpp"
#include "../../includes/InputSet/TopologyFileSpace/topologybondtype.hpp"
#include "../../includes/InputSet/TopologyFileSpace/topologyangle.hpp"
#include "../../includes/InputSet/TopologyFileSpace/topologyangletype.hpp"
#include "../../includes/InputSet/TopologyFileSpace/topologydihedral.hpp"
#include "../../includes/InputSet/TopologyFileSpace/topologydihedraltype.hpp"
#include "../../includes/InputSet/TopologyFileSpace/topologyatompair.hpp"
#include "../../includes/InputSet/CoordinateFileSpace/coordinatefile.hpp"
#include "../../includes/ParameterSet/PrepFileSpace/prepfile.hpp"
#include "../../includes/ParameterSet/PrepFileSpace/prepfileresidue.hpp"
#include "../../includes/ParameterSet/PrepFileSpace/prepfileatom.hpp"
#include "../../includes/InputSet/PdbFileSpace/pdbfile.hpp"
#include "../../includes/InputSet/PdbFileSpace/pdbtitlecard.hpp"
#include "../../includes/InputSet/PdbFileSpace/pdbmodelcard.hpp"
#include "../../includes/InputSet/PdbFileSpace/pdbmodel.hpp"
#include "../../includes/InputSet/PdbFileSpace/pdbmodelresidueset.hpp"
#include "../../includes/InputSet/PdbFileSpace/pdbatomcard.hpp"
#include "../../includes/InputSet/PdbFileSpace/pdbheterogenatomcard.hpp"
#include "../../includes/InputSet/PdbFileSpace/pdbatom.hpp"
#include "../../includes/InputSet/PdbFileSpace/pdbconnectcard.hpp"
#include "../../includes/InputSet/PdbFileSpace/pdblinkcard.hpp"
#include "../../includes/InputSet/PdbFileSpace/pdblink.hpp"
#include "../../includes/InputSet/PdbFileSpace/pdblinkresidue.hpp"
#include "../../includes/InputSet/PdbFileSpace/pdbfileprocessingexception.hpp"
#include "../../includes/InputSet/PdbqtFileSpace/pdbqtfile.hpp"
#include "../../includes/InputSet/PdbqtFileSpace/pdbqtatom.hpp"
#include "../../includes/InputSet/PdbqtFileSpace/pdbqtmodel.hpp"
#include "../../includes/InputSet/PdbqtFileSpace/pdbqtmodelcard.hpp"
#include "../../includes/InputSet/PdbqtFileSpace/pdbqtatomcard.hpp"
#include "../../includes/InputSet/PdbqtFileSpace/pdbqtmodelresidueset.hpp"
#include "../../includes/ParameterSet/LibraryFileSpace/libraryfile.hpp"
#include "../../includes/ParameterSet/LibraryFileSpace/libraryfileatom.hpp"
#include "../../includes/ParameterSet/LibraryFileSpace/libraryfileresidue.hpp"
#include "../../includes/ParameterSet/ParameterFileSpace/parameterfile.hpp"
#include "../../includes/ParameterSet/ParameterFileSpace/parameterfilebond.hpp"
#include "../../includes/ParameterSet/ParameterFileSpace/parameterfileangle.hpp"
#include "../../includes/ParameterSet/ParameterFileSpace/parameterfiledihedral.hpp"
#include "../../includes/ParameterSet/ParameterFileSpace/parameterfiledihedralterm.hpp"
#include "../../includes/ParameterSet/ParameterFileSpace/parameterfileatom.hpp"
#include "../../includes/ParameterSet/PrepFileSpace/prepfile.hpp"
#include "../../includes/ParameterSet/PrepFileSpace/prepfileresidue.hpp"
#include "../../includes/ParameterSet/PrepFileSpace/prepfileatom.hpp"
#include "../../includes/utils.hpp"
#include "../../includes/common.hpp"
#include "../../includes/GeometryTopology/grid.hpp"
#include "../../includes/GeometryTopology/cell.hpp"

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
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////

Assembly::Assembly() : description_(""), model_index_(0), sequence_number_(1), id_("1")
{
    residues_ = ResidueVector();
    assemblies_ = AssemblyVector();
}

Assembly::Assembly(vector<string> file_paths, gmml::InputFileType type)
{
    source_file_type_ = type;
    description_ = "";
    model_index_ = 0;
    sequence_number_ = 1;
    id_ = "1";
    switch(type)
    {
        case gmml::PDB:
            source_file_ = file_paths.at(0);
            residues_ = ResidueVector();
            BuildAssemblyFromPdbFile(source_file_);
            assemblies_ = AssemblyVector();
            break;
        case gmml::PDBQT:
            source_file_ = file_paths.at(0);
            residues_ = ResidueVector();
            BuildAssemblyFromPdbqtFile(source_file_);
            assemblies_ = AssemblyVector();
            break;
        case gmml::TOP:
            source_file_ = file_paths.at(0);
            residues_ = ResidueVector();
            BuildAssemblyFromTopologyFile(source_file_);
            assemblies_ = AssemblyVector();
            break;
        case gmml::LIB:
            source_file_ = file_paths.at(0);
            residues_ = ResidueVector();
            BuildAssemblyFromLibraryFile(source_file_);
            assemblies_ = AssemblyVector();
            break;
        case gmml::PREP:
            source_file_ = file_paths.at(0);
            residues_ = ResidueVector();
            BuildAssemblyFromPrepFile(source_file_);
            assemblies_ = AssemblyVector();
            break;
        case gmml::TOP_CRD:
            source_file_ = file_paths.at(0)+";"+file_paths.at(1);
            residues_ = ResidueVector();
            BuildAssemblyFromTopologyCoordinateFile(file_paths.at(0), file_paths.at(1));
            assemblies_ = AssemblyVector();
            break;
        case gmml::MULTIPLE:
            break;
        case gmml::UNKNOWN:
            break;
    }
}

Assembly::Assembly(Assembly *assembly) : description_(""), model_index_(0), sequence_number_(1), id_("1")
{
    source_file_ = assembly->GetSourceFile();
    assemblies_ = AssemblyVector();
    AssemblyVector assemblies = assembly->GetAssemblies();
    for(AssemblyVector::iterator it = assemblies.begin(); it != assemblies.end(); it++)
        assemblies_.push_back(new Assembly(*it));

    residues_ = ResidueVector();
    ResidueVector residues = assembly->GetResidues();
    for(ResidueVector::iterator it = residues.begin(); it != residues.end(); it++)
        residues_.push_back(new Residue(*it));
}

Assembly::Assembly(vector<vector<string> > file_paths, vector<gmml::InputFileType> types)
{
    stringstream name;
    stringstream source_file;
    sequence_number_ = 1;
    id_ = "1";
    for(unsigned int i = 0; i < file_paths.size(); i++)
    {
        vector<string> file = file_paths.at(i);
        gmml::InputFileType input_type = types.at(i);
        Assembly* assembly = new Assembly(file, input_type);
        assembly->SetSequenceNumber(i + 1);
        stringstream ss;
        ss << id_ << "." << i + 1;
        assembly->SetId(ss.str());
        assemblies_.push_back(assembly);
        if(i < file_paths.size() - 1)
        {
            name << assembly->GetName() << "-";
            source_file << assembly->GetSourceFile() << "#";
        }
        else
        {
            name << assembly->GetName();
            source_file << assembly->GetSourceFile();
        }
    }
    source_file_type_ = gmml::MULTIPLE;
    name_ = name.str();
    model_index_ = 0;
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////

string Assembly::GetName()
{
    return name_;
}
Assembly::AssemblyVector Assembly::GetAssemblies()
{
    return assemblies_;
}
Assembly::ResidueVector Assembly::GetResidues()
{
    return residues_;
}
string Assembly::GetChemicalType()
{
    return chemical_type_;
}
int Assembly::GetSequenceNumber()
{
    return sequence_number_;
}
string Assembly::GetId()
{
    return id_;
}
string Assembly::GetDescription()
{
    return description_;
}
string Assembly::GetSourceFile()
{
    return source_file_;
}
gmml::InputFileType Assembly::GetSourceFileType()
{
    return source_file_type_;
}
int Assembly::GetModelIndex()
{
    return model_index_;
}
Assembly::AtomVector Assembly::GetAllAtomsOfAssembly()
{
    AtomVector all_atoms_of_assembly = AtomVector();
    AssemblyVector assemblies = this->GetAssemblies();
    for(AssemblyVector::iterator it = assemblies.begin(); it != assemblies.end(); it++)
    {
        Assembly* assembly = (*it);
        AtomVector atoms_of_assembly = assembly->GetAllAtomsOfAssembly();
        for(AtomVector::iterator it1 = atoms_of_assembly.begin(); it1 != atoms_of_assembly.end(); it1++)
        {
            Atom* atom = (*it1);
            all_atoms_of_assembly.push_back(atom);
        }
    }
    ResidueVector residues = this->GetResidues();
    for(ResidueVector::iterator it = residues.begin(); it != residues.end(); it++)
    {
        Residue* residue = (*it);
        AtomVector atoms = residue->GetAtoms();
        for(AtomVector::iterator it1 = atoms.begin(); it1 != atoms.end(); it1++)
        {
            Atom* atom = (*it1);
            all_atoms_of_assembly.push_back(atom);
        }
    }
    return all_atoms_of_assembly;
}
Assembly::AtomVector Assembly::GetAllAtomsOfAssemblyExceptProteinWaterResiduesAtoms()
{
    AtomVector all_atoms_of_assembly = AtomVector();
    AssemblyVector assemblies = this->GetAssemblies();
    for(AssemblyVector::iterator it = assemblies.begin(); it != assemblies.end(); it++)
    {
        Assembly* assembly = (*it);
        AtomVector atoms_of_assembly = assembly->GetAllAtomsOfAssembly();
        for(AtomVector::iterator it1 = atoms_of_assembly.begin(); it1 != atoms_of_assembly.end(); it1++)
        {
            Atom* atom = (*it1);
            all_atoms_of_assembly.push_back(atom);
        }
    }
    ResidueVector residues = this->GetResidues();
    for(ResidueVector::iterator it = residues.begin(); it != residues.end(); it++)
    {
        Residue* residue = (*it);
        if(residue->GetName().compare("HOH") != 0)
        {
            AtomVector atoms = residue->GetAtoms();
            for(AtomVector::iterator it1 = atoms.begin(); it1 != atoms.end(); it1++)
            {
                Atom* atom = (*it1);
                if(atom->GetDescription().find("Het;") != string::npos)
                    all_atoms_of_assembly.push_back(atom);
            }
        }
    }
    return all_atoms_of_assembly;
}
Assembly::ResidueVector Assembly::GetAllResiduesOfAssembly()
{
    ResidueVector all_residues_of_assembly = ResidueVector();
    AssemblyVector assemblies = this->GetAssemblies();
    for(AssemblyVector::iterator it = assemblies.begin(); it != assemblies.end(); it++)
    {
        Assembly* assembly = (*it);
        ResidueVector residues_of_assembly = assembly->GetAllResiduesOfAssembly();
        for(ResidueVector::iterator it1 = residues_of_assembly.begin(); it1 != residues_of_assembly.end(); it1++)
        {
            all_residues_of_assembly.push_back(*it1);
        }
    }
    ResidueVector residues = this->GetResidues();
    for(ResidueVector::iterator it = residues.begin(); it != residues.end(); it++)
    {
        all_residues_of_assembly.push_back(*it);
    }
    return all_residues_of_assembly;
}
Assembly::CoordinateVector Assembly::GetAllCoordinates()
{
    CoordinateVector coordinates = CoordinateVector();
    for(AssemblyVector::iterator it = this->assemblies_.begin(); it != this->assemblies_.end(); it++)
    {
        Assembly* assembly = (*it);
        CoordinateVector assembly_coordinate = assembly->GetAllCoordinates();
        if(assembly_coordinate.size() == 0)
        {
            cout << "Central data structure is not complete in order for generating this type of file: Missing coordinate(s)" << endl;
            gmml::log(__LINE__, __FILE__, gmml::ERR, "Central data structure is not complete in order for generating this type of file: Missing coordinate(s)");
            return CoordinateVector();
        }
        for(CoordinateVector::iterator it1 = assembly_coordinate.begin(); it1 != assembly_coordinate.end(); it1++)
        {
            coordinates.push_back(*it1);
        }
    }
    for(ResidueVector::iterator it = this->residues_.begin(); it != this->residues_.end(); it++)
    {
        Residue* residue = (*it);
        AtomVector residue_atoms = residue->GetAtoms();
        for(AtomVector::iterator it1 = residue_atoms.begin(); it1 != residue_atoms.end(); it1++)
        {
            Atom* atom = (*it1);
            if(atom->GetCoordinates().size() == 0)
            {
                cout << "Central data structure is not complete in order for generating this type of file: Missing coordinate(s)" << endl;
                gmml::log(__LINE__, __FILE__, gmml::ERR, "Central data structure is not complete in order for generating this type of file: Missing coordinate(s)");
                return CoordinateVector();
            }
            else
            {
                coordinates.push_back(atom->GetCoordinates()[model_index_]);
            }
        }
    }
    return coordinates;
}
Assembly::NoteVector Assembly::GetNotes()
{
    return notes_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void Assembly::SetName(string name)
{
    name_ = name;
}
void Assembly::SetAssemblies(AssemblyVector assemblies)
{
    assemblies_.clear();
    for(AssemblyVector::iterator it = assemblies.begin(); it != assemblies.end(); it++)
    {
        assemblies_.push_back(*it);
    }
}
void Assembly::AddAssembly(Assembly *assembly)
{
    stringstream ss;
    ss << this->name_ << "-" << assembly->GetName();
    this->name_ = ss.str();
    stringstream sss;
    sss << this->source_file_ << "#" << assembly->GetSourceFile();
    this->source_file_ = sss.str();
    source_file_type_ = gmml::MULTIPLE;
    assembly->SetSequenceNumber(assemblies_.size() + 1);
    stringstream ssss;
    ssss << this->id_ << "." << assemblies_.size() + 1;
    assembly->UpdateIds(ssss.str());
    assembly->SetId(ssss.str());
    this->assemblies_.push_back(assembly);
}
void Assembly::UpdateIds(string new_id)
{
    for(AssemblyVector::iterator  it = assemblies_.begin(); it != assemblies_.end(); it++)
    {
        Assembly* assembly = *it;
        stringstream ss;
        ss << new_id << "." << assembly->GetId().substr(assembly->GetId().find_first_of(this->id_) + this->id_.size() + 1);
        (*it)->UpdateIds(ss.str());
        (*it)->SetId(ss.str());
    }
    for(ResidueVector::iterator it = residues_.begin(); it != residues_.end(); it++)
    {
        vector<string> id_tokens = Split((*it)->GetId(), "_");
        stringstream ss;
        for(int i = 0; i < id_tokens.size() - 1; i++)
        {
            ss << id_tokens.at(i) << "_";
        }
        ss << new_id;
        (*it)->SetId(ss.str());

        AtomVector atoms = (*it)->GetAtoms();
        for(AtomVector::iterator it1 = atoms.begin(); it1 != atoms.end(); it1++)
        {
            vector<string> id_tokens = Split((*it1)->GetId(), "_");
            stringstream ss;
            for(int i = 0; i < id_tokens.size() - 1; i++)
            {
                ss << id_tokens.at(i) << "_";
            }
            ss << new_id;
            (*it1)->SetId(ss.str());
        }
    }
}
void Assembly::SetResidues(ResidueVector residues)
{
    residues_.clear();
    for(ResidueVector::iterator it = residues.begin(); it != residues.end(); it++)
    {
        residues_.push_back(*it);
    }
}
void Assembly::AddResidue(Residue *residue)
{
    residues_.push_back(residue);
}
void Assembly::SetChemicalType(string chemical_type)
{
    chemical_type_ = chemical_type;
}
void Assembly::SetSequenceNumber(int sequence_number)
{
    sequence_number_ = sequence_number;
}
void Assembly::SetId(string id)
{
    id_ = id;
}
void Assembly::SetDescription(string description)
{
    description_ = description;
}
void Assembly::SetSourceFile(string source_file)
{
    source_file_ = source_file;
}
void Assembly::SetSourceFileType(gmml::InputFileType source_file_type)
{
    source_file_type_ = source_file_type;
}
void Assembly::SetModelIndex(int model_index)
{
    model_index_ = model_index;
}
void Assembly::SetNotes(NoteVector notes)
{
    notes_.clear();
    for(NoteVector::iterator it = notes.begin(); it != notes.end(); it++)
    {
        notes_.push_back(*it);
    }
}
void Assembly::AddNote(Note *note)
{
    notes_.push_back(note);
}

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

void Assembly::AttachResidues(Residue *residue, Residue *parent_residue, int branch_index, string parameter_file)
{
    ///Translate all atoms of the attached residue to place them in proper position with respect to the tail atom of the parent residue/assembly
    this->SetAttachedResidueBond(residue, parent_residue, branch_index, parameter_file);

    ///Rotate all atoms of the attached residue to set the proper bond angle between the attached residue and the parent residue
    this->SetAttachedResidueAngle(residue, parent_residue, branch_index, parameter_file);

    ///Rotate all atoms of the attached residue to set the proper Phi, Psi and Omega torsion angles
    this->SetAttachedResidueTorsion(residue, parent_residue, branch_index);
}

void Assembly::RemoveHydrogenAtAttachedPosition(Residue *residue, int branch_index)
{
    Atom* oxygen = residue->GetTailAtoms().at(branch_index);
    if(oxygen != NULL)
    {
        int oxygen_index = 1;
        if(oxygen->GetName().size() > 1 && isdigit(oxygen->GetName().at(1)))
            oxygen_index = ConvertString<int>(ConvertT<char>(oxygen->GetName().at(1)));

        Atom* hydrogen = NULL;

        AtomVector oxygen_neighbors = oxygen->GetNode()->GetNodeNeighbors();
        for(AtomVector::iterator it = oxygen_neighbors.begin(); it != oxygen_neighbors.end(); it++)
        {
            Atom* neighbor = *it;
            if(neighbor->GetName().at(0) == 'H' &&
                    (neighbor->GetName().size() > 1 && isdigit(neighbor->GetName().at(1)) &&
                     ConvertString<int>(ConvertT<char>(neighbor->GetName().at(1))) == oxygen_index))
            {
                hydrogen = neighbor;
                break;
            }

        }
        residue->RemoveAtom(hydrogen);
    }
}

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

void Assembly::AdjustCharge(Residue *residue, Residue *parent_residue, int branch_index)
{
    if(residue->GetName().compare("SO3") == 0)
    {
      parent_residue->GetTailAtoms().at(branch_index)->MolecularDynamicAtom::SetCharge(
                  parent_residue->GetTailAtoms().at(branch_index)->MolecularDynamicAtom::GetCharge() + 0.031);
    }
    else if(residue->GetName().compare("MEX") == 0 || residue->GetName().compare("ACX") == 0)
    {
        Atom* oxygen = parent_residue->GetTailAtoms().at(branch_index);
        Atom* carbon = NULL;
        int oxygen_index = 1;
        if(oxygen->GetName().size() > 1 && isdigit(oxygen->GetName().at(1)))
            oxygen_index = ConvertString<int>(ConvertT<char>(oxygen->GetName().at(1)));

        AtomVector oxygen_neighbors = oxygen->GetNode()->GetNodeNeighbors();
        for(AtomVector::iterator it = oxygen_neighbors.begin(); it != oxygen_neighbors.end(); it++)
        {
            Atom* neighbor = *it;
            if(neighbor->GetName().at(0) == 'C' &&
                    (neighbor->GetName().size() > 1 && isdigit(neighbor->GetName().at(1)) &&
                     ConvertString<int>(ConvertT<char>(neighbor->GetName().at(1))) == oxygen_index))
            {
                carbon = neighbor;
                break;
            }
        }
        if(carbon != NULL)
        {
            if(residue->GetName().compare("MEX") == 0)
                carbon->MolecularDynamicAtom::SetCharge(carbon->MolecularDynamicAtom::GetCharge() - 0.039);
            if(residue->GetName().compare("ACX") == 0)
                carbon->MolecularDynamicAtom::SetCharge(carbon->MolecularDynamicAtom::GetCharge() + 0.008);
        }
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

PdbFile* Assembly::BuildPdbFileStructureFromAssembly(int link_card_direction, int connect_card_existance)
{
    cout << "Creating PDB file" << endl;
    gmml::log(__LINE__, __FILE__, gmml::INF, "Creating PDB file ...");
    PdbFile* pdb_file = new PdbFile();
    PdbTitleCard* title_card = new PdbTitleCard();
    title_card->SetTitle("Generated by GMML");
    // Set pdb_file title card
    pdb_file->SetTitle(title_card);

    PdbModelCard* model_card = new PdbModelCard();
    PdbModelCard::PdbModelMap models = PdbModelCard::PdbModelMap();
    PdbModel* model = new PdbModel();
    model->SetModelSerialNumber(1);
    PdbModelResidueSet* residue_set = new PdbModelResidueSet();
    int serial_number = 1;
    int sequence_number = 1;

    AssemblytoPdbSequenceNumberMap assembly_to_sequence_number_map = AssemblytoPdbSequenceNumberMap();
    AssemblytoPdbSerialNumberMap assembly_to_serial_number_map = AssemblytoPdbSerialNumberMap();
    ExtractPdbModelCardFromAssembly(residue_set, serial_number, sequence_number, model_index_, assembly_to_sequence_number_map,
                                    assembly_to_serial_number_map);

    PdbLinkCard* link_card = new PdbLinkCard();
    ExtractPdbLinkCardFromAssembly(link_card, model_index_, assembly_to_sequence_number_map, link_card_direction);
    link_card->SetRecordName("LINK");
    pdb_file->SetLinks(link_card);

    if(connect_card_existance == 1)
    {
        PdbConnectCard* connect_card = new PdbConnectCard();
        ExtractPdbConnectCardFromAssembly(connect_card, assembly_to_serial_number_map);
        pdb_file->SetConnectivities(connect_card);
    }
    model->SetModelResidueSet(residue_set);
    models[1] = model;
    model_card->SetModels(models);
    pdb_file->SetModels(model_card);

    cout << "PDB file created" << endl;
    gmml::log(__LINE__, __FILE__, gmml::INF, "PDB file created");
    return pdb_file;
}

PdbqtFile* Assembly::BuildPdbqtFileStructureFromAssembly()
{
    cout << "Creating PDBQT file" << endl;
    gmml::log(__LINE__, __FILE__, gmml::INF, "Creating PDBQT file ...");
    PdbqtFile* pdbqt_file = new PdbqtFile();

    PdbqtModelCard* model_card = new PdbqtModelCard();
    PdbqtModelCard::PdbqtModelMap models = PdbqtModelCard::PdbqtModelMap();
    PdbqtModel* model = new PdbqtModel();
    model->SetModelSerialNumber(1);
    PdbqtModelResidueSet* residue_set = new PdbqtModelResidueSet();
    int serial_number = 1;
    int sequence_number = 1;

    ExtractPdbqtModelCardFromAssembly(residue_set, serial_number, sequence_number, model_index_);

    model->SetModelResidueSet(residue_set);
    models[1] = model;
    model_card->SetModels(models);
    pdbqt_file->SetModels(model_card);
    cout << "PDBQT file created" << endl;
    gmml::log(__LINE__, __FILE__, gmml::INF, "PDBQT file created");

    return pdbqt_file;
}

void Assembly::ExtractPdbModelCardFromAssembly(PdbModelResidueSet* residue_set, int &serial_number, int &sequence_number, int model_number,
                                               AssemblytoPdbSequenceNumberMap& assembly_to_pdb_sequence_number_map,
                                               AssemblytoPdbSerialNumberMap& assembly_to_pdb_serial_number_map)
{
    for(AssemblyVector::iterator it = this->assemblies_.begin(); it != this->assemblies_.end(); it++)
    {
        Assembly* assembly = (*it);
        AssemblyVector assemblies = assembly->GetAssemblies();
        for(AssemblyVector::iterator it1 = assemblies.begin(); it1 != assemblies.end(); it1++)
        {
            ExtractPdbModelCardFromAssembly(residue_set, serial_number, sequence_number, model_number, assembly_to_pdb_sequence_number_map,
                                            assembly_to_pdb_serial_number_map);
        }
        PdbAtomCard* atom_card = new PdbAtomCard();
        PdbHeterogenAtomCard* het_atom_card = new PdbHeterogenAtomCard();
        PdbAtomCard::PdbAtomMap atom_map = PdbAtomCard::PdbAtomMap();
        PdbAtomCard::PdbAtomOrderVector atom_vector = PdbAtomCard::PdbAtomOrderVector();
        PdbHeterogenAtomCard::PdbHeterogenAtomMap het_atom_map = PdbHeterogenAtomCard::PdbHeterogenAtomMap();
        PdbHeterogenAtomCard::PdbHeterogenAtomOrderVector het_atom_vector = PdbHeterogenAtomCard::PdbHeterogenAtomOrderVector();
        ResidueVector residues = assembly->GetResidues();
        for(ResidueVector::iterator it1 = residues.begin(); it1 != residues.end(); it1++)
        {
            Residue* residue = (*it1);
            AtomVector atoms = residue->GetAtoms();
            for(AtomVector::iterator it2 = atoms.begin(); it2 != atoms.end(); it2++)
            {
                Atom* atom = (*it2);
                vector<string> dscr = Split(atom->GetDescription(), ";");
                //                stringstream ss;
                //                ss << setw(2) << fixed << setprecision(1) << atom->MolecularDynamicAtom::GetCharge();
                string residue_name = atom->GetResidue()->GetName();
                if(residue_name.compare("TIP3PBOX") == 0 || residue_name.compare("TIP5PBOX") == 0)
                    residue_name = "HOH";
                PdbAtom* pdb_atom = new PdbAtom(serial_number, atom->GetName(), ' ', residue_name, ' ', sequence_number, ' ',
                                                *((atom->GetCoordinates()).at(model_number)), dNotSet, dNotSet, atom->GetElementSymbol(), "");//ss.str());
                vector<string> atom_id_tokens = Split(atom->GetId(), "_");
                pdb_atom->SetAtomChainId(ConvertString<char>(atom_id_tokens.at(3)));
                pdb_atom->SetAtomInsertionCode(ConvertString<char>(atom_id_tokens.at(5)));
                pdb_atom->SetAtomAlternateLocation(ConvertString<char>(atom_id_tokens.at(6)));
                assembly_to_pdb_sequence_number_map[gmml::ConvertString<int>(atom_id_tokens.at(4))] = sequence_number;
                assembly_to_pdb_serial_number_map[ConvertString<int>(atom_id_tokens.at(1))] = serial_number;

                if(find(dscr.begin(), dscr.end(), "Atom") != dscr.end())
                {
                    atom_map[serial_number] = pdb_atom;
                    atom_vector.push_back(pdb_atom);
                    serial_number++;
                }
                else if(find(dscr.begin(), dscr.end(), "Het") != dscr.end())
                {
                    het_atom_map[serial_number] = pdb_atom;
                    het_atom_vector.push_back(pdb_atom);
                    serial_number++;
                }
                else
                {
                    atom_map[serial_number] = pdb_atom;
                    atom_vector.push_back(pdb_atom);
                    serial_number++;
                }
            }
            sequence_number++;
        }
        atom_card->SetAtoms(atom_map);
        atom_card->SetOrderedAtoms(atom_vector);
        het_atom_card->SetHeterogenAtoms(het_atom_map);
        het_atom_card->SetOrderedHeterogenAtoms(het_atom_vector);
        residue_set->AddAtom(atom_card);
        residue_set->AddHeterogenAtom(het_atom_card);
    }
    PdbAtomCard* atom_card = new PdbAtomCard();
    PdbHeterogenAtomCard* het_atom_card = new PdbHeterogenAtomCard();
    PdbAtomCard::PdbAtomMap atom_map = PdbAtomCard::PdbAtomMap();
    PdbAtomCard::PdbAtomOrderVector atom_vector = PdbAtomCard::PdbAtomOrderVector();
    PdbHeterogenAtomCard::PdbHeterogenAtomMap het_atom_map = PdbHeterogenAtomCard::PdbHeterogenAtomMap();
    PdbHeterogenAtomCard::PdbHeterogenAtomOrderVector het_atom_vector = PdbHeterogenAtomCard::PdbHeterogenAtomOrderVector();
    for(ResidueVector::iterator it1 = residues_.begin(); it1 != residues_.end(); it1++)
    {
        Residue* residue = (*it1);
        AtomVector atoms = residue->GetAtoms();
        for(AtomVector::iterator it2 = atoms.begin(); it2 != atoms.end(); it2++)
        {
            Atom* atom = (*it2);
            vector<string> dscr = Split(atom->GetDescription(), ";");
            //            stringstream ss;
            //            ss << setw(2) << fixed << setprecision(1) << atom->MolecularDynamicAtom::GetCharge();
            string residue_name = atom->GetResidue()->GetName();
            if(residue_name.compare("TIP3PBOX") == 0 || residue_name.compare("TIP5PBOX") == 0)
                residue_name = "HOH";
            PdbAtom* pdb_atom = new PdbAtom(serial_number, atom->GetName(), ' ', residue_name, ' ', sequence_number, ' ',
                                            *((atom->GetCoordinates()).at(model_number)), dNotSet, dNotSet, atom->GetElementSymbol(), "");//ss.str());
            vector<string> atom_id_tokens = Split(atom->GetId(), "_");
            pdb_atom->SetAtomChainId(ConvertString<char>(atom_id_tokens.at(3)));
            pdb_atom->SetAtomInsertionCode(ConvertString<char>(atom_id_tokens.at(5)));
            pdb_atom->SetAtomAlternateLocation(ConvertString<char>(atom_id_tokens.at(6)));
            assembly_to_pdb_sequence_number_map[gmml::ConvertString<int>(atom_id_tokens.at(4))] = sequence_number;
            assembly_to_pdb_serial_number_map[ConvertString<int>(atom_id_tokens.at(1))] = serial_number;

            if(find(dscr.begin(), dscr.end(), "Atom") != dscr.end())
            {

                atom_map[serial_number] = pdb_atom;
                atom_vector.push_back(pdb_atom);
                serial_number++;
            }
            else if(find(dscr.begin(), dscr.end(), "Het") != dscr.end())
            {
                het_atom_map[serial_number] = pdb_atom;
                het_atom_vector.push_back(pdb_atom);
                serial_number++;
            }
            else
            {
                atom_map[serial_number] = pdb_atom;
                atom_vector.push_back(pdb_atom);
                serial_number++;
            }
        }
        sequence_number++;
    }
    atom_card->SetAtoms(atom_map);
    atom_card->SetOrderedAtoms(atom_vector);
    het_atom_card->SetHeterogenAtoms(het_atom_map);
    het_atom_card->SetOrderedHeterogenAtoms(het_atom_vector);
    residue_set->AddAtom(atom_card);
    residue_set->AddHeterogenAtom(het_atom_card);
}

void Assembly::ExtractPdbqtModelCardFromAssembly(PdbqtModelResidueSet* residue_set, int &serial_number, int &sequence_number, int model_number)
{
    PdbqtAtomCard* atom_card = new PdbqtAtomCard();
    PdbqtAtomCard::PdbqtAtomMap atom_map = PdbqtAtomCard::PdbqtAtomMap();
    for(AssemblyVector::iterator it = this->assemblies_.begin(); it != this->assemblies_.end(); it++)
    {
        Assembly* assembly = (*it);
        AssemblyVector assemblies = assembly->GetAssemblies();
        for(AssemblyVector::iterator it1 = assemblies.begin(); it1 != assemblies.end(); it1++)
        {
            ExtractPdbqtModelCardFromAssembly(residue_set, serial_number, sequence_number, model_number);
        }
        ResidueVector residues = assembly->GetResidues();
        for(ResidueVector::iterator it1 = residues.begin(); it1 != residues.end(); it1++)
        {
            Residue* residue = (*it1);
            AtomVector atoms = residue->GetAtoms();
            for(AtomVector::iterator it2 = atoms.begin(); it2 != atoms.end(); it2++)
            {
                Atom* atom = (*it2);
                vector<string> dscr = Split(atom->GetDescription(), ";");
                if(find(dscr.begin(), dscr.end(), "Atom") != dscr.end())
                {
                    PdbqtAtom* pdb_atom = new PdbqtAtom(serial_number, atom->GetName(), ' ', atom->GetResidue()->GetName(), ' ', sequence_number, ' ',
                                                        *((atom->GetCoordinates()).at(model_number)), dNotSet, dNotSet, atom->MolecularDynamicAtom::GetCharge(),
                                                        atom->GetAtomType(), "ATOM");

                    atom_map[serial_number] = pdb_atom;
                    serial_number++;
                }
                else if(find(dscr.begin(), dscr.end(), "Het") != dscr.end())
                {
                    PdbqtAtom* pdb_atom = new PdbqtAtom(serial_number, atom->GetName(), ' ', atom->GetResidue()->GetName(), ' ', sequence_number, ' ',
                                                        *((atom->GetCoordinates()).at(model_number)), dNotSet, dNotSet, atom->MolecularDynamicAtom::GetCharge(),
                                                        atom->GetAtomType(), "HETATOM");
                    atom_map[serial_number] = pdb_atom;
                    serial_number++;
                }
                else
                {
                    PdbqtAtom* pdb_atom = new PdbqtAtom(serial_number, atom->GetName(), ' ', atom->GetResidue()->GetName(), ' ', sequence_number, ' ',
                                                        *((atom->GetCoordinates()).at(model_number)), dNotSet, dNotSet, atom->MolecularDynamicAtom::GetCharge(),
                                                        atom->GetAtomType(), "ATOM");

                    atom_map[serial_number] = pdb_atom;
                    serial_number++;
                }
            }
            sequence_number++;
        }
    }
    for(ResidueVector::iterator it1 = residues_.begin(); it1 != residues_.end(); it1++)
    {
        Residue* residue = (*it1);
        AtomVector atoms = residue->GetAtoms();
        for(AtomVector::iterator it2 = atoms.begin(); it2 != atoms.end(); it2++)
        {
            Atom* atom = (*it2);
            vector<string> dscr = Split(atom->GetDescription(), ";");

            if(find(dscr.begin(), dscr.end(), "Atom") != dscr.end())
            {
                PdbqtAtom* pdb_atom = new PdbqtAtom(serial_number, atom->GetName(), ' ', atom->GetResidue()->GetName(), ' ', sequence_number, ' ',
                                                    *((atom->GetCoordinates()).at(model_number)), dNotSet, dNotSet, atom->MolecularDynamicAtom::GetCharge(),
                                                    atom->GetAtomType(), "ATOM");
                atom_map[serial_number] = pdb_atom;
                serial_number++;
            }
            else if(find(dscr.begin(), dscr.end(), "Het") != dscr.end())
            {
                PdbqtAtom* pdb_atom = new PdbqtAtom(serial_number, atom->GetName(), ' ', atom->GetResidue()->GetName(), ' ', sequence_number, ' ',
                                                    *((atom->GetCoordinates()).at(model_number)), dNotSet, dNotSet, atom->MolecularDynamicAtom::GetCharge(),
                                                    atom->GetAtomType(), "HETATOM");
                atom_map[serial_number] = pdb_atom;
                serial_number++;
            }
            else
            {
                PdbqtAtom* pdb_atom = new PdbqtAtom(serial_number, atom->GetName(), ' ', atom->GetResidue()->GetName(), ' ', sequence_number, ' ',
                                                    *((atom->GetCoordinates()).at(model_number)), dNotSet, dNotSet, atom->MolecularDynamicAtom::GetCharge(),
                                                    atom->GetAtomType(), "ATOM");
                atom_map[serial_number] = pdb_atom;
                serial_number++;
            }
        }
        sequence_number++;
    }
    atom_card->SetAtoms(atom_map);
    residue_set->SetAtoms(atom_card);
}

void Assembly::ExtractPdbLinkCardFromAssembly(PdbLinkCard* link_card, int model_index, AssemblytoPdbSequenceNumberMap assembly_to_pdb_sequence_number,
                                              int link_card_direction)
{
    PdbLinkCard::LinkVector link_vector = PdbLinkCard::LinkVector();
    vector<string> visited_links = vector<string>();
    AtomVector all_atoms = GetAllAtomsOfAssembly();
    for(AtomVector::iterator it = all_atoms.begin(); it != all_atoms.end(); it++)
    {
        Atom* atom = (*it);
        Residue* residue = atom->GetResidue();
        AtomNode* node = atom->GetNode();
        if(node != NULL)
        {
            AtomVector neighbors = node->GetNodeNeighbors();

            if((link_card_direction == 1 && (atom->GetName().find("C") != string::npos) &&
                (find(visited_links.begin(), visited_links.end(), atom->GetId()) == visited_links.end())) ||
                    (link_card_direction == -1 && (atom->GetName().find("O") != string::npos) &&
                     (find(visited_links.begin(), visited_links.end(), atom->GetId()) == visited_links.end())))
            {
                for(AtomVector::iterator it1 = neighbors.begin(); it1 != neighbors.end(); it1++)
                {
                    Atom* neighbor = (*it1);
                    Residue* neighbor_residue = neighbor->GetResidue();
                    if(find(visited_links.begin(), visited_links.end(), neighbor->GetId()) == visited_links.end())
                    {
                        if(residue->GetId().compare(neighbor_residue->GetId()) != 0)
                        {
                            visited_links.push_back(atom->GetId());
                            visited_links.push_back(neighbor->GetId());
                            PdbLinkResidue* link_residue1 = new PdbLinkResidue();
                            vector<string> atom_id_tokens = Split(atom->GetId(), "_");
                            link_residue1->SetAtomName(atom->GetName());
                            link_residue1->SetResidueChainId(ConvertString<char>(atom_id_tokens.at(3)));
                            link_residue1->SetResidueSequenceNumber(assembly_to_pdb_sequence_number[ConvertString<int>(atom_id_tokens.at(4))]);
                            link_residue1->SetResidueInsertionCode(ConvertString<char>(atom_id_tokens.at(5)));
                            link_residue1->SetAlternateLocationIndicator(ConvertString<char>(atom_id_tokens.at(6)));
                            link_residue1->SetResidueName(residue->GetName());
                            PdbLinkResidue* link_residue2 = new PdbLinkResidue();
                            vector<string> neighbor_id_tokens = Split(neighbor->GetId(), "_");
                            link_residue2->SetAtomName(neighbor->GetName());
                            link_residue2->SetResidueChainId(ConvertString<char>(neighbor_id_tokens.at(3)));
                            link_residue2->SetResidueSequenceNumber(assembly_to_pdb_sequence_number[ConvertString<int>(neighbor_id_tokens.at(4))]);
                            link_residue2->SetResidueInsertionCode(ConvertString<char>(neighbor_id_tokens.at(5)));
                            link_residue2->SetAlternateLocationIndicator(ConvertString<char>(neighbor_id_tokens.at(6)));
                            link_residue2->SetResidueName(neighbor_residue->GetName());

                            PdbLink* pdb_link = new PdbLink();
                            pdb_link->AddResidue(link_residue1);
                            pdb_link->AddResidue(link_residue2);
                            double distance = atom->GetCoordinates().at(model_index)->Distance(*(neighbor->GetCoordinates().at(model_index)));
                            pdb_link->SetLinkLength(distance);

                            link_vector.push_back(pdb_link);
                        }
                    }
                }
            }
        }
    }
    link_card->SetResidueLinks(link_vector);
}

void Assembly::ExtractPdbConnectCardFromAssembly(PdbConnectCard *connect_card, AssemblytoPdbSerialNumberMap assembly_to_pdb_serial_number)
{
    PdbConnectCard::BondedAtomsSerialNumbersMap bonded_atoms_serial_number_map = PdbConnectCard::BondedAtomsSerialNumbersMap();
    AtomVector all_atoms = this->GetAllAtomsOfAssembly();
    for(AtomVector::iterator it = all_atoms.begin(); it != all_atoms.end(); it++)
    {
        Atom* atom = *it;
        AtomNode* node = atom->GetNode();
        if(node != NULL)
        {
            AtomVector neighbors = node->GetNodeNeighbors();
            vector<string> atom_id_tokens = Split(atom->GetId(), "_");
            int atom_serial_number = assembly_to_pdb_serial_number[ConvertString<int>(atom_id_tokens.at(1))];
            bonded_atoms_serial_number_map[atom_serial_number] = vector<int>();
            for(AtomVector::iterator it1 = neighbors.begin(); it1 != neighbors.end(); it1++)
            {
                Atom* neighbor = *it1;
                vector<string> neighbor_id_tokens = Split(neighbor->GetId(), "_");
                bonded_atoms_serial_number_map[atom_serial_number].push_back(assembly_to_pdb_serial_number[ConvertString<int>(neighbor_id_tokens.at(1))]);
            }
        }
    }
    connect_card->SetBondedAtomsSerialNumbers(bonded_atoms_serial_number_map);
}

PrepFile* Assembly::BuildPrepFileStructureFromAssembly(string parameter_file_path)
{
    cout << "Creating prep file ..." << endl;
    gmml::log(__LINE__, __FILE__, gmml::INF, "Creating prep file ...");
    PrepFile* prep_file = new PrepFile();
    ResidueVector assembly_residues = this->GetAllResiduesOfAssembly();
    PrepFile::ResidueMap prep_residues = PrepFile::ResidueMap();
    ParameterFile* parameter_file = new ParameterFile(parameter_file_path);
    ParameterFileSpace::ParameterFile::DihedralMap dihedrals = parameter_file->GetDihedrals();
    for(ResidueVector::iterator it = assembly_residues.begin(); it != assembly_residues.end(); it++)
    {
        Residue* assembly_residue = *it;
        vector<string> inserted_improper_dihedral_types = vector<string>();
        vector<vector<string> > inserted_improper_dihedrals = vector<vector<string> >();
        PrepFileResidue* prep_residue = new PrepFileResidue();
        PrepFileResidue::PrepFileAtomVector prep_atoms = PrepFileResidue::PrepFileAtomVector();
        prep_residue->SetTitle(assembly_residue->GetName());
        prep_residue->SetName(assembly_residue->GetName());
        prep_residue->SetGeometryType(PrepFileSpace::kGeometryCorrect);
        prep_residue->SetCoordinateType(PrepFileSpace::kINT);
        prep_residue->SetDummyAtomOmission(PrepFileSpace::kOmit);
        prep_residue->SetDummyAtomType("DU");
        prep_residue->SetDummyAtomPosition(PrepFileSpace::kPositionBeg);
        prep_residue->SetOutputFormat(PrepFileSpace::kFormatted);

        AtomVector assembly_atoms = assembly_residue->GetAtoms();
        CoordinateVector cartesian_coordinate_list = CoordinateVector();
        int atom_index = 1;
        PrepFileResidue::Loop loops = PrepFileResidue::Loop();
        vector<int> bond_index = vector<int>();
        for(int i = 0; i < DEFAULT_DUMMY_ATOMS; i ++)
        {
            PrepFileAtom* dummy_atom = new PrepFileAtom();
            dummy_atom->SetIndex(atom_index);
            dummy_atom->SetName("DUMM");
            dummy_atom->SetType("DU");
            dummy_atom->SetTopologicalType(gmml::kTopTypeM);
            dummy_atom->SetBondIndex(i);
            bond_index.push_back(i);
            dummy_atom->SetAngleIndex(i-1);
            dummy_atom->SetDihedral(i-2);
            if(i <= 0)
            {
                dummy_atom->SetBondLength(0.0);
                dummy_atom->SetAngle(0.0);
                dummy_atom->SetDihedral(0.0);
                cartesian_coordinate_list.push_back(new Coordinate(0.0, 0.0, 0.0)); //-5.0, -5.0, -5.0
            }
            else if(i <= 1)
            {
                dummy_atom->SetBondLength(dCutOff);
                dummy_atom->SetAngle(0.0);
                dummy_atom->SetDihedral(0.0);
                cartesian_coordinate_list.push_back(new Coordinate(dCutOff, 0.0, 0.0)); //1.522, 0.0, 0.0 ; -4.5, -5.0, -5.0
            }
            else
            {
                dummy_atom->SetBondLength(dCutOff);
                dummy_atom->SetAngle(90.0);
                dummy_atom->SetDihedral(0.0);
                cartesian_coordinate_list.push_back(new Coordinate(dCutOff, dCutOff, 0.0)); //2.43347, 1.34044, 0.0 ; -4.7, -4.2, -4.0
            }
            dummy_atom->SetCharge(0.0);
            dummy_atom->SetTopologicalType(kTopTypeM);
            atom_index++;
            prep_atoms.push_back(dummy_atom);
        }
        vector<TopologicalType> residue_topological_types = GetAllTopologicalTypesOfAtomsOfResidue(assembly_atoms, loops, bond_index);
        for(AtomVector::iterator it1 = assembly_atoms.begin(); it1 != assembly_atoms.end(); it1++)
        {
            CoordinateVector coordinate_list = CoordinateVector();

            Atom* assembly_atom = (*it1);
            PrepFileAtom* prep_atom = new PrepFileAtom();
            prep_atom->SetIndex(atom_index);

            // Set topological type
            // Set bond index, angle index, dihedral index
            // Set parent, grandparent and great_grandparent of the current atom
            // Call coordinate conversion function
            prep_atom->SetBondIndex(bond_index.at(atom_index - 1));
            prep_atom->SetAngleIndex(bond_index.at(bond_index.at(atom_index - 1) - 1));
            prep_atom->SetDihedralIndex(bond_index.at(bond_index.at(bond_index.at(atom_index - 1) - 1) - 1));
            prep_atom->SetCharge(assembly_atom->MolecularDynamicAtom::GetCharge());
            prep_atom->SetName(assembly_atom->GetName());
            prep_atom->SetTopologicalType(residue_topological_types.at(atom_index - DEFAULT_DUMMY_ATOMS - 1));
            prep_atom->SetType(assembly_atom->GetAtomType());
            int parent_index = prep_atom->GetBondIndex() - 1;
            int grandparent_index = prep_atom->GetAngleIndex() - 1;
            int great_grandparent_index = prep_atom->GetDihedralIndex() - 1;
            coordinate_list.push_back(cartesian_coordinate_list.at(great_grandparent_index));
            coordinate_list.push_back(cartesian_coordinate_list.at(grandparent_index));
            coordinate_list.push_back(cartesian_coordinate_list.at(parent_index));
            cartesian_coordinate_list.push_back(assembly_atom->GetCoordinates().at(model_index_));

            Coordinate* internal_coordinate = gmml::ConvertCartesianCoordinate2InternalCoordinate(assembly_atom->GetCoordinates().at(model_index_),
                                                                                                  coordinate_list);
            prep_atom->SetBondLength(internal_coordinate->GetX());
            prep_atom->SetAngle(internal_coordinate->GetY());
            prep_atom->SetDihedral(internal_coordinate->GetZ());

            atom_index++;
            prep_atoms.push_back(prep_atom);

            ExtractPrepImproperDihedralTypesFromAssembly(assembly_atom, inserted_improper_dihedral_types, dihedrals);
            ExtractPrepImproperDihedralsFromAssembly(assembly_atom, inserted_improper_dihedrals, inserted_improper_dihedral_types, dihedrals);
        }
        prep_residue->SetImproperDihedrals(inserted_improper_dihedrals);
        prep_residue->SetLoops(loops);
        prep_residue->SetAtoms(prep_atoms);
        prep_residue->SetCharge(prep_residue->CalculatePrepResidueCharge());
        prep_residues[assembly_residue->GetName()] = prep_residue;
    }
    //    prep_file->SetPath();
    prep_file->SetResidues(prep_residues);
    return prep_file;
}

void Assembly::ExtractPrepImproperDihedralTypesFromAssembly(Atom *assembly_atom, vector<string>
                                                            & inserted_improper_dihedral_types, ParameterFile::DihedralMap& dihedrals)
{
    AtomNode* atom_node = assembly_atom->GetNode();
    if(atom_node != NULL)
    {
        AtomVector neighbors = atom_node->GetNodeNeighbors();
        if(neighbors.size() == 3)
        {
            Atom* neighbor1 = neighbors.at(0);
            Atom* neighbor2 = neighbors.at(1);
            Atom* neighbor3 = neighbors.at(2);
            vector<vector<string> > all_improper_dihedrals_atom_type_permutations = CreateAllAtomTypePermutationsforImproperDihedralType(neighbor1->GetAtomType(), neighbor2->GetAtomType(),
                                                                                                                                         neighbor3->GetAtomType(), assembly_atom->GetAtomType());
            for(vector<vector<string> >::iterator it = all_improper_dihedrals_atom_type_permutations.begin(); it != all_improper_dihedrals_atom_type_permutations.end(); it++)
            {
                vector<string> improper_dihedral_permutation = (*it);
                if(dihedrals[improper_dihedral_permutation] != NULL)
                {
                    stringstream ss;
                    ss << improper_dihedral_permutation.at(0) << "_" << improper_dihedral_permutation.at(1) << "_" << improper_dihedral_permutation.at(2) << "_" << improper_dihedral_permutation.at(3);
                    if(find(inserted_improper_dihedral_types.begin(), inserted_improper_dihedral_types.end(), ss.str()) == inserted_improper_dihedral_types.end())
                    {
                        inserted_improper_dihedral_types.push_back(ss.str());
                        break;
                    }
                }
            }
        }
    }
}

void Assembly::ExtractPrepImproperDihedralsFromAssembly(Atom *assembly_atom, vector<vector<string> >& inserted_improper_dihedrals, vector<string> inserted_improper_dihedral_types,
                                                        ParameterFile::DihedralMap &dihedrals)
{
    AtomNode* atom_node = assembly_atom->GetNode();
    if(atom_node != NULL)
    {
        AtomVector neighbors = atom_node->GetNodeNeighbors();
        if(neighbors.size() == 3)
        {
            Atom* neighbor1 = neighbors.at(0);
            Atom* neighbor2 = neighbors.at(1);
            Atom* neighbor3 = neighbors.at(2);
            vector<vector<string> > all_improper_dihedrals_atom_type_permutations = CreateAllAtomTypePermutationsforImproperDihedralType(neighbor1->GetAtomType(), neighbor2->GetAtomType(),
                                                                                                                                         neighbor3->GetAtomType(), assembly_atom->GetAtomType());
            for(vector<vector<string> >::iterator it = all_improper_dihedrals_atom_type_permutations.begin(); it != all_improper_dihedrals_atom_type_permutations.end(); it++)
            {
                vector<string> improper_dihedral_permutation = (*it);
                stringstream sss;
                sss << improper_dihedral_permutation.at(0) << "_" << improper_dihedral_permutation.at(1) << "_" << improper_dihedral_permutation.at(2) << "_" << improper_dihedral_permutation.at(3);
                if(find(inserted_improper_dihedral_types.begin(), inserted_improper_dihedral_types.end(), sss.str()) != inserted_improper_dihedral_types.end())
                {
                    vector<string> dihedral1 = vector<string>();
                    vector<string> dihedral2 = vector<string>();
                    vector<string> dihedral3 = vector<string>();
                    vector<string> reverse_dihedral1 = vector<string>();
                    vector<string> reverse_dihedral2 = vector<string>();
                    vector<string> reverse_dihedral3 = vector<string>();

                    dihedral1.push_back(neighbor1->GetName());
                    dihedral1.push_back(neighbor2->GetName());
                    dihedral1.push_back(assembly_atom->GetName());
                    dihedral1.push_back(neighbor3->GetName());
                    reverse_dihedral1.push_back(neighbor3->GetName());
                    reverse_dihedral1.push_back(assembly_atom->GetName());
                    reverse_dihedral1.push_back(neighbor2->GetName());
                    reverse_dihedral1.push_back(neighbor1->GetName());

                    dihedral2.push_back(neighbor1->GetName());
                    dihedral2.push_back(assembly_atom->GetName());
                    dihedral2.push_back(neighbor3->GetName());
                    dihedral2.push_back(neighbor2->GetName());
                    reverse_dihedral2.push_back(neighbor2->GetName());
                    reverse_dihedral2.push_back(neighbor3->GetName());
                    reverse_dihedral2.push_back(assembly_atom->GetName());
                    reverse_dihedral2.push_back(neighbor1->GetName());

                    dihedral3.push_back(neighbor1->GetName());
                    dihedral3.push_back(neighbor3->GetName());
                    dihedral3.push_back(assembly_atom->GetName());
                    dihedral3.push_back(neighbor2->GetName());
                    reverse_dihedral3.push_back(neighbor2->GetName());
                    reverse_dihedral3.push_back(assembly_atom->GetName());
                    reverse_dihedral3.push_back(neighbor3->GetName());
                    reverse_dihedral3.push_back(neighbor1->GetName());

                    if(find(inserted_improper_dihedrals.begin(), inserted_improper_dihedrals.end(), dihedral1) == inserted_improper_dihedrals.end() &&
                            find(inserted_improper_dihedrals.begin(), inserted_improper_dihedrals.end(), dihedral2) == inserted_improper_dihedrals.end() &&
                            find(inserted_improper_dihedrals.begin(), inserted_improper_dihedrals.end(), dihedral3) == inserted_improper_dihedrals.end() &&
                            find(inserted_improper_dihedrals.begin(), inserted_improper_dihedrals.end(), reverse_dihedral1) == inserted_improper_dihedrals.end() &&
                            find(inserted_improper_dihedrals.begin(), inserted_improper_dihedrals.end(), reverse_dihedral2) == inserted_improper_dihedrals.end() &&
                            find(inserted_improper_dihedrals.begin(), inserted_improper_dihedrals.end(), reverse_dihedral3) == inserted_improper_dihedrals.end())
                    {
                        int permutation_index = distance(all_improper_dihedrals_atom_type_permutations.begin(), it);
                        ParameterFileDihedral* parameter_file_dihedral = NULL;
                        parameter_file_dihedral = dihedrals[improper_dihedral_permutation];

                        if(parameter_file_dihedral != NULL)
                        {
                            if(permutation_index % 6 == 0)
                            {
                                inserted_improper_dihedrals.push_back(dihedral1);
                            }
                            if(permutation_index % 6 == 2)
                            {
                                inserted_improper_dihedrals.push_back(dihedral2);
                            }
                            if(permutation_index % 6 == 4)
                            {
                                inserted_improper_dihedrals.push_back(dihedral3);
                            }
                            if(permutation_index % 6 == 1)
                            {
                                inserted_improper_dihedrals.push_back(reverse_dihedral1);
                            }
                            if(permutation_index % 6 == 3)
                            {
                                inserted_improper_dihedrals.push_back(reverse_dihedral2);
                            }
                            if(permutation_index % 6 == 5)
                            {
                                inserted_improper_dihedrals.push_back(reverse_dihedral3);
                            }
                            break;
                        }
                    }
                }
            }
        }
    }
}


vector<TopologicalType> Assembly::GetAllTopologicalTypesOfAtomsOfResidue(AtomVector assembly_atoms, PrepFileResidue::Loop& loops,
                                                                         vector<int> & bond_index, int dummy_atoms)
{
    vector<TopologicalType> topological_types = vector<TopologicalType>();
    vector<bool> visited = vector<bool>();
    vector<int> stack = vector<int>();
    for(AtomVector::iterator it = assembly_atoms.begin(); it != assembly_atoms.end(); it++)
    {
        topological_types.push_back(kTopTypeM);
        bond_index.push_back(0);
        visited.push_back(false);
    }

    for(AtomVector::iterator it = assembly_atoms.begin(); it != assembly_atoms.end(); it++)
    {
        Atom* assembly_atom = *it;
        int index = distance(assembly_atoms.begin(), it);

        // Finding all neighbors' indices of current atom in the assembly
        AtomNode* node = assembly_atom->GetNode();
        vector<int> neighbors_atom_index = vector<int>();
        if(node != NULL)
        {
            AtomVector atom_neighbors = node->GetNodeNeighbors();
            for(AtomVector::iterator it1 = atom_neighbors.begin(); it1 != atom_neighbors.end(); it1++)
            {
                Atom* atom_neighbor = *it1;
                for(AtomVector::iterator it2 = assembly_atoms.begin(); it2 != assembly_atoms.end(); it2++)
                {
                    Atom* atom = *it2;
                    if(atom->GetId().compare(atom_neighbor->GetId()) == 0)
                    {
                        int atom_index = distance(assembly_atoms.begin(), it2);
                        neighbors_atom_index.push_back(atom_index);
                        break;
                    }
                }
            }
        }

        int number_of_visited_neighbors = 0;
        for(int i = 0; i < neighbors_atom_index.size(); i++)
            if(visited.at(neighbors_atom_index.at(i)))
                number_of_visited_neighbors++;

        if(stack.empty())                   // Stack is empty
        {
            topological_types.at(index) = kTopTypeM;
            visited.at(index) = true;
            bond_index.at(index + dummy_atoms) = dummy_atoms;
            int neighbors = neighbors_atom_index.size();
            for(int i = 0; i < neighbors; i++)
                stack.push_back(index);
        }
        else                                // Stack is not empty
        {
            //AtomVector atom_neighbors = atom->GetNode()->GetNodeNeighbors();
            int top_stack_index = stack.size() - 1;
            int top_stack_atom_index = stack.at(top_stack_index);
            int neighbors = neighbors_atom_index.size();
            bool is_top_neighbor = false;
            if(find(neighbors_atom_index.begin(), neighbors_atom_index.end(), top_stack_atom_index) != neighbors_atom_index.end())
                is_top_neighbor = true;
            if(is_top_neighbor)             // Current atom is a neighbor of the atom on top of the stack
            {
                switch(neighbors)
                {
                    case 1:                 // Terminal atom
                        topological_types.at(index) = kTopTypeE;
                        visited.at(index) = true;
                        bond_index.at(index + dummy_atoms) = top_stack_atom_index + dummy_atoms + 1;
                        stack.pop_back();
                        break;
                    case 2:
                        switch(number_of_visited_neighbors)
                        {
                            case 1:
                                if(topological_types.at(top_stack_atom_index) == kTopTypeM)
                                    if(top_stack_index != 0 && stack.at(top_stack_index - 1) == top_stack_atom_index)
                                        topological_types.at(index) = kTopTypeS;
                                    else
                                        topological_types.at(index) = kTopTypeM;
                                else
                                    topological_types.at(index) = kTopTypeS;
                                stack.pop_back();
                                stack.push_back(index);
                                break;
                            case 2:
                                topological_types.at(index) = kTopTypeE;
                                stack.pop_back();
                                int loop_back_atom_index = 0;
                                for(int i = 0; i < neighbors_atom_index.size(); i++)
                                    if(neighbors_atom_index.at(i) != top_stack_atom_index && visited.at(neighbors_atom_index.at(i)) == true)
                                        loop_back_atom_index = neighbors_atom_index.at(i);
                                for(int i = 0; i < stack.size(); i++)
                                    if(stack.at(i) == loop_back_atom_index)
                                    {
                                        stack.erase(stack.begin() + i);
                                        break;
                                    }
                                loops[index + dummy_atoms + 1] = loop_back_atom_index + dummy_atoms + 1;
                                break;
                        }
                        visited.at(index) = true;
                        bond_index.at(index + dummy_atoms) = top_stack_atom_index + dummy_atoms + 1;

                        break;
                    case 3:
                        switch(number_of_visited_neighbors)
                        {
                            case 1:
                                if(topological_types.at(top_stack_atom_index) == kTopTypeM)
                                    if(top_stack_index != 0 && stack.at(top_stack_index - 1) == top_stack_atom_index)
                                        topological_types.at(index) = kTopTypeB;
                                    else
                                        topological_types.at(index) = kTopTypeM;
                                else
                                    topological_types.at(index) = kTopTypeB;
                                stack.pop_back();
                                stack.push_back(index);
                                stack.push_back(index);
                                break;
                            case 2:
                                topological_types.at(index) = kTopTypeS;
                                stack.pop_back();
                                int loop_back_atom_index = 0;
                                for(int i = 0; i < neighbors_atom_index.size(); i++)
                                    if(neighbors_atom_index.at(i) != top_stack_atom_index && visited.at(neighbors_atom_index.at(i)) == true)
                                        loop_back_atom_index = neighbors_atom_index.at(i);
                                for(int i = 0; i < stack.size(); i++)
                                    if(stack.at(i) == loop_back_atom_index)
                                    {
                                        stack.erase(stack.begin() + i);
                                        break;
                                    }
                                loops[index + dummy_atoms + 1] = loop_back_atom_index + dummy_atoms + 1;
                                stack.push_back(index);
                                break;
                        }
                        visited.at(index) = true;
                        bond_index.at(index + dummy_atoms) = top_stack_atom_index + dummy_atoms + 1;
                        break;
                    case 4:
                        switch(number_of_visited_neighbors)
                        {
                            case 1:
                                if(topological_types.at(top_stack_atom_index) == kTopTypeM)
                                    if(top_stack_index != 0 && stack.at(top_stack_index - 1) == top_stack_atom_index)
                                        topological_types.at(index) = kTopType3;
                                    else
                                        topological_types.at(index) = kTopTypeM;
                                else
                                    topological_types.at(index) = kTopType3;
                                stack.pop_back();
                                stack.push_back(index);
                                stack.push_back(index);
                                stack.push_back(index);
                                break;
                            case 2:
                                topological_types.at(index) = kTopTypeB;
                                stack.pop_back();
                                int loop_back_atom_index = 0;
                                for(int i = 0; i < neighbors_atom_index.size(); i++)
                                    if(neighbors_atom_index.at(i) != top_stack_atom_index && visited.at(neighbors_atom_index.at(i)) == true)
                                        loop_back_atom_index = neighbors_atom_index.at(i);
                                for(int i = 0; i < stack.size(); i++)
                                    if(stack.at(i) == loop_back_atom_index)
                                    {
                                        stack.erase(stack.begin() + i);
                                        break;
                                    }
                                loops[index + dummy_atoms + 1] = loop_back_atom_index + dummy_atoms + 1;
                                stack.push_back(index);
                                stack.push_back(index);
                                break;
                        }
                        visited.at(index) = true;
                        bond_index.at(index + dummy_atoms) = top_stack_atom_index + dummy_atoms + 1;
                        break;
                }
            }
            else                            // Current atom is not a neighbor of the atom on top of the stack
            {
                int stack_neighbor_index = -1;
                for(int i = stack.size() - 1; i >= 0; i--)
                {
                    if(find(neighbors_atom_index.begin(), neighbors_atom_index.end(), stack.at(i)) != neighbors_atom_index.end())
                    {
                        stack_neighbor_index = i;
                        break;
                    }
                }
                if(stack_neighbor_index == -1)
                {
                    //                    cout << "EMPTY" << endl;
                    topological_types.at(index) = kTopTypeM;
                    visited.at(index) = true;
                    bond_index.at(index + dummy_atoms) = dummy_atoms;
                    int neighbors = neighbors_atom_index.size();
                    for(int i = 0; i < neighbors; i++)
                        stack.push_back(index);
                }
                else
                {
                    int stack_neighbor_atom_index = stack.at(stack_neighbor_index);
                    switch(neighbors)
                    {
                        case 1:                 // Terminal atom
                            topological_types.at(index) = kTopTypeE;
                            visited.at(index) = true;
                            bond_index.at(index + dummy_atoms) = stack_neighbor_atom_index + dummy_atoms + 1;
                            stack.erase(stack.begin() + stack_neighbor_index);
                            break;
                        case 2:
                            switch(number_of_visited_neighbors)
                            {
                                case 1:
                                    if(topological_types.at(stack_neighbor_atom_index) == kTopTypeM)
                                        if(stack_neighbor_index != 0 && stack.at(stack_neighbor_index - 1) == stack_neighbor_atom_index)
                                            topological_types.at(index) = kTopTypeS;
                                        else
                                            topological_types.at(index) = kTopTypeM;
                                    else
                                        topological_types.at(index) = kTopTypeS;
                                    break;
                                    stack.erase(stack.begin() + stack_neighbor_index);
                                    stack.push_back(index);
                                case 2:
                                    topological_types.at(index) = kTopTypeE;
                                    stack.erase(stack.begin() + stack_neighbor_index);
                                    int loop_back_atom_index = 0;
                                    for(int i = 0; i < neighbors_atom_index.size(); i++)
                                        if(neighbors_atom_index.at(i) != top_stack_atom_index && visited.at(neighbors_atom_index.at(i)) == true)
                                            loop_back_atom_index = neighbors_atom_index.at(i);
                                    for(int i = 0; i < stack.size(); i++)
                                        if(stack.at(i) == loop_back_atom_index)
                                        {
                                            stack.erase(stack.begin() + i);
                                            break;
                                        }
                                    loops[index + dummy_atoms + 1] = loop_back_atom_index + dummy_atoms + 1;
                                    break;
                            }
                            visited.at(index) = true;
                            bond_index.at(index + dummy_atoms) = stack_neighbor_atom_index + dummy_atoms + 1;
                            break;
                        case 3:
                            switch(number_of_visited_neighbors)
                            {
                                case 1:
                                    if(topological_types.at(stack_neighbor_atom_index) == kTopTypeM)
                                        if(stack_neighbor_index != 0 && stack.at(stack_neighbor_index - 1) == stack_neighbor_atom_index)
                                            topological_types.at(index) = kTopTypeB;
                                        else
                                            topological_types.at(index) = kTopTypeM;
                                    else
                                        topological_types.at(index) = kTopTypeB;
                                    stack.erase(stack.begin() + stack_neighbor_index);
                                    stack.push_back(index);
                                    stack.push_back(index);
                                    break;
                                case 2:
                                    topological_types.at(index) = kTopTypeS;
                                    stack.erase(stack.begin() + stack_neighbor_index);
                                    int loop_back_atom_index = 0;
                                    for(int i = 0; i < neighbors_atom_index.size(); i++)
                                        if(neighbors_atom_index.at(i) != top_stack_atom_index && visited.at(neighbors_atom_index.at(i)) == true)
                                            loop_back_atom_index = neighbors_atom_index.at(i);
                                    for(int i = 0; i < stack.size(); i++)
                                        if(stack.at(i) == loop_back_atom_index)
                                        {
                                            stack.erase(stack.begin() + i);
                                            break;
                                        }
                                    loops[index + dummy_atoms + 1] = loop_back_atom_index + dummy_atoms + 1;
                                    stack.push_back(index);
                                    break;
                            }
                            visited.at(index) = true;
                            bond_index.at(index + dummy_atoms) = stack_neighbor_atom_index + dummy_atoms + 1;
                            break;
                        case 4:
                            switch(number_of_visited_neighbors)
                            {
                                case 1:
                                    if(topological_types.at(stack_neighbor_atom_index) == kTopTypeM)
                                        if(stack_neighbor_index != 0 && stack.at(stack_neighbor_index - 1) == stack_neighbor_atom_index)
                                            topological_types.at(index) = kTopType3;
                                        else
                                            topological_types.at(index) = kTopTypeM;
                                    else
                                        topological_types.at(index) = kTopType3;
                                    stack.erase(stack.begin() + stack_neighbor_index);
                                    stack.push_back(index);
                                    stack.push_back(index);
                                    stack.push_back(index);
                                    break;
                                case 2:
                                    topological_types.at(index) = kTopTypeB;
                                    stack.erase(stack.begin() + stack_neighbor_index);
                                    int loop_back_atom_index = 0;
                                    for(int i = 0; i < neighbors_atom_index.size(); i++)
                                        if(neighbors_atom_index.at(i) != top_stack_atom_index && visited.at(neighbors_atom_index.at(i)) == true)
                                            loop_back_atom_index = neighbors_atom_index.at(i);
                                    for(int i = 0; i < stack.size(); i++)
                                        if(stack.at(i) == loop_back_atom_index)
                                        {
                                            stack.erase(stack.begin() + i);
                                            break;
                                        }
                                    loops[index + dummy_atoms + 1] = loop_back_atom_index + dummy_atoms + 1;
                                    stack.push_back(index);
                                    stack.push_back(index);
                                    break;
                            }
                            visited.at(index) = true;
                            bond_index.at(index + dummy_atoms) = stack_neighbor_atom_index + dummy_atoms + 1;
                            break;
                    }
                }
            }
        }
    }
    return topological_types;
}

TopologyFile* Assembly::BuildTopologyFileStructureFromAssembly(string parameter_file_path, string ion_parameter_file_path)
{
    cout << "Creating topology file ..." << endl;
    gmml::log(__LINE__, __FILE__, gmml::INF, "Creating topology file ...");
    TopologyFile* topology_file = new TopologyFile();

    TopologyAssembly* topology_assembly = new TopologyAssembly();
    ResidueVector assembly_residues = this->GetAllResiduesOfAssembly();
    int number_of_excluded_atoms = 0;
    int residue_counter = 0;
    int atom_counter = 1;
    stringstream ss;
    int bond_type_counter = 0;
    vector<vector<string> > inserted_bond_types = vector<vector<string> >();
    vector<vector<string> > inserted_bonds = vector<vector<string> >();
    int angle_type_counter = 0;
    vector<vector<string> > inserted_angle_types = vector<vector<string> >();
    vector<vector<string> > inserted_angles = vector<vector<string> >();
    int dihedral_type_counter = 0;
    vector<string>  inserted_dihedral_types = vector<string>();
    vector<vector<string> > inserted_dihedrals = vector<vector<string> >();
    TopologyFile::TopologyAtomPairMap pairs = TopologyFile::TopologyAtomPairMap();
    int pair_count = 1;
    vector<string> inserted_pairs = vector<string>();
    vector<string> excluded_atom_list = vector<string>();
    ParameterFile* parameter_file = new ParameterFile(parameter_file_path);
    ParameterFileSpace::ParameterFile::BondMap bonds = parameter_file->GetBonds();
    ParameterFileSpace::ParameterFile::AngleMap angles = parameter_file->GetAngles();
    ParameterFileSpace::ParameterFile::DihedralMap dihedrals = parameter_file->GetDihedrals();
    ParameterFileSpace::ParameterFile::AtomTypeMap atom_types_map = parameter_file->GetAtomTypes();
    ParameterFile* ion_parameter_file = NULL;
    ParameterFileSpace::ParameterFile::AtomTypeMap ion_atom_types_map = ParameterFile::AtomTypeMap();
    if(ion_parameter_file_path.compare("") != 0)
    {
        ion_parameter_file = new ParameterFile(ion_parameter_file_path, gmml::IONICMOD);
        ion_atom_types_map = ion_parameter_file->GetAtomTypes();
    }
    for(ResidueVector::iterator it = assembly_residues.begin(); it != assembly_residues.end(); it++)
    {
        Residue* assembly_residue = *it;
        TopologyResidue* topology_residue = new TopologyResidue();
        residue_counter++;
        topology_residue->SetIndex(residue_counter);
        string residue_name = assembly_residue->GetName();
        if(residue_name.compare("TIP3") == 0 || residue_name.compare("TIP3PBOX") == 0)
            residue_name = "TP3";
        if(residue_name.compare("TIP5") == 0 || residue_name.compare("TIP5PBOX") == 0)
            residue_name = "TP5";
        topology_residue->SetResidueName(residue_name);
        if(distance(assembly_residues.begin(), it) == (int)assembly_residues.size()-1)
            ss << assembly_residue->GetName();
        else
            ss << assembly_residue->GetName() << "-";
        topology_residue->SetStartingAtomIndex(atom_counter);
        AtomVector assembly_atoms = assembly_residue->GetAtoms();
        PrepFileResidue::Loop loops = PrepFileResidue::Loop();
        vector<int> bond_index = vector<int>();
        int atom_index = 1;

        vector<TopologicalType> residue_topological_types = GetAllTopologicalTypesOfAtomsOfResidue(assembly_atoms, loops, bond_index, 0);

        for(AtomVector::iterator it1 = assembly_atoms.begin(); it1 != assembly_atoms.end(); it1++)
        {
            Atom* assembly_atom = (*it1);
            stringstream key1;
            key1 << assembly_atom->GetId();
            TopologyAtom* topology_atom = new TopologyAtom();
            topology_atom->SetAtomName(assembly_atom->GetName());
            topology_atom->SetAtomCharge(assembly_atom->MolecularDynamicAtom::GetCharge() * CHARGE_DIVIDER);
            topology_atom->SetAtomicNumber(iNotSet);
            topology_atom->SetAtomMass(assembly_atom->GetMass());
            topology_atom->SetResidueName(assembly_residue->GetName());
            topology_atom->SetType(assembly_atom->GetAtomType());
            if(atom_types_map.find(assembly_atom->GetAtomType()) == atom_types_map.end())
                cout << assembly_atom->GetAtomType() << " atom type is not found in parameter file" << endl;
            topology_atom->SetTreeChainClasification("0");
            topology_atom->SetRadii(dNotSet);
            topology_atom->SetScreen(dNotSet);
            topology_atom->SetIndex(atom_counter);
            topology_atom->SetTreeChainClasification(gmml::ConvertTopologicalType2String(residue_topological_types.at(atom_index - 1)));
            topology_residue->AddAtom(topology_atom);
            atom_counter++;
            atom_index++;

            ///Bond Types, Bonds
            AtomNode* atom_node = assembly_atom->GetNode();
            if(atom_node != NULL)
            {
                AtomVector neighbors = atom_node->GetNodeNeighbors();
                for(AtomVector::iterator it2 = neighbors.begin(); it2 != neighbors.end(); it2++)
                {
                    Atom* neighbor = (*it2);
                    stringstream key2;
                    key2 << neighbor->GetId();
                    ExtractTopologyBondTypesFromAssembly(inserted_bond_types, assembly_atom, neighbor, bonds, bond_type_counter, topology_file);
                    ExtractTopologyBondsFromAssembly(inserted_bonds, inserted_bond_types, assembly_atom, neighbor, topology_file);

                    ///Excluded Atoms
                    stringstream first_order_interaction;
                    stringstream reverse_first_order_interaction;
                    first_order_interaction << key1.str() << "-" << key2.str();
                    reverse_first_order_interaction << key2.str() << "-" << key1.str();
                    if(find(excluded_atom_list.begin(), excluded_atom_list.end(), first_order_interaction.str()) == excluded_atom_list.end() &&
                            find(excluded_atom_list.begin(), excluded_atom_list.end(), reverse_first_order_interaction.str()) == excluded_atom_list.end())
                    {
                        excluded_atom_list.push_back(first_order_interaction.str());
                        vector<string> key_tokens = Split(key2.str(),"_");
                        topology_atom->AddExcludedAtom(key_tokens.at(2) + "(" + key_tokens.at(4) + "):" + key_tokens.at(0) + "(" + key_tokens.at(1) +")");
                    }

                    ///Angle Types, Angle
                    AtomNode* neighbor_node = neighbor->GetNode();
                    AtomVector neighbors_of_neighbor = neighbor_node->GetNodeNeighbors();

                    for(AtomVector::iterator it3 = neighbors_of_neighbor.begin(); it3 != neighbors_of_neighbor.end(); it3++)
                    {
                        Atom* neighbor_of_neighbor = (*it3);
                        stringstream key3;
                        key3 << neighbor_of_neighbor->GetId();
                        if(key1.str().compare(key3.str()) != 0)
                        {
                            ExtractTopologyAngleTypesFromAssembly(assembly_atom, neighbor, neighbor_of_neighbor, inserted_angle_types, angle_type_counter,
                                                                  topology_file, angles);
                            ExtractTopologyAnglesFromAssembly(assembly_atom, neighbor, neighbor_of_neighbor, inserted_angles, inserted_angle_types, topology_file);

                            ///Excluded Atoms
                            stringstream second_order_interaction;
                            stringstream reverse_second_order_interaction;
                            second_order_interaction << key1.str() << "-" << key3.str();
                            reverse_second_order_interaction << key3.str() << "-" << key1.str();
                            if(find(excluded_atom_list.begin(), excluded_atom_list.end(), second_order_interaction.str()) == excluded_atom_list.end() &&
                                    find(excluded_atom_list.begin(), excluded_atom_list.end(), reverse_second_order_interaction.str()) == excluded_atom_list.end())
                            {
                                excluded_atom_list.push_back(second_order_interaction.str());
                                vector<string> key_tokens = Split(key3.str(),"_");
                                topology_atom->AddExcludedAtom(key_tokens.at(2) + "(" + key_tokens.at(4) + "):" + key_tokens.at(0) + "(" + key_tokens.at(1) +")");
                            }

                            //Dihedral Types, Dihedrals
                            AtomNode* neighbor_of_neighbor_node = neighbor_of_neighbor->GetNode();
                            AtomVector neighbors_of_neighbor_neighbor = neighbor_of_neighbor_node->GetNodeNeighbors();
                            for(AtomVector::iterator it4 =  neighbors_of_neighbor_neighbor.begin(); it4 != neighbors_of_neighbor_neighbor.end(); it4++)
                            {
                                Atom* neighbor_of_neighbor_of_neighbor = (*it4);
                                stringstream key4;
                                key4 << neighbor_of_neighbor_of_neighbor->GetId();
                                if(key2.str().compare(key4.str()) != 0)
                                {
                                    ExtractTopologyDihedralTypesFromAssembly(assembly_atom, neighbor, neighbor_of_neighbor, neighbor_of_neighbor_of_neighbor,
                                                                             inserted_dihedral_types, dihedral_type_counter, topology_file, dihedrals);
                                    ExtractTopologyDihedralsFromAssembly(assembly_atom, neighbor, neighbor_of_neighbor, neighbor_of_neighbor_of_neighbor,
                                                                         inserted_dihedrals, inserted_dihedral_types, dihedrals, topology_file);

                                    ///Excluded Atoms
                                    stringstream third_order_interaction;
                                    stringstream reverse_third_order_interaction;
                                    third_order_interaction << key1.str() << "-" << key4.str();
                                    reverse_third_order_interaction << key4.str() << "-" << key1.str();
                                    if(find(excluded_atom_list.begin(), excluded_atom_list.end(), third_order_interaction.str()) == excluded_atom_list.end() &&
                                            find(excluded_atom_list.begin(), excluded_atom_list.end(), reverse_third_order_interaction.str()) == excluded_atom_list.end())
                                    {
                                        excluded_atom_list.push_back(third_order_interaction.str());
                                        vector<string> key_tokens = Split(key4.str(),"_");
                                        topology_atom->AddExcludedAtom(key_tokens.at(2) + "(" + key_tokens.at(4) + "):" + key_tokens.at(0) + "(" + key_tokens.at(1) +")");
                                    }
                                }
                            }
                        }
                    }
                }
            }
            topology_atom->GetExcludedAtoms().size() == 0 ? number_of_excluded_atoms++ :
                                                            number_of_excluded_atoms += topology_atom->GetExcludedAtoms().size();
        }
        topology_assembly->AddResidue(topology_residue);
    }

    AtomVector all_atoms = this->GetAllAtomsOfAssembly();
    for(AtomVector::iterator it = all_atoms.begin(); it != all_atoms.end(); it++)
    {
        Atom* assembly_atom = *it;
        ///Pairs
        for(AtomVector::iterator it2 = all_atoms.begin(); it2 != all_atoms.end(); it2++)
        {
            Atom* pair_assembly_atom = (*it2);
            string atom_type1 = assembly_atom->GetAtomType();
            string atom_type2 = pair_assembly_atom->GetAtomType();
            vector<string> pair_vector = vector<string>();
            pair_vector.push_back(atom_type1);
            pair_vector.push_back(atom_type2);
            stringstream sss;
            sss << atom_type1 << "-" << atom_type2;
            stringstream reverse_sss;
            reverse_sss << atom_type2 << "-" << atom_type1;
            if(find(inserted_pairs.begin(), inserted_pairs.end(), sss.str()) == inserted_pairs.end() &&
                    find(inserted_pairs.begin(), inserted_pairs.end(), reverse_sss.str()) == inserted_pairs.end())
            {
                TopologyAtomPair* topology_atom_pair = new TopologyAtomPair();
                if(atom_types_map.find(atom_type1) != atom_types_map.end() && atom_types_map.find(atom_type2) != atom_types_map.end())
                {
                    ParameterFileAtom* parameter_atom1 = atom_types_map[atom_type1];
                    ParameterFileAtom* parameter_atom2 = atom_types_map[atom_type2];
                    double epsilon = sqrt(parameter_atom1->GetWellDepth() * parameter_atom2->GetWellDepth());
                    double sigma = pow(parameter_atom1->GetRadius() + parameter_atom2->GetRadius(), 6);
                    double coefficient_a = epsilon * sigma * sigma;
                    double coefficient_b = 2 * epsilon * sigma;
                    topology_atom_pair->SetCoefficientA(coefficient_a);
                    topology_atom_pair->SetCoefficientB(coefficient_b);
                    topology_atom_pair->SetPairType(sss.str());
                    topology_atom_pair->SetIndex(pair_count);
                    pair_count++;
                    inserted_pairs.push_back(sss.str());
                    pairs[sss.str()] = topology_atom_pair;
                }
                else if(!ion_atom_types_map.empty() && ion_atom_types_map.find(atom_type1) != ion_atom_types_map.end() &&
                        ion_atom_types_map.find(atom_type2) != ion_atom_types_map.end())
                {
                    ParameterFileAtom* parameter_atom1 = ion_atom_types_map[atom_type1];
                    ParameterFileAtom* parameter_atom2 = ion_atom_types_map[atom_type2];
                    double epsilon = sqrt(parameter_atom1->GetWellDepth() * parameter_atom2->GetWellDepth());
                    double sigma = pow(parameter_atom1->GetRadius() + parameter_atom2->GetRadius(), 6);
                    double coefficient_a = epsilon * sigma * sigma;
                    double coefficient_b = 2 * epsilon * sigma;
                    topology_atom_pair->SetCoefficientA(coefficient_a);
                    topology_atom_pair->SetCoefficientB(coefficient_b);
                    topology_atom_pair->SetPairType(sss.str());
                    topology_atom_pair->SetIndex(pair_count);
                    pair_count++;
                    inserted_pairs.push_back(sss.str());
                    pairs[sss.str()] = topology_atom_pair;
                }
                else if(atom_types_map.find(atom_type1) != atom_types_map.end() && !ion_atom_types_map.empty() &&
                        ion_atom_types_map.find(atom_type2) != ion_atom_types_map.end())
                {
                    ParameterFileAtom* parameter_atom1 = atom_types_map[atom_type1];
                    ParameterFileAtom* parameter_atom2 = ion_atom_types_map[atom_type2];
                    double epsilon = sqrt(parameter_atom1->GetWellDepth() * parameter_atom2->GetWellDepth());
                    double sigma = pow(parameter_atom1->GetRadius() + parameter_atom2->GetRadius(), 6);
                    double coefficient_a = epsilon * sigma * sigma;
                    double coefficient_b = 2 * epsilon * sigma;
                    topology_atom_pair->SetCoefficientA(coefficient_a);
                    topology_atom_pair->SetCoefficientB(coefficient_b);
                    topology_atom_pair->SetPairType(sss.str());
                    topology_atom_pair->SetIndex(pair_count);
                    pair_count++;
                    inserted_pairs.push_back(sss.str());
                    pairs[sss.str()] = topology_atom_pair;
                }
                else if(!ion_atom_types_map.empty() && ion_atom_types_map.find(atom_type1) != ion_atom_types_map.end() &&
                        atom_types_map.find(atom_type2) != atom_types_map.end())
                {
                    ParameterFileAtom* parameter_atom1 = ion_atom_types_map[atom_type1];
                    ParameterFileAtom* parameter_atom2 = atom_types_map[atom_type2];
                    double epsilon = sqrt(parameter_atom1->GetWellDepth() * parameter_atom2->GetWellDepth());
                    double sigma = pow(parameter_atom1->GetRadius() + parameter_atom2->GetRadius(), 6);
                    double coefficient_a = epsilon * sigma * sigma;
                    double coefficient_b = 2 * epsilon * sigma;
                    topology_atom_pair->SetCoefficientA(coefficient_a);
                    topology_atom_pair->SetCoefficientB(coefficient_b);
                    topology_atom_pair->SetPairType(sss.str());
                    topology_atom_pair->SetIndex(pair_count);
                    pair_count++;
                    inserted_pairs.push_back(sss.str());
                    pairs[sss.str()] = topology_atom_pair;
                }
            }
        }
    }

    topology_assembly->SetAssemblyName(ss.str());
    topology_file->SetAtomPairs(pairs);
    topology_file->SetAssembly(topology_assembly);

    // Set headers
    topology_file->SetNumberOfAtoms(this->CountNumberOfAtoms());
    topology_file->SetNumberOfTypes(this->CountNumberOfAtomTypes());
    topology_file->SetNumberOfBondsIncludingHydrogen(this->CountNumberOfBondsIncludingHydrogen(parameter_file_path));
    topology_file->SetNumberOfBondsExcludingHydrogen(this->CountNumberOfBondsExcludingHydrogen(parameter_file_path));
    topology_file->SetNumberOfAnglesIncludingHydrogen(this->CountNumberOfAnglesIncludingHydrogen(parameter_file_path));
    topology_file->SetNumberOfAnglesExcludingHydrogen(this->CountNumberOfAnglesExcludingHydrogen(parameter_file_path));
    topology_file->SetNumberOfDihedralsIncludingHydrogen(this->CountNumberOfDihedralsIncludingHydrogen(parameter_file_path));
    topology_file->SetNumberOfDihedralsExcludingHydrogen(this->CountNumberOfDihedralsExcludingHydrogen(parameter_file_path));
    topology_file->SetNumberOfExcludedAtoms(number_of_excluded_atoms);
    topology_file->SetNumberOfResidues(this->CountNumberOfResidues());
    topology_file->SetTotalNumberOfBonds(this->CountNumberOfBondsExcludingHydrogen(parameter_file_path));
    topology_file->SetTotalNumberOfAngles(this->CountNumberOfAnglesExcludingHydrogen(parameter_file_path));
    topology_file->SetTotalNumberOfDihedrals(this->CountNumberOfDihedralsExcludingHydrogen(parameter_file_path));
    topology_file->SetNumberOfBondTypes(this->CountNumberOfBondTypes(parameter_file_path));
    topology_file->SetNumberOfAngleTypes(this->CountNumberOfAngleTypes(parameter_file_path));
    topology_file->SetNumberOfDihedralTypes(this->CountNumberOfDihedralTypes(parameter_file_path));
    topology_file->SetNumberOfAtomsInLargestResidue(this->CountMaxNumberOfAtomsInLargestResidue());

    return topology_file;
}

void Assembly::ExtractTopologyBondTypesFromAssembly(vector<vector<string> > &inserted_bond_types, Atom* assembly_atom, Atom* neighbor, ParameterFileSpace::ParameterFile::BondMap &bonds,
                                                    int &bond_type_counter, TopologyFile* topology_file)
{
    vector<string> atom_pair_type = vector<string>();
    vector<string> reverse_atom_pair_type = vector<string>();
    stringstream key2;
    key2 << neighbor->GetId();
    atom_pair_type.push_back(assembly_atom->GetAtomType());
    atom_pair_type.push_back(neighbor->GetAtomType());
    reverse_atom_pair_type.push_back(neighbor->GetAtomType());
    reverse_atom_pair_type.push_back(assembly_atom->GetAtomType());

    if(find(inserted_bond_types.begin(), inserted_bond_types.end(), atom_pair_type) == inserted_bond_types.end() &&
            find(inserted_bond_types.begin(), inserted_bond_types.end(), reverse_atom_pair_type) == inserted_bond_types.end())
    {
        ParameterFileBond* parameter_file_bond;
        if(bonds.find(atom_pair_type) != bonds.end())
        {
            parameter_file_bond = bonds[atom_pair_type];
            inserted_bond_types.push_back(atom_pair_type);
        }
        else if(bonds.find(reverse_atom_pair_type) != bonds.end())
        {
            parameter_file_bond = bonds[reverse_atom_pair_type];
            inserted_bond_types.push_back(reverse_atom_pair_type);
        }
        else
        {
            //            stringstream ss;
            //            ss << atom_pair_type.at(0) << "-" << atom_pair_type.at(1) << " bond type does not exist in the parameter files";
            //            cout << ss.str() << endl;
            //            gmml::log(__LINE__, __FILE__, gmml::ERR, ss.str());
            return;
        }
        TopologyBondType* topology_bond_type = new TopologyBondType();
        topology_bond_type->SetForceConstant(parameter_file_bond->GetForceConstant());
        topology_bond_type->SetEquilibriumValue(parameter_file_bond->GetLength());
        topology_bond_type->SetIndex(bond_type_counter);
        bond_type_counter++;
        topology_file->AddBondType(topology_bond_type);
    }
}

void Assembly::ExtractTopologyBondsFromAssembly(vector<vector<string> > &inserted_bonds, vector<vector<string> > &inserted_bond_types,
                                                Atom *assembly_atom, Atom *neighbor, TopologyFileSpace::TopologyFile* topology_file)
{
    vector<string> atom_pair_type = vector<string>();
    vector<string> reverse_atom_pair_type = vector<string>();
    stringstream key2;
    key2 << neighbor->GetId();
    atom_pair_type.push_back(assembly_atom->GetAtomType());
    atom_pair_type.push_back(neighbor->GetAtomType());
    reverse_atom_pair_type.push_back(neighbor->GetAtomType());
    reverse_atom_pair_type.push_back(assembly_atom->GetAtomType());

    vector<string> atom_pair_name = vector<string>();;
    vector<string> reverse_atom_pair_name = vector<string>();;
    atom_pair_name.push_back(assembly_atom->GetName() + "(" + Split(assembly_atom->GetId(),"_").at(1) + ")");
    atom_pair_name.push_back(neighbor->GetName() + "(" + Split(neighbor->GetId(),"_").at(1) + ")");
    reverse_atom_pair_name.push_back(neighbor->GetName() + "(" + Split(neighbor->GetId(),"_").at(1) + ")");
    reverse_atom_pair_name.push_back(assembly_atom->GetName() + "(" + Split(assembly_atom->GetId(),"_").at(1) + ")");
    vector<string> residue_names = vector<string>();;
    vector<string> reverse_residue_names = vector<string>();;
    residue_names.push_back(assembly_atom->GetResidue()->GetName()+"("+Split(assembly_atom->GetResidue()->GetId(),"_").at(2)+")");
    residue_names.push_back(neighbor->GetResidue()->GetName()+"("+Split(neighbor->GetResidue()->GetId(),"_").at(2)+")");
    reverse_residue_names.push_back(neighbor->GetResidue()->GetName()+"("+Split(neighbor->GetResidue()->GetId(),"_").at(2)+")");
    reverse_residue_names.push_back(assembly_atom->GetResidue()->GetName()+"("+Split(assembly_atom->GetResidue()->GetId(),"_").at(2)+")");
    vector<string> bond = vector<string>();
    vector<string> reverse_bond = vector<string>();
    stringstream ss;
    ss << residue_names.at(0) << ":" << atom_pair_name.at(0);
    stringstream ss1;
    ss1 << residue_names.at(1) << ":" << atom_pair_name.at(1);
    bond.push_back(ss.str());
    bond.push_back(ss1.str());
    reverse_bond.push_back(ss1.str());
    reverse_bond.push_back(ss.str());


    if(find(inserted_bonds.begin(), inserted_bonds.end(), bond) == inserted_bonds.end() &&
            find(inserted_bonds.begin(), inserted_bonds.end(), reverse_bond) == inserted_bonds.end())
    {
        TopologyBond* topology_bond;
        if(find(inserted_bonds.begin(), inserted_bonds.end(), bond) == inserted_bonds.end())
        {
            topology_bond = new TopologyBond(atom_pair_name, residue_names);
            inserted_bonds.push_back(bond);
        }
        else if (find(inserted_bonds.begin(), inserted_bonds.end(), reverse_bond) == inserted_bonds.end())
        {
            topology_bond = new TopologyBond(reverse_atom_pair_name, reverse_residue_names);
            inserted_bonds.push_back(reverse_bond);
        }

        if((assembly_atom->GetName().substr(0,1).compare("H") == 0 ||
            (assembly_atom->GetName().substr(1,1).compare("H") == 0 && isdigit(ConvertString<char>(assembly_atom->GetName().substr(0,1)))))
                || (neighbor->GetName().substr(0,1).compare("H") == 0 ||
                    (neighbor->GetName().substr(1,1).compare("H") == 0 && isdigit(ConvertString<char>(neighbor->GetName().substr(0,1))))))
            topology_bond->SetIncludingHydrogen(true);
        else
            topology_bond->SetIncludingHydrogen(false);

        int index = 0;
        if(find(inserted_bond_types.begin(), inserted_bond_types.end(), atom_pair_type) != inserted_bond_types.end())
            index = distance(inserted_bond_types.begin(), find(inserted_bond_types.begin(), inserted_bond_types.end(), atom_pair_type));
        else if(find(inserted_bond_types.begin(), inserted_bond_types.end(), reverse_atom_pair_type) != inserted_bond_types.end())
            index = distance(inserted_bond_types.begin(), find(inserted_bond_types.begin(), inserted_bond_types.end(), reverse_atom_pair_type));
        else
        {
            //            stringstream ss;
            //            ss << atom_pair_type.at(0) << "-" << atom_pair_type.at(1) << " bond type does not exist in the parameter files";
            //            cout << ss.str() << endl;
            //            gmml::log(__LINE__, __FILE__, gmml::ERR, ss.str());
            return;
        }
        topology_bond->SetBondType(topology_file->GetBondTypeByIndex(index));
        topology_file->AddBond(topology_bond);
    }
}



void Assembly::ExtractTopologyAngleTypesFromAssembly(Atom* assembly_atom, Atom* neighbor, Atom* neighbor_of_neighbor, vector<vector<string> > &inserted_angle_types,
                                                     int &angle_type_counter, TopologyFile* topology_file, ParameterFileSpace::ParameterFile::AngleMap &angles)
{
    vector<string> angle_type = vector<string>();
    vector<string> reverse_angle_type = vector<string>();
    angle_type.push_back(assembly_atom->GetAtomType());
    angle_type.push_back(neighbor->GetAtomType());
    angle_type.push_back(neighbor_of_neighbor->GetAtomType());
    reverse_angle_type.push_back(neighbor_of_neighbor->GetAtomType());
    reverse_angle_type.push_back(neighbor->GetAtomType());
    reverse_angle_type.push_back(assembly_atom->GetAtomType());

    if(find(inserted_angle_types.begin(), inserted_angle_types.end(), angle_type) == inserted_angle_types.end() &&
            find(inserted_angle_types.begin(), inserted_angle_types.end(), reverse_angle_type) == inserted_angle_types.end())
    {
        ParameterFileAngle* parameter_file_angle;
        if(angles.find(angle_type) != angles.end())
        {
            parameter_file_angle = angles[angle_type];
            inserted_angle_types.push_back(angle_type);
        }
        else if(angles.find(reverse_angle_type) != angles.end())
        {
            parameter_file_angle = angles[reverse_angle_type];
            inserted_angle_types.push_back(reverse_angle_type);
        }
        else
        {
            //            stringstream ss;
            //            ss << angle_type.at(0) << "-" << angle_type.at(1) << "-" << angle_type.at(2) << " angle type does not exist in the parameter files";
            //            cout << ss.str() << endl;
            //            gmml::log(__LINE__, __FILE__, gmml::ERR, ss.str());
            return;
        }
        TopologyAngleType* topology_angle_type = new TopologyAngleType();
        topology_angle_type->SetForceConstant(parameter_file_angle->GetForceConstant());
        topology_angle_type->SetEquilibriumValue(parameter_file_angle->GetAngle());
        topology_angle_type->SetIndex(angle_type_counter);
        angle_type_counter++;
        topology_file->AddAngleType(topology_angle_type);
    }
}

void Assembly::ExtractTopologyAnglesFromAssembly(Atom* assembly_atom, Atom* neighbor, Atom* neighbor_of_neighbor, vector<vector<string> > &inserted_angles,
                                                 vector<vector<string> > &inserted_angle_types, TopologyFile* topology_file)
{
    vector<string> angle_type = vector<string>();
    vector<string> reverse_angle_type = vector<string>();
    angle_type.push_back(assembly_atom->GetAtomType());
    angle_type.push_back(neighbor->GetAtomType());
    angle_type.push_back(neighbor_of_neighbor->GetAtomType());
    reverse_angle_type.push_back(neighbor_of_neighbor->GetAtomType());
    reverse_angle_type.push_back(neighbor->GetAtomType());
    reverse_angle_type.push_back(assembly_atom->GetAtomType());

    vector<string> angle_atom_names = vector<string>();
    vector<string> reverse_angle_atom_names = vector<string>();
    angle_atom_names.push_back(assembly_atom->GetName() + "(" + Split(assembly_atom->GetId(),"_").at(1) + ")");
    angle_atom_names.push_back(neighbor->GetName() + "(" + Split(neighbor->GetId(),"_").at(1) + ")");
    angle_atom_names.push_back(neighbor_of_neighbor->GetName() + "(" + Split(neighbor_of_neighbor->GetId(),"_").at(1) + ")");
    reverse_angle_atom_names.push_back(neighbor_of_neighbor->GetName() + "(" + Split(neighbor_of_neighbor->GetId(),"_").at(1) + ")");
    reverse_angle_atom_names.push_back(neighbor->GetName() + "(" + Split(neighbor->GetId(),"_").at(1) + ")");
    reverse_angle_atom_names.push_back(assembly_atom->GetName() + "(" + Split(assembly_atom->GetId(),"_").at(1) + ")");

    vector<string> residue_names = vector<string>();
    vector<string> reverse_residue_names = vector<string>();
    residue_names.push_back(assembly_atom->GetResidue()->GetName()+"("+Split(assembly_atom->GetResidue()->GetId(),"_").at(2)+")");
    residue_names.push_back(neighbor->GetResidue()->GetName()+"("+Split(neighbor->GetResidue()->GetId(),"_").at(2)+")");
    residue_names.push_back(neighbor_of_neighbor->GetResidue()->GetName()+"("+Split(neighbor_of_neighbor->GetResidue()->GetId(),"_").at(2)+")");
    reverse_residue_names.push_back(neighbor_of_neighbor->GetResidue()->GetName()+"("+Split(neighbor_of_neighbor->GetResidue()->GetId(),"_").at(2)+")");
    reverse_residue_names.push_back(neighbor->GetResidue()->GetName()+"("+Split(neighbor->GetResidue()->GetId(),"_").at(2)+")");
    reverse_residue_names.push_back(assembly_atom->GetResidue()->GetName()+"("+Split(assembly_atom->GetResidue()->GetId(),"_").at(2)+")");
    vector<string> angle = vector<string>();
    vector<string> reverse_angle = vector<string>();
    stringstream ss;
    ss << residue_names.at(0) << ":" << angle_atom_names.at(0);
    stringstream ss1;
    ss1 << residue_names.at(1) << ":" << angle_atom_names.at(1);
    stringstream ss2;
    ss2 << residue_names.at(2) << ":" << angle_atom_names.at(2);
    angle.push_back(ss.str());
    angle.push_back(ss1.str());
    angle.push_back(ss2.str());
    reverse_angle.push_back(ss2.str());
    reverse_angle.push_back(ss1.str());
    reverse_angle.push_back(ss.str());

    if(find(inserted_angles.begin(), inserted_angles.end(), angle) == inserted_angles.end() &&
            find(inserted_angles.begin(), inserted_angles.end(), reverse_angle) == inserted_angles.end())
    {
        TopologyAngle* topology_angle;
        if(find(inserted_angles.begin(), inserted_angles.end(), angle) == inserted_angles.end())
        {
            topology_angle = new TopologyAngle(angle_atom_names, residue_names);
            inserted_angles.push_back(angle);
        }
        else if (find(inserted_angles.begin(), inserted_angles.end(), reverse_angle) == inserted_angles.end())
        {
            topology_angle = new TopologyAngle(reverse_angle_atom_names, reverse_residue_names);
            inserted_angles.push_back(reverse_angle);
        }
        if((assembly_atom->GetName().substr(0,1).compare("H") == 0 ||
            (assembly_atom->GetName().substr(1,1).compare("H") == 0 && isdigit(ConvertString<char>(assembly_atom->GetName().substr(0,1)))))
                || (neighbor->GetName().substr(0,1).compare("H") == 0 ||
                    (neighbor->GetName().substr(1,1).compare("H") == 0 && isdigit(ConvertString<char>(neighbor->GetName().substr(0,1)))))
                ||(neighbor_of_neighbor->GetName().substr(0,1).compare("H") == 0 ||
                   (neighbor_of_neighbor->GetName().substr(1,1).compare("H") == 0 && isdigit(ConvertString<char>(neighbor_of_neighbor->GetName().substr(0,1))))))
            topology_angle->SetIncludingHydrogen(true);
        else
            topology_angle->SetIncludingHydrogen(false);

        int index = 0;
        if(find(inserted_angle_types.begin(), inserted_angle_types.end(), angle_type) != inserted_angle_types.end())
            index = distance(inserted_angle_types.begin(), find(inserted_angle_types.begin(), inserted_angle_types.end(), angle_type));
        else if(find(inserted_angle_types.begin(), inserted_angle_types.end(), reverse_angle_type) != inserted_angle_types.end())
            index = distance(inserted_angle_types.begin(), find(inserted_angle_types.begin(), inserted_angle_types.end(), reverse_angle_type));
        else
        {
            //            stringstream ss;
            //            ss << angle_type.at(0) << "-" << angle_type.at(1) << "-" << angle_type.at(2) << " angle type does not exist in the parameter files";
            //            cout << ss.str() << endl;
            //            gmml::log(__LINE__, __FILE__, gmml::ERR, ss.str());
            return;
        }
        topology_angle->SetAnlgeType(topology_file->GetAngleTypeByIndex(index));
        topology_file->AddAngle(topology_angle);
    }
}

void Assembly::ExtractTopologyDihedralTypesFromAssembly(Atom *assembly_atom, Atom *neighbor, Atom *neighbor_of_neighbor, Atom *neighbor_of_neighbor_of_neighbor,
                                                        vector<string>& inserted_dihedral_types, int &dihedral_type_counter, TopologyFile *topology_file, ParameterFile::DihedralMap& dihedrals)
{
    vector<vector<string> > all_atom_type_permutations = CreateAllAtomTypePermutationsforDihedralType(assembly_atom->GetAtomType(), neighbor->GetAtomType(),
                                                                                                      neighbor_of_neighbor->GetAtomType(), neighbor_of_neighbor_of_neighbor->GetAtomType());
    bool is_found = false;
    for(vector<vector<string> >::iterator it = all_atom_type_permutations.begin(); it != all_atom_type_permutations.end(); it++)
    {
        vector<string> atom_types = (*it);
        if(dihedrals[atom_types] != NULL)
        {
            stringstream ss;
            ss << atom_types.at(0) << "_" << atom_types.at(1) << "_" << atom_types.at(2) << "_" << atom_types.at(3);
            if(find(inserted_dihedral_types.begin(), inserted_dihedral_types.end(), ss.str()) == inserted_dihedral_types.end())
            {
                ParameterFileDihedral* parameter_file_dihedral = dihedrals[atom_types];
                vector<ParameterFileDihedralTerm> dihedral_terms = parameter_file_dihedral->GetTerms();
                for(vector<ParameterFileDihedralTerm>::iterator it1 = dihedral_terms.begin(); it1 != dihedral_terms.end(); it1++)
                {
                    inserted_dihedral_types.push_back(ss.str());
                    ParameterFileDihedralTerm parameter_file_dihedral_term = (*it1);
                    TopologyDihedralType* topology_dihedral_type = new TopologyDihedralType();
                    topology_dihedral_type->SetIndex(dihedral_type_counter);
                    dihedral_type_counter++;
                    topology_dihedral_type->SetForceConstant(parameter_file_dihedral_term.GetForceConstant());
                    topology_dihedral_type->SetPeriodicity(fabs(parameter_file_dihedral_term.GetPeriodicity()));
                    topology_dihedral_type->SetPhase(parameter_file_dihedral_term.GetPhase());
                    topology_dihedral_type->SetScee(parameter_file_dihedral->GetScee());
                    topology_dihedral_type->SetScnb(parameter_file_dihedral->GetScnb());
                    topology_file->AddDihedralType(topology_dihedral_type);
                }
                is_found = true;
                break;
            }
        }
    }
    if(!is_found)
    {
        //        stringstream ss;
        //        ss << all_atom_type_permutations.at(0).at(0) << "-" << all_atom_type_permutations.at(0).at(1) << "-" << all_atom_type_permutations.at(0).at(2) << "-"
        //           << all_atom_type_permutations.at(0).at(3) << " dihedral type (or any other permutation of it) does not exist in the parameter files";
        //        cout << ss.str() << endl;
        //        gmml::log(__LINE__, __FILE__, gmml::ERR, ss.str());
    }

    ///Improper Dihedrals
    AtomNode* atom_node = assembly_atom->GetNode();
    AtomVector neighbors = atom_node->GetNodeNeighbors();
    if(neighbors.size() == 3)
    {
        Atom* neighbor1 = neighbors.at(0);
        Atom* neighbor2 = neighbors.at(1);
        Atom* neighbor3 = neighbors.at(2);
        vector<vector<string> > all_improper_dihedrals_atom_type_permutations = CreateAllAtomTypePermutationsforImproperDihedralType(neighbor1->GetAtomType(), neighbor2->GetAtomType(),
                                                                                                                                     neighbor3->GetAtomType(), assembly_atom->GetAtomType());
        bool is_improper_found = false;
        for(vector<vector<string> >::iterator it = all_improper_dihedrals_atom_type_permutations.begin(); it != all_improper_dihedrals_atom_type_permutations.end(); it++)
        {
            vector<string> improper_dihedral_permutation = (*it);
            if(dihedrals[improper_dihedral_permutation] != NULL)
            {
                stringstream ss;
                ss << improper_dihedral_permutation.at(0) << "_" << improper_dihedral_permutation.at(1) << "_" << improper_dihedral_permutation.at(2) << "_" << improper_dihedral_permutation.at(3);
                if(find(inserted_dihedral_types.begin(), inserted_dihedral_types.end(), ss.str()) == inserted_dihedral_types.end())
                {
                    ParameterFileDihedral* parameter_file_dihedral = dihedrals[improper_dihedral_permutation];
                    vector<ParameterFileDihedralTerm> dihedral_terms = parameter_file_dihedral->GetTerms();
                    for(vector<ParameterFileDihedralTerm>::iterator it1 = dihedral_terms.begin(); it1 != dihedral_terms.end(); it1++)
                    {
                        inserted_dihedral_types.push_back(ss.str());
                        ParameterFileDihedralTerm parameter_file_dihedral_term = (*it1);
                        TopologyDihedralType* topology_dihedral_type = new TopologyDihedralType();
                        topology_dihedral_type->SetIndex(dihedral_type_counter);
                        dihedral_type_counter++;
                        topology_dihedral_type->SetForceConstant(parameter_file_dihedral_term.GetForceConstant());
                        topology_dihedral_type->SetPeriodicity(fabs(parameter_file_dihedral_term.GetPeriodicity()));
                        topology_dihedral_type->SetPhase(parameter_file_dihedral_term.GetPhase());
                        topology_dihedral_type->SetScee(parameter_file_dihedral->GetScee());
                        topology_dihedral_type->SetScnb(parameter_file_dihedral->GetScnb());
                        topology_file->AddDihedralType(topology_dihedral_type);
                    }
                    is_improper_found = true;
                    break;
                }
            }
        }
        if(!is_improper_found)
        {
            //            stringstream ss;
            //            ss << all_improper_dihedrals_atom_type_permutations.at(0).at(0) << "-" << all_improper_dihedrals_atom_type_permutations.at(0).at(1) << "-"
            //               << all_improper_dihedrals_atom_type_permutations.at(0).at(2) << "-" << all_improper_dihedrals_atom_type_permutations.at(0).at(3)
            //               << " improer dihedral type (or any other permutation of it) does not exist in the parameter files";
            //            cout << ss.str() << endl;
            //            gmml::log(__LINE__, __FILE__, gmml::ERR, ss.str());
        }
    }
}

void Assembly::ExtractTopologyDihedralsFromAssembly(Atom *assembly_atom, Atom *neighbor, Atom *neighbor_of_neighbor, Atom *neighbor_of_neighbor_of_neighbor,
                                                    vector<vector<string> >& inserted_dihedrals, vector<string>& inserted_dihedral_types,
                                                    ParameterFile::DihedralMap &dihedrals, TopologyFile *topology_file)
{
    vector<vector<string> > all_atom_type_permutations = CreateAllAtomTypePermutationsforDihedralType(assembly_atom->GetAtomType(), neighbor->GetAtomType(),
                                                                                                      neighbor_of_neighbor->GetAtomType(), neighbor_of_neighbor_of_neighbor->GetAtomType());
    for(vector<vector<string> >::iterator it = all_atom_type_permutations.begin(); it != all_atom_type_permutations.end(); it++)
    {
        vector<string> atom_types = (*it);
        stringstream sss;
        sss << atom_types.at(0) << "_" << atom_types.at(1) << "_" << atom_types.at(2) << "_" << atom_types.at(3);
        if(find(inserted_dihedral_types.begin(), inserted_dihedral_types.end(), sss.str()) != inserted_dihedral_types.end())
        {
            vector<string> dihedral_atom_names = vector<string>();
            vector<string> reverse_dihedral_atom_names = vector<string>();
            dihedral_atom_names.push_back(assembly_atom->GetName() + "(" + Split(assembly_atom->GetId(),"_").at(1) + ")");
            dihedral_atom_names.push_back(neighbor->GetName() + "(" + Split(neighbor->GetId(),"_").at(1) + ")");
            dihedral_atom_names.push_back(neighbor_of_neighbor->GetName() + "(" + Split(neighbor_of_neighbor->GetId(),"_").at(1) + ")");
            dihedral_atom_names.push_back(neighbor_of_neighbor_of_neighbor->GetName() + "(" + Split(neighbor_of_neighbor_of_neighbor->GetId(),"_").at(1) + ")");
            reverse_dihedral_atom_names.push_back(neighbor_of_neighbor_of_neighbor->GetName() + "(" + Split(neighbor_of_neighbor_of_neighbor->GetId(),"_").at(1) + ")");
            reverse_dihedral_atom_names.push_back(neighbor_of_neighbor->GetName() + "(" + Split(neighbor_of_neighbor->GetId(),"_").at(1) + ")");
            reverse_dihedral_atom_names.push_back(neighbor->GetName() + "(" + Split(neighbor->GetId(),"_").at(1) + ")");
            reverse_dihedral_atom_names.push_back(assembly_atom->GetName() + "(" + Split(assembly_atom->GetId(),"_").at(1) + ")");

            vector<string> residue_names = vector<string>();
            vector<string> reverse_residue_names = vector<string>();
            residue_names.push_back(assembly_atom->GetResidue()->GetName()+"("+Split(assembly_atom->GetResidue()->GetId(),"_").at(2)+")");
            residue_names.push_back(neighbor->GetResidue()->GetName()+"("+Split(neighbor->GetResidue()->GetId(),"_").at(2)+")");
            residue_names.push_back(neighbor_of_neighbor->GetResidue()->GetName()+"("+Split(neighbor_of_neighbor->GetResidue()->GetId(),"_").at(2)+")");
            residue_names.push_back(neighbor_of_neighbor_of_neighbor->GetResidue()->GetName()+"("+Split(neighbor_of_neighbor_of_neighbor->GetResidue()->GetId(),"_").at(2)+")");
            reverse_residue_names.push_back(neighbor_of_neighbor_of_neighbor->GetResidue()->GetName()+"("+Split(neighbor_of_neighbor_of_neighbor->GetResidue()->GetId(),"_").at(2)+")");
            reverse_residue_names.push_back(neighbor_of_neighbor->GetResidue()->GetName()+"("+Split(neighbor_of_neighbor->GetResidue()->GetId(),"_").at(2)+")");
            reverse_residue_names.push_back(neighbor->GetResidue()->GetName()+"("+Split(neighbor->GetResidue()->GetId(),"_").at(2)+")");
            reverse_residue_names.push_back(assembly_atom->GetResidue()->GetName()+"("+Split(assembly_atom->GetResidue()->GetId(),"_").at(2)+")");

            vector<string> dihedral = vector<string>();
            vector<string> reverse_dihedral = vector<string>();
            stringstream ss;
            ss << residue_names.at(0) << ":" << dihedral_atom_names.at(0);
            stringstream ss1;
            ss1 << residue_names.at(1) << ":" << dihedral_atom_names.at(1);
            stringstream ss2;
            ss2 << residue_names.at(2) << ":" << dihedral_atom_names.at(2);
            stringstream ss3;
            ss3 << residue_names.at(3) << ":" << dihedral_atom_names.at(3);
            dihedral.push_back(ss.str());
            dihedral.push_back(ss1.str());
            dihedral.push_back(ss2.str());
            dihedral.push_back(ss3.str());
            reverse_dihedral.push_back(ss3.str());
            reverse_dihedral.push_back(ss2.str());
            reverse_dihedral.push_back(ss1.str());
            reverse_dihedral.push_back(ss.str());

            if(find(inserted_dihedrals.begin(), inserted_dihedrals.end(), dihedral) == inserted_dihedrals.end() &&
                    find(inserted_dihedrals.begin(), inserted_dihedrals.end(), reverse_dihedral) == inserted_dihedrals.end())
            {
                ParameterFileDihedral* parameter_file_dihedral = dihedrals[atom_types];
                vector<ParameterFileDihedralTerm> dihedral_terms = parameter_file_dihedral->GetTerms();
                for(vector<ParameterFileDihedralTerm>::iterator it1 = dihedral_terms.begin(); it1 != dihedral_terms.end(); it1++)
                {
                    TopologyDihedral* topology_dihedral = new TopologyDihedral();
                    topology_dihedral->SetIsImproper(false);
                    topology_dihedral->SetIgnoredGroupInteraction(false);///not sure
                    if((assembly_atom->GetName().substr(0,1).compare("H") == 0 ||
                        (assembly_atom->GetName().substr(1,1).compare("H") == 0 && isdigit(ConvertString<char>(assembly_atom->GetName().substr(0,1)))))
                            || (neighbor->GetName().substr(0,1).compare("H") == 0 ||
                                (neighbor->GetName().substr(1,1).compare("H") == 0 && isdigit(ConvertString<char>(neighbor->GetName().substr(0,1)))))
                            ||(neighbor_of_neighbor->GetName().substr(0,1).compare("H") == 0 ||
                               (neighbor_of_neighbor->GetName().substr(1,1).compare("H") == 0 && isdigit(ConvertString<char>(neighbor_of_neighbor->GetName().substr(0,1)))))
                            ||(neighbor_of_neighbor_of_neighbor->GetName().substr(0,1).compare("H") == 0 ||
                               (neighbor_of_neighbor_of_neighbor->GetName().substr(1,1).compare("H") == 0 && isdigit(ConvertString<char>(neighbor_of_neighbor_of_neighbor->GetName().substr(0,1))))))
                        topology_dihedral->SetIncludingHydrogen(true);
                    else
                        topology_dihedral->SetIncludingHydrogen(false);

                    if(atom_types.at(0).compare(assembly_atom->GetAtomType()) == 0 ||
                            atom_types.at(0).compare("X") == 0 && atom_types.at(3).compare(neighbor_of_neighbor_of_neighbor->GetAtomType()) == 0 ||
                            atom_types.at(0).compare("X") == 0 && atom_types.at(3).compare("X") == 0 && atom_types.at(1).compare(neighbor->GetAtomType()) == 0 ||
                            atom_types.at(0).compare("X") == 0 && atom_types.at(3).compare("X") == 0 && atom_types.at(2).compare(neighbor_of_neighbor->GetAtomType()) == 0 )
                    {
                        topology_dihedral->SetResidueNames(residue_names);
                        topology_dihedral->SetDihedrals(dihedral_atom_names);
                    }
                    else
                    {
                        topology_dihedral->SetResidueNames(reverse_residue_names);
                        topology_dihedral->SetDihedrals(reverse_dihedral_atom_names);
                    }

                    int index = 0;
                    if(find(inserted_dihedral_types.begin(), inserted_dihedral_types.end(), sss.str()) != inserted_dihedral_types.end())
                        index = distance(inserted_dihedral_types.begin(), find(inserted_dihedral_types.begin(), inserted_dihedral_types.end(), sss.str())) +
                                distance(dihedral_terms.begin(), it1);
                    topology_dihedral->SetDihedralType(topology_file->GetDihedralTypeByIndex(index));
                    topology_file->AddDihedral(topology_dihedral);
                }
                if(atom_types.at(0).compare(assembly_atom->GetAtomType()) == 0 ||
                        (atom_types.at(0).compare("X") == 0 && atom_types.at(3).compare(neighbor_of_neighbor_of_neighbor->GetAtomType()) == 0) ||
                        (atom_types.at(0).compare("X") == 0 && atom_types.at(3).compare("X") == 0 && atom_types.at(1).compare(neighbor->GetAtomType()) == 0) ||
                        (atom_types.at(0).compare("X") == 0 && atom_types.at(3).compare("X") == 0 && atom_types.at(2).compare(neighbor_of_neighbor->GetAtomType()) == 0) )
                    inserted_dihedrals.push_back(dihedral);
                else
                    inserted_dihedrals.push_back(reverse_dihedral);

                break;
            }
        }
    }
    /**/
    ///Improper Dihedrals
    AtomNode* atom_node = assembly_atom->GetNode();
    AtomVector neighbors = atom_node->GetNodeNeighbors();
    if(neighbors.size() == 3)
    {
        Atom* neighbor1 = neighbors.at(0);
        Atom* neighbor2 = neighbors.at(1);
        Atom* neighbor3 = neighbors.at(2);
        vector<vector<string> > all_improper_dihedrals_atom_type_permutations = CreateAllAtomTypePermutationsforImproperDihedralType(neighbor1->GetAtomType(), neighbor2->GetAtomType(),
                                                                                                                                     neighbor3->GetAtomType(), assembly_atom->GetAtomType());
        for(vector<vector<string> >::iterator it = all_improper_dihedrals_atom_type_permutations.begin(); it != all_improper_dihedrals_atom_type_permutations.end(); it++)
        {
            vector<string> improper_dihedral_permutation = (*it);
            stringstream sss;
            sss << improper_dihedral_permutation.at(0) << "_" << improper_dihedral_permutation.at(1) << "_" << improper_dihedral_permutation.at(2) << "_" << improper_dihedral_permutation.at(3);
            if(find(inserted_dihedral_types.begin(), inserted_dihedral_types.end(), sss.str()) != inserted_dihedral_types.end())
            {
                vector<string> dihedral_atom_names1 = vector<string>();
                dihedral_atom_names1.push_back(neighbor1->GetName() + "(" + Split(neighbor1->GetId(),"_").at(1) + ")");
                dihedral_atom_names1.push_back(neighbor2->GetName() + "(" + Split(neighbor2->GetId(),"_").at(1) + ")");
                dihedral_atom_names1.push_back(assembly_atom->GetName() + "(" + Split(assembly_atom->GetId(),"_").at(1) + ")");
                dihedral_atom_names1.push_back(neighbor3->GetName() + "(" + Split(neighbor3->GetId(),"_").at(1) + ")");
                //                vector<string> dihedral_atom_names2 = vector<string>();
                //                dihedral_atom_names2.push_back(neighbor1->GetName() + "(" + Split(neighbor1->GetId(),"_").at(1) + ")");
                //                dihedral_atom_names2.push_back(assembly_atom->GetName() + "(" + Split(assembly_atom->GetId(),"_").at(1) + ")");
                //                dihedral_atom_names2.push_back(neighbor3->GetName() + "(" + Split(neighbor3->GetId(),"_").at(1) + ")");
                //                dihedral_atom_names2.push_back(neighbor2->GetName() + "(" + Split(neighbor2->GetId(),"_").at(1) + ")");
                vector<string> dihedral_atom_names3 = vector<string>();
                dihedral_atom_names3.push_back(neighbor1->GetName() + "(" + Split(neighbor1->GetId(),"_").at(1) + ")");
                dihedral_atom_names3.push_back(neighbor3->GetName() + "(" + Split(neighbor3->GetId(),"_").at(1) + ")");
                dihedral_atom_names3.push_back(assembly_atom->GetName() + "(" + Split(assembly_atom->GetId(),"_").at(1) + ")");
                dihedral_atom_names3.push_back(neighbor2->GetName() + "(" + Split(neighbor2->GetId(),"_").at(1) + ")");

                //                vector<string> reverse_dihedral_atom_names1 = vector<string>();
                //                reverse_dihedral_atom_names1.push_back(neighbor3->GetName() + "(" + Split(neighbor3->GetId(),"_").at(1) + ")");
                //                reverse_dihedral_atom_names1.push_back(assembly_atom->GetName() + "(" + Split(assembly_atom->GetId(),"_").at(1) + ")");
                //                reverse_dihedral_atom_names1.push_back(neighbor2->GetName() + "(" + Split(neighbor2->GetId(),"_").at(1) + ")");
                //                reverse_dihedral_atom_names1.push_back(neighbor1->GetName() + "(" + Split(neighbor1->GetId(),"_").at(1) + ")");
                vector<string> reverse_dihedral_atom_names2 = vector<string>();
                reverse_dihedral_atom_names2.push_back(neighbor2->GetName() + "(" + Split(neighbor2->GetId(),"_").at(1) + ")");
                reverse_dihedral_atom_names2.push_back(neighbor3->GetName() + "(" + Split(neighbor3->GetId(),"_").at(1) + ")");
                reverse_dihedral_atom_names2.push_back(assembly_atom->GetName() + "(" + Split(assembly_atom->GetId(),"_").at(1) + ")");
                reverse_dihedral_atom_names2.push_back(neighbor1->GetName() + "(" + Split(neighbor1->GetId(),"_").at(1) + ")");
                //                vector<string> reverse_dihedral_atom_names3 = vector<string>();
                //                reverse_dihedral_atom_names3.push_back(neighbor2->GetName() + "(" + Split(neighbor2->GetId(),"_").at(1) + ")");
                //                reverse_dihedral_atom_names3.push_back(assembly_atom->GetName() + "(" + Split(assembly_atom->GetId(),"_").at(1) + ")");
                //                reverse_dihedral_atom_names3.push_back(neighbor3->GetName() + "(" + Split(neighbor3->GetId(),"_").at(1) + ")");
                //                reverse_dihedral_atom_names3.push_back(neighbor1->GetName() + "(" + Split(neighbor1->GetId(),"_").at(1) + ")");

                vector<string> residue_names1 = vector<string>();
                residue_names1.push_back(neighbor1->GetResidue()->GetName()+"("+Split(neighbor1->GetResidue()->GetId(),"_").at(2)+")");
                residue_names1.push_back(neighbor2->GetResidue()->GetName()+"("+Split(neighbor2->GetResidue()->GetId(),"_").at(2)+")");
                residue_names1.push_back(assembly_atom->GetResidue()->GetName()+"("+Split(assembly_atom->GetResidue()->GetId(),"_").at(2)+")");
                residue_names1.push_back(neighbor3->GetResidue()->GetName()+"("+Split(neighbor3->GetResidue()->GetId(),"_").at(2)+")");
                //                vector<string> residue_names2 = vector<string>();
                //                residue_names2.push_back(neighbor1->GetName()+"("+Split(neighbor1->GetResidue()->GetId(),"_").at(2)+")");
                //                residue_names2.push_back(assembly_atom->GetResidue()->GetName()+"("+Split(assembly_atom->GetResidue()->GetId(),"_").at(2)+")");
                //                residue_names2.push_back(neighbor3->GetResidue()->GetName()+"("+Split(neighbor3->GetResidue()->GetId(),"_").at(2)+")");
                //                residue_names2.push_back(neighbor2->GetResidue()->GetName()+"("+Split(neighbor2->GetResidue()->GetId(),"_").at(2)+")");
                vector<string> residue_names3 = vector<string>();
                residue_names3.push_back(neighbor1->GetResidue()->GetName()+"("+Split(neighbor1->GetResidue()->GetId(),"_").at(2)+")");
                residue_names3.push_back(neighbor3->GetResidue()->GetName()+"("+Split(neighbor3->GetResidue()->GetId(),"_").at(2)+")");
                residue_names3.push_back(assembly_atom->GetResidue()->GetName()+"("+Split(assembly_atom->GetResidue()->GetId(),"_").at(2)+")");
                residue_names3.push_back(neighbor2->GetResidue()->GetName()+"("+Split(neighbor2->GetResidue()->GetId(),"_").at(2)+")");

                //                vector<string> reverse_residue_names1 = vector<string>();
                //                reverse_residue_names1.push_back(neighbor3->GetResidue()->GetName()+"("+Split(neighbor3->GetResidue()->GetId(),"_").at(2)+")");
                //                reverse_residue_names1.push_back(assembly_atom->GetResidue()->GetName()+"("+Split(assembly_atom->GetResidue()->GetId(),"_").at(2)+")");
                //                reverse_residue_names1.push_back(neighbor2->GetResidue()->GetName()+"("+Split(neighbor2->GetResidue()->GetId(),"_").at(2)+")");
                //                reverse_residue_names1.push_back(neighbor1->GetResidue()->GetName()+"("+Split(neighbor1->GetResidue()->GetId(),"_").at(2)+")");
                vector<string> reverse_residue_names2 = vector<string>();
                reverse_residue_names2.push_back(neighbor2->GetResidue()->GetName()+"("+Split(neighbor2->GetResidue()->GetId(),"_").at(2)+")");
                reverse_residue_names2.push_back(neighbor3->GetResidue()->GetName()+"("+Split(neighbor3->GetResidue()->GetId(),"_").at(2)+")");
                reverse_residue_names2.push_back(assembly_atom->GetResidue()->GetName()+"("+Split(assembly_atom->GetResidue()->GetId(),"_").at(2)+")");
                reverse_residue_names2.push_back(neighbor1->GetName()+"("+Split(neighbor1->GetResidue()->GetId(),"_").at(2)+")");
                //                vector<string> reverse_residue_names3 = vector<string>();
                //                reverse_residue_names3.push_back(neighbor2->GetResidue()->GetName()+"("+Split(neighbor2->GetResidue()->GetId(),"_").at(2)+")");
                //                reverse_residue_names3.push_back(assembly_atom->GetResidue()->GetName()+"("+Split(assembly_atom->GetResidue()->GetId(),"_").at(2)+")");
                //                reverse_residue_names3.push_back(neighbor3->GetResidue()->GetName()+"("+Split(neighbor3->GetResidue()->GetId(),"_").at(2)+")");
                //                reverse_residue_names3.push_back(neighbor1->GetResidue()->GetName()+"("+Split(neighbor1->GetResidue()->GetId(),"_").at(2)+")");

                vector<string> dihedral1 = vector<string>();
                vector<string> dihedral2 = vector<string>();
                vector<string> dihedral3 = vector<string>();
                vector<string> reverse_dihedral1 = vector<string>();
                vector<string> reverse_dihedral2 = vector<string>();
                vector<string> reverse_dihedral3 = vector<string>();
                stringstream ss;
                ss << residue_names1.at(2) << ":" << dihedral_atom_names1.at(2);
                stringstream ss1;
                ss1 << residue_names1.at(0) << ":" << dihedral_atom_names1.at(0);
                stringstream ss2;
                ss2 << residue_names1.at(1) << ":" << dihedral_atom_names1.at(1);
                stringstream ss3;
                ss3 << residue_names1.at(3) << ":" << dihedral_atom_names1.at(3);

                dihedral1.push_back(ss1.str());
                dihedral1.push_back(ss2.str());
                dihedral1.push_back(ss.str());
                dihedral1.push_back(ss3.str());
                reverse_dihedral1.push_back(ss3.str());
                reverse_dihedral1.push_back(ss.str());
                reverse_dihedral1.push_back(ss2.str());
                reverse_dihedral1.push_back(ss1.str());

                dihedral2.push_back(ss1.str());
                dihedral2.push_back(ss.str());
                dihedral2.push_back(ss3.str());
                dihedral2.push_back(ss2.str());
                reverse_dihedral2.push_back(ss2.str());
                reverse_dihedral2.push_back(ss3.str());
                reverse_dihedral2.push_back(ss.str());
                reverse_dihedral2.push_back(ss1.str());

                dihedral3.push_back(ss1.str());
                dihedral3.push_back(ss3.str());
                dihedral3.push_back(ss.str());
                dihedral3.push_back(ss2.str());
                reverse_dihedral3.push_back(ss2.str());
                reverse_dihedral3.push_back(ss.str());
                reverse_dihedral3.push_back(ss3.str());
                reverse_dihedral3.push_back(ss1.str());

                if(find(inserted_dihedrals.begin(), inserted_dihedrals.end(), dihedral1) == inserted_dihedrals.end() &&
                        find(inserted_dihedrals.begin(), inserted_dihedrals.end(), dihedral2) == inserted_dihedrals.end() &&
                        find(inserted_dihedrals.begin(), inserted_dihedrals.end(), dihedral3) == inserted_dihedrals.end() &&
                        find(inserted_dihedrals.begin(), inserted_dihedrals.end(), reverse_dihedral1) == inserted_dihedrals.end() &&
                        find(inserted_dihedrals.begin(), inserted_dihedrals.end(), reverse_dihedral2) == inserted_dihedrals.end() &&
                        find(inserted_dihedrals.begin(), inserted_dihedrals.end(), reverse_dihedral3) == inserted_dihedrals.end())
                {
                    int permutation_index = distance(all_improper_dihedrals_atom_type_permutations.begin(), it);
                    ParameterFileDihedral* parameter_file_dihedral = dihedrals[improper_dihedral_permutation];
                    vector<ParameterFileDihedralTerm> dihedral_terms = parameter_file_dihedral->GetTerms();
                    for(vector<ParameterFileDihedralTerm>::iterator it1 = dihedral_terms.begin(); it1 != dihedral_terms.end(); it1++)
                    {
                        TopologyDihedral* topology_dihedral = new TopologyDihedral();
                        topology_dihedral->SetIsImproper(true);
                        topology_dihedral->SetIgnoredGroupInteraction(false);///not sure

                        if((assembly_atom->GetName().substr(0,1).compare("H") == 0 ||
                            (assembly_atom->GetName().substr(1,1).compare("H") == 0 && isdigit(ConvertString<char>(assembly_atom->GetName().substr(0,1)))))
                                || (neighbor1->GetName().substr(0,1).compare("H") == 0 ||
                                    (neighbor1->GetName().substr(1,1).compare("H") == 0 && isdigit(ConvertString<char>(neighbor1->GetName().substr(0,1)))))
                                ||(neighbor2->GetName().substr(0,1).compare("H") == 0 ||
                                   (neighbor2->GetName().substr(1,1).compare("H") == 0 && isdigit(ConvertString<char>(neighbor2->GetName().substr(0,1)))))
                                ||(neighbor3->GetName().substr(0,1).compare("H") == 0 ||
                                   (neighbor3->GetName().substr(1,1).compare("H") == 0 && isdigit(ConvertString<char>(neighbor3->GetName().substr(0,1))))))
                            topology_dihedral->SetIncludingHydrogen(true);
                        else
                            topology_dihedral->SetIncludingHydrogen(false);

                        if(permutation_index % 6 == 0)
                        {
                            topology_dihedral->SetResidueNames(residue_names1);
                            topology_dihedral->SetDihedrals(dihedral_atom_names1);
                        }
                        //                        if(permutation_index % 6 == 2)
                        //                        {
                        //                            topology_dihedral->SetResidueNames(residue_names2);
                        //                            topology_dihedral->SetDihedrals(dihedral_atom_names2);
                        //                        }
                        if(permutation_index % 6 == 4)
                        {
                            topology_dihedral->SetResidueNames(residue_names3);
                            topology_dihedral->SetDihedrals(dihedral_atom_names3);
                        }
                        //                        if(permutation_index % 6 == 1)
                        //                        {
                        //                            topology_dihedral->SetResidueNames(reverse_residue_names1);
                        //                            topology_dihedral->SetDihedrals(reverse_dihedral_atom_names1);
                        //                        }
                        if(permutation_index % 6 == 3)
                        {
                            topology_dihedral->SetResidueNames(reverse_residue_names2);
                            topology_dihedral->SetDihedrals(reverse_dihedral_atom_names2);
                        }
                        //                        if(permutation_index % 6 == 5)
                        //                        {
                        //                            topology_dihedral->SetResidueNames(reverse_residue_names3);
                        //                            topology_dihedral->SetDihedrals(reverse_dihedral_atom_names3);
                        //                        }

                        int index = 0;
                        if(find(inserted_dihedral_types.begin(), inserted_dihedral_types.end(), sss.str()) != inserted_dihedral_types.end())
                            index = distance(inserted_dihedral_types.begin(), find(inserted_dihedral_types.begin(), inserted_dihedral_types.end(), sss.str())) +
                                    distance(dihedral_terms.begin(), it1);
                        topology_dihedral->SetDihedralType(topology_file->GetDihedralTypeByIndex(index));
                        topology_file->AddDihedral(topology_dihedral);
                    }
                    if(permutation_index % 6 == 0)
                    {
                        inserted_dihedrals.push_back(dihedral1);
                    }
                    if(permutation_index % 6 == 2)
                    {
                        inserted_dihedrals.push_back(dihedral2);
                    }
                    if(permutation_index % 6 == 4)
                    {
                        inserted_dihedrals.push_back(dihedral3);
                    }
                    if(permutation_index % 6 == 1)
                    {
                        inserted_dihedrals.push_back(reverse_dihedral1);
                    }
                    if(permutation_index % 6 == 3)
                    {
                        inserted_dihedrals.push_back(reverse_dihedral2);
                    }
                    if(permutation_index % 6 == 5)
                    {
                        inserted_dihedrals.push_back(reverse_dihedral3);
                    }
                    break;
                }
            }
        }
    }
    /**/
}

CoordinateFile* Assembly::BuildCoordinateFileStructureFromAssembly()
{
    cout << "Creating coordinate file ..." << endl;
    gmml::log(__LINE__, __FILE__, gmml::INF, "Creating coordinate file ...");
    vector<Coordinate*> coordinates = this->GetAllCoordinates();
    CoordinateFile* coordinate_file = new CoordinateFile();
    coordinate_file->SetCoordinates(coordinates);
    coordinate_file->SetNumberOfCoordinates(coordinates.size());
    string title = "Generated by GMML";
    coordinate_file->SetTitle(title);
    return coordinate_file;
}

LibraryFile* Assembly::BuildLibraryFileStructureFromAssembly()
{
    cout << "Creating library file ..." << endl;
    gmml::log(__LINE__, __FILE__, gmml::INF, "Creating library file ...");
    LibraryFile* library_file = new LibraryFile();
    LibraryFile::ResidueMap residue_map = LibraryFile::ResidueMap();
    ResidueVector residues_of_assembly = this->GetAllResiduesOfAssembly();
    for(ResidueVector::iterator it = residues_of_assembly.begin(); it != residues_of_assembly.end(); it++)
    {
        Residue* assembly_residue = *it;
        //        cout << assembly_residue->GetId() << endl;
        int residue_index = distance(residues_of_assembly.begin(), it) + 1;
        AtomVector assembly_residue_atoms = assembly_residue->GetAtoms();
        LibraryFileResidue* library_residue = new LibraryFileResidue();
        library_residue->SetName(assembly_residue->GetName());
        AtomVector head_atoms = assembly_residue->GetHeadAtoms();
        AtomVector tail_atoms = assembly_residue->GetTailAtoms();
        if(!head_atoms.empty())
        {
            int head_atom_index = distance(assembly_residue_atoms.begin(), find(assembly_residue_atoms.begin(), assembly_residue_atoms.end(), head_atoms.at(0))) + 1;
            library_residue->SetHeadAtomIndex(head_atom_index);
        }
        if(!tail_atoms.empty())
        {
            int tail_atom_index = distance(assembly_residue_atoms.begin(), find(assembly_residue_atoms.begin(), assembly_residue_atoms.end(), tail_atoms.at(0))) + 1;
            library_residue->SetTailAtomIndex(tail_atom_index);
        }
        int order = 1;
        for(AtomVector::iterator it1 = assembly_residue_atoms.begin(); it1 != assembly_residue_atoms.end(); it1++)
        {
            Atom* residue_atom = (*it1);
            int atom_index = distance(assembly_residue_atoms.begin(), it1) + 1;
            vector<int> bonded_atom_indices = vector<int>();
            AtomNode* atom_node = residue_atom->GetNode();
            if(atom_node != NULL)
            {
                AtomVector atom_neighbours = atom_node->GetNodeNeighbors();
                for(AtomVector::iterator it2 = atom_neighbours.begin(); it2 != atom_neighbours.end(); it2++)
                {
                    Atom* atom = *it2;
                    int bonded_atom_index = distance(assembly_residue_atoms.begin(), find(assembly_residue_atoms.begin(), assembly_residue_atoms.end(),
                                                                                          atom));
                    bonded_atom_indices.push_back(bonded_atom_index);
                }
            }
            LibraryFileAtom* atom = new LibraryFileAtom(residue_atom->GetAtomType(), residue_atom->GetName(), residue_index, atom_index,
                                                        gmml::iNotSet, residue_atom->MolecularDynamicAtom::GetCharge(),
                                                        *(residue_atom->GetCoordinates()[assembly_residue->GetAssembly()->GetModelIndex()]), bonded_atom_indices,
                    order);
            order++;
            library_residue->AddAtom(atom);
        }
        residue_map[library_residue->GetName()] = library_residue;
    }
    library_file->SetResidues(residue_map);
    return library_file;
}

void Assembly::BuildStructure(gmml::BuildingStructureOption building_option, vector<string> options, vector<string> file_paths)
{
    switch(building_option)
    {
        case gmml::DISTANCE:
            if(options.size() == 1)
            {
                stringstream ss(description_);
                ss << "Building option: Distance;";
                description_ = ss.str();
                vector<string> tokens = gmml::Split(options.at(0), ":");
                if(tokens.at(0).compare("cutoff") == 0)
                {
                    double cutoff = gmml::ConvertString<double>(tokens.at(1));
                    this->BuildStructureByDistance(cutoff);
                }
                if(tokens.at(0).compare("model_index") == 0)
                {
                    double cutoff = gmml::dCutOff;
                    int model_index = gmml::ConvertString<int>(tokens.at(1));
                    this->BuildStructureByDistance(cutoff, model_index);
                }
            }
            if(options.size() == 2)
            {
                stringstream ss(description_);
                ss << "Building option: Distance;";
                description_ = ss.str();
                vector<string> cutoff_tokens = gmml::Split(options.at(0), ":");
                vector<string> model_tokens = gmml::Split(options.at(1), ":");
                double cutoff = gmml::dCutOff;
                int model_index = 0;
                if(cutoff_tokens.at(0).compare("cutoff") == 0)
                {
                    cutoff = gmml::ConvertString<double>(cutoff_tokens.at(1));
                }
                if(model_tokens.at(0).compare("model_index") == 0)
                {
                    model_index = gmml::ConvertString<int>(model_tokens.at(1));
                }
                BuildStructureByDistance(cutoff, model_index);
            }
            else
            {
                BuildStructureByDistance();
            }
            break;
        case gmml::ORIGINAL:
        {
            stringstream ss(description_);
            ss << "Building option: Original;";
            ss << "File type: " << gmml::ConvertAssemblyInputFileType2String(this->GetSourceFileType()) << ";"
               << "File path: " << this->GetSourceFile() << ";";
            description_ = ss.str();
            this->BuildStructureByOriginalFileBondingInformation();
            break;
        }
        case gmml::DATABASE:
            vector<gmml::InputFileType> types = vector<gmml::InputFileType>();
            for(unsigned int i = 0; i < options.size(); i++)
            {
                vector<string> tokens = gmml::Split(options.at(i), ":");
                if(tokens.at(0).compare("type") == 0)
                {
                    gmml::InputFileType type = gmml::ConvertString2AssemblyInputFileType(tokens.at(1));
                    types.push_back(type);
                }
            }
            if(types.size() == file_paths.size())
            {
                stringstream ss(description_);
                ss << "Building option: Distance;";
                for(unsigned int i = 0; i < types.size(); i++)
                {
                    ss << "File type: " << types.at(i) << ";" << "File path: " << file_paths.at(i) << ";";
                }
                description_ = ss.str();
                this->BuildStructureByDatabaseFilesBondingInformation(types, file_paths);
            }
            break;
    }
}
pthread_mutex_t mutex1 = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t mutex2 = PTHREAD_MUTEX_INITIALIZER;
///7*7 half matrix 6^2/2 + 4, 2 types of thread complete tile, half tile _ 6*6 complete tiles, 22 threads
///or first chunk small next one bigger and so on

void* BuildStructureByDistanceThread(void* args){
    DistanceCalculationThreadArgument* arg = (DistanceCalculationThreadArgument*)args;
    double cutoff = arg->cutoff;
    int model_index = arg->model_index;
    int ti = arg->thread_index;
    int t = arg->number_of_threads;

    //    cout << "Thread" << ti << " start" << endl;
    Assembly::AtomVector all_atoms_of_assembly = arg->a->GetAllAtomsOfAssembly();
    int atoms_size = all_atoms_of_assembly.size();
    int i = ti * (atoms_size/t);

    for(Assembly::AtomVector::iterator it = all_atoms_of_assembly.begin() + ti * (atoms_size/t); ; it++)
    {
        if(ti + 1 < t)
        {
            if(it == all_atoms_of_assembly.begin() + (ti+1) * (atoms_size/t))
                break;
        }
        else
        {
            if(it == all_atoms_of_assembly.end() - 1)
            {
                if((*it)->GetNode() == NULL)
                {
                    Atom* atom = (*it);
                    AtomNode* atom_node = new AtomNode();
                    atom_node->SetAtom(atom);
                    atom->SetNode(atom_node);
                }
                break;
            }
        }
        Atom* atom = (*it);
        AtomNode* atom_node;
        pthread_mutex_lock(&mutex1);
        if(atom->GetNode() == NULL)
        {
            atom_node = new AtomNode();
            atom_node->SetAtom(atom);
            atom->SetNode(atom_node);
        }
        else
            atom_node = atom->GetNode();
        atom_node->SetId(i);
        //        cout << "Thread" << ti << " atom id " << i << endl;
        i++;
        pthread_mutex_unlock(&mutex1);
        for(Assembly::AtomVector::iterator it1 = it + 1; it1 != all_atoms_of_assembly.end(); it1++)
        {
            Atom* neighbor_atom = (*it1);
            // X distance
            if(atom->GetCoordinates().at(model_index)->GetX() - neighbor_atom->GetCoordinates().at(model_index)->GetX() < cutoff)
            {
                // Y distance
                if(atom->GetCoordinates().at(model_index)->GetY() - neighbor_atom->GetCoordinates().at(model_index)->GetY() < cutoff)
                {
                    // Z distance
                    if(atom->GetCoordinates().at(model_index)->GetZ() - neighbor_atom->GetCoordinates().at(model_index)->GetZ() < cutoff)
                    {
                        if((atom->GetCoordinates().at(model_index)->Distance(*(neighbor_atom->GetCoordinates().at(model_index)))) < cutoff)
                        {
                            AtomNode* neighbor_node;
                            pthread_mutex_lock(&mutex1);
                            if (neighbor_atom->GetNode() == NULL)
                            {
                                neighbor_node = new AtomNode();
                                neighbor_node->SetAtom(neighbor_atom);
                            }
                            else
                                neighbor_node = neighbor_atom->GetNode();
                            atom_node->AddNodeNeighbor(neighbor_atom);
                            neighbor_node->AddNodeNeighbor(atom);
                            neighbor_atom->SetNode(neighbor_node);
                            pthread_mutex_unlock(&mutex1);
                        }
                    }
                }
            }
        }
        atom->SetNode(atom_node);
    }

    //    cout << "Thread" << ti << " END" << endl;
    pthread_exit((void*) ti);
}

void* BuildStructureByDistanceByOptimizedThread(void* args){

    DistanceCalculationThreadArgument* arg = (DistanceCalculationThreadArgument*)args;
    double cutoff = arg->cutoff;
    int model_index = arg->model_index;
    int ti = arg->thread_index;
    int t = arg->number_of_threads;


    //    cout << "Thread" << ti << " start" << endl;
    Assembly::AtomVector all_atoms_of_assembly = arg->a->GetAllAtomsOfAssembly();
    int atoms_size = all_atoms_of_assembly.size();

    int increase_factor = 0;
    increase_factor = atoms_size / ((t * t+1) / 2);
    int end_index = 0;
    int begin_index = 0;
    int i = 1;
    for(i = 1; i < ti+1; i++)
    {
        begin_index = begin_index + (i*increase_factor);
    }
    end_index = begin_index + (i*increase_factor);

    for(Assembly::AtomVector::iterator it = all_atoms_of_assembly.begin() + begin_index; ; it++)
    {
        if(ti + 1 < t)///if it's not the last thread
        {
            if(it == all_atoms_of_assembly.begin() + end_index)
                break;
        }
        else
        {
            if(it == all_atoms_of_assembly.end() - 1)
            {
                if((*it)->GetNode() == NULL)
                {
                    Atom* atom = (*it);
                    AtomNode* atom_node = new AtomNode();
                    atom_node->SetAtom(atom);
                    atom->SetNode(atom_node);
                }
                break;
            }
        }
        int index = distance(all_atoms_of_assembly.begin(), it);
        Atom* atom = (*it);
        AtomNode* atom_node;
        pthread_mutex_lock(&mutex1);
        if(atom->GetNode() == NULL)
        {
            atom_node = new AtomNode();
            atom_node->SetAtom(atom);
            atom->SetNode(atom_node);
        }
        else
            atom_node = atom->GetNode();
        atom_node->SetId(index);
        //        cout << "Thread" << ti << " atom id " << i << endl;
        pthread_mutex_unlock(&mutex1);
        for(Assembly::AtomVector::iterator it1 = it + 1; it1 != all_atoms_of_assembly.end(); it1++)
        {
            Atom* neighbor_atom = (*it1);
            // X distance
            if(atom->GetCoordinates().at(model_index)->GetX() - neighbor_atom->GetCoordinates().at(model_index)->GetX() < cutoff)
            {
                // Y distance
                if(atom->GetCoordinates().at(model_index)->GetY() - neighbor_atom->GetCoordinates().at(model_index)->GetY() < cutoff)
                {
                    // Z distance
                    if(atom->GetCoordinates().at(model_index)->GetZ() - neighbor_atom->GetCoordinates().at(model_index)->GetZ() < cutoff)
                    {
                        if((atom->GetCoordinates().at(model_index)->Distance(*(neighbor_atom->GetCoordinates().at(model_index)))) < cutoff)
                        {
                            AtomNode* neighbor_node;
                            pthread_mutex_lock(&mutex1);
                            if (neighbor_atom->GetNode() == NULL)
                            {
                                neighbor_node = new AtomNode();
                                neighbor_node->SetAtom(neighbor_atom);
                            }
                            else
                                neighbor_node = neighbor_atom->GetNode();
                            atom_node->AddNodeNeighbor(neighbor_atom);
                            neighbor_node->AddNodeNeighbor(atom);
                            neighbor_atom->SetNode(neighbor_node);
                            pthread_mutex_unlock(&mutex1);
                        }
                    }
                }
            }
        }
        atom->SetNode(atom_node);
    }

    //    cout << "Thread" << ti << " END" << endl;
    pthread_exit((void*) ti);
}

void* BuildStructureByDistanceByMatrixThread(void* args){

    DistanceCalculationByMatrixThreadArgument* arg = (DistanceCalculationByMatrixThreadArgument*)args;
    int ti = arg->thread_index;
    double cutoff = arg->cutoff;
    int model_index = arg->model_index;
    Assembly::AtomVector* first_chunk = arg->first_chunk;
    Assembly::AtomVector* second_chunk = arg->second_chunk;

    for(Assembly::AtomVector::iterator it = first_chunk->begin(); it != first_chunk->end(); it++)
    {
        int index = distance(first_chunk->begin(), it);
        Atom* atom = (*it);
        AtomNode* atom_node;
        pthread_mutex_lock(&mutex1);
        if(atom->GetNode() == NULL)
        {
            atom_node = new AtomNode();
            atom_node->SetAtom(atom);
            atom->SetNode(atom_node);
        }
        else
            atom_node = atom->GetNode();
        atom_node->SetId(index);
        pthread_mutex_unlock(&mutex1);
        for(Assembly::AtomVector::iterator it1 = second_chunk->begin(); it1 != second_chunk->end(); it1++)
        {
            Atom* neighbor_atom = (*it1);
            // X distance
            if(atom->GetCoordinates().at(model_index)->GetX() - neighbor_atom->GetCoordinates().at(model_index)->GetX() < cutoff)
            {
                // Y distance
                if(atom->GetCoordinates().at(model_index)->GetY() - neighbor_atom->GetCoordinates().at(model_index)->GetY() < cutoff)
                {
                    // Z distance
                    if(atom->GetCoordinates().at(model_index)->GetZ() - neighbor_atom->GetCoordinates().at(model_index)->GetZ() < cutoff)
                    {
                        if((atom->GetCoordinates().at(model_index)->Distance(*(neighbor_atom->GetCoordinates().at(model_index)))) < cutoff)
                        {
                            AtomNode* neighbor_node;
                            pthread_mutex_lock(&mutex1);
                            if (neighbor_atom->GetNode() == NULL)
                            {
                                neighbor_node = new AtomNode();
                                neighbor_node->SetAtom(neighbor_atom);
                            }
                            else
                                neighbor_node = neighbor_atom->GetNode();
                            atom_node->AddNodeNeighbor(neighbor_atom);
                            neighbor_node->AddNodeNeighbor(atom);
                            neighbor_atom->SetNode(neighbor_node);
                            pthread_mutex_unlock(&mutex1);
                        }
                    }
                }
            }
        }
        atom->SetNode(atom_node);
    }
    pthread_exit((void*) ti);
}
void* BuildStructureByDistanceByMatrixDiameterThread(void* args){

    DistanceCalculationByMatrixThreadArgument* arg = (DistanceCalculationByMatrixThreadArgument*)args;
    int ti = arg->thread_index;
    double cutoff = arg->cutoff;
    int model_index = arg->model_index;
    Assembly::AtomVector* first_chunk = arg->first_chunk;
    Assembly::AtomVector* second_chunk = arg->second_chunk;
    int j = ti * 10;

    vector<Assembly::AtomVector*> chunks = vector<Assembly::AtomVector*>();
    chunks.push_back(first_chunk);
    chunks.push_back(second_chunk);
    for(int i = 0; i < 2; i++)
    {
        Assembly::AtomVector* chunk = chunks.at(i);
        if(chunk->size() != 0)
        {
            for(Assembly::AtomVector::iterator it = chunk->begin(); it != chunk->end(); it++)
            {
                Atom* atom = (*it);
                //                    cout << "chunk" << i << " start " << atom->GetId() << endl;
                AtomNode* atom_node;
                pthread_mutex_lock(&mutex1);
                if(atom->GetNode() == NULL)
                {
                    atom_node = new AtomNode();
                    atom_node->SetAtom(atom);
                    atom->SetNode(atom_node);
                }
                else
                    atom_node = atom->GetNode();
                atom_node->SetId(j);
                j++;
                if(it != chunk->end())
                {
                    pthread_mutex_unlock(&mutex1);
                    for(Assembly::AtomVector::iterator it1 = it+1; it1 != chunk->end(); it1++)
                    {
                        Atom* neighbor_atom = (*it1);
                        // X distance
                        if(atom->GetCoordinates().at(model_index)->GetX() - neighbor_atom->GetCoordinates().at(model_index)->GetX() < cutoff)
                        {
                            // Y distance
                            if(atom->GetCoordinates().at(model_index)->GetY() - neighbor_atom->GetCoordinates().at(model_index)->GetY() < cutoff)
                            {
                                // Z distance
                                if(atom->GetCoordinates().at(model_index)->GetZ() - neighbor_atom->GetCoordinates().at(model_index)->GetZ() < cutoff)
                                {
                                    if((atom->GetCoordinates().at(model_index)->Distance(*(neighbor_atom->GetCoordinates().at(model_index)))) < cutoff)
                                    {
                                        AtomNode* neighbor_node;
                                        pthread_mutex_lock(&mutex1);
                                        if (neighbor_atom->GetNode() == NULL)
                                        {
                                            neighbor_node = new AtomNode();
                                            neighbor_node->SetAtom(neighbor_atom);
                                        }
                                        else
                                            neighbor_node = neighbor_atom->GetNode();
                                        atom_node->AddNodeNeighbor(neighbor_atom);
                                        neighbor_node->AddNodeNeighbor(atom);
                                        neighbor_atom->SetNode(neighbor_node);
                                        pthread_mutex_unlock(&mutex1);
                                    }
                                }
                            }
                        }
                    }
                }
                atom->SetNode(atom_node);
            }
        }
    }
    pthread_exit((void*) ti);
}

// MATRIX VERSION
/*
void Assembly::BuildStructureByDistance(int number_of_threads, double cutoff, int model_index)
{
    cout << "Building structure by distance ..." << endl;
    gmml::log(__LINE__, __FILE__, gmml::INF, "Building structure by distance ...");
    model_index_ = model_index;

    Assembly::AtomVector all_atoms_of_assembly = this->GetAllAtomsOfAssembly();
    int atoms_size = all_atoms_of_assembly.size();

    number_of_threads = 5;
    int matrix_size = 3;
//    number_of_threads = 2;
//    int matrix_size = 2;
    pthread_t threads[number_of_threads];
    DistanceCalculationByMatrixThreadArgument arg[number_of_threads];
    int k = 0;
    for(int i = 0; i < matrix_size; i++)
    {
        AtomVector* first_chunk = new AtomVector();
        for(AtomVector::iterator it = all_atoms_of_assembly.begin() + i * (atoms_size/matrix_size); ; it++)
        {
            if(i + 1 < matrix_size)
            {
                if(it == all_atoms_of_assembly.begin() + (i+1) * (atoms_size/matrix_size))
                    break;
            }
            else
            {
                if(it == all_atoms_of_assembly.end())
                    break;
            }
            Atom* atom = (*it);
            first_chunk->push_back(atom);
        }
        for(int j = i+1; j < matrix_size; j++)
        {
            AtomVector* second_chunk = new AtomVector();
            for(AtomVector::iterator it = all_atoms_of_assembly.begin() + (j * (atoms_size/matrix_size)); ; it++)
            {
                if(j + 1 < matrix_size)
                {
                    if(it == all_atoms_of_assembly.begin() + ((j+1)*(atoms_size/matrix_size)) )
                        break;
                }
                else
                {
                    if(it == all_atoms_of_assembly.end())
                        break;
                }
                Atom* atom = (*it);
                second_chunk->push_back(atom);
            }
//            cout << "thread=" << k << ", first chunk size = " << first_chunk->size() << ", starting from " << i * (atoms_size/matrix_size);
//            cout << ", 2nd chunk size = " << second_chunk->size() << ", starting from " << (j * ((atoms_size/matrix_size) )) << endl;
            arg[k] = DistanceCalculationByMatrixThreadArgument(k, model_index, cutoff, first_chunk, second_chunk);
            pthread_create(&threads[k], NULL, &BuildStructureByDistanceByMatrixThread, &arg[k]);
            k++;
        }
    }

    bool no_atoms_for_second_chunk = false;
    for(int i = 0; i < matrix_size; i++)///DIAMETER threads
    {
        AtomVector* first_chunk = new AtomVector();
        AtomVector* second_chunk = new AtomVector();
        for(AtomVector::iterator it = all_atoms_of_assembly.begin() + (i * (atoms_size/matrix_size)); ; it++)
        {
            if(i + 1 < 4)
            {
                if(it == all_atoms_of_assembly.begin() + ((i+1) * (atoms_size/matrix_size)))
                    break;
            }
            else
            {
                if(it == all_atoms_of_assembly.end())
                {
                    no_atoms_for_second_chunk = true;
                    break;
                }
            }
            Atom* atom = (*it);
            first_chunk->push_back(atom);
        }
//        cout << "thread=" << k << ", first chunk size = " << first_chunk->size() << ", starting from " << i * (atoms_size/matrix_size);
        i++;
        if(!no_atoms_for_second_chunk)
        {
            for(AtomVector::iterator it = all_atoms_of_assembly.begin() + (i * (atoms_size/matrix_size)); ; it++)
            {
                if(i + 1 < matrix_size)
                {
                    if(it == all_atoms_of_assembly.begin() + (i+1) * (atoms_size/matrix_size))
                        break;
                }
                else
                    break;

                Atom* atom = (*it);
                second_chunk->push_back(atom);
            }
        }
        //cout << ", 2nd chunk size = " << second_chunk->size() << ", starting from " << (i * (atoms_size/matrix_size)) << endl;
        arg[k] = DistanceCalculationByMatrixThreadArgument(k, model_index, cutoff, first_chunk, second_chunk);
            pthread_create(&threads[k], NULL, &BuildStructureByDistanceByMatrixDiameterThread, &arg[k]);
        k++;
    }
    for(int i = 0; i < number_of_threads; i++)
    {
        pthread_join(threads[i], NULL);
    }
}
*/

///First and second version
void Assembly::BuildStructureByDistance(int number_of_threads, double cutoff, int model_index)
{
    cout << "Building structure by distance ..." << endl;
    gmml::log(__LINE__, __FILE__, gmml::INF, "Building structure by distance ...");
    model_index_ = model_index;

    pthread_t threads[number_of_threads];
    DistanceCalculationThreadArgument arg[number_of_threads];
    for(int i = 0; i < number_of_threads; i++)
    {
        arg[i] = DistanceCalculationThreadArgument(i, number_of_threads, model_index, cutoff, this);
        //        pthread_create(&threads[i], NULL, &BuildStructureByDistanceThread, &arg[i]); ///First version. Workload of threads are not equal
        pthread_create(&threads[i], NULL, &BuildStructureByDistanceByOptimizedThread, &arg[i]); ///Second version. Workload of threads are roughly equal.
    }
    for(int i = 0; i < number_of_threads; i++)
    {
        pthread_join(threads[i], NULL);
    }
}

void Assembly::BuildStructureByOriginalFileBondingInformation()
{
    gmml::InputFileType type = this->GetSourceFileType();
    switch(type)
    {
        case gmml::PDB:
            this->BuildStructureByPDBFileInformation();
            break;
        case gmml::PDBQT:
            break;
        case gmml::TOP:
            this->BuildStructureByTOPFileInformation();
            break;
        case gmml::LIB:
            this->BuildStructureByLIBFileInformation();
            break;
        case gmml::PREP:
            this->BuildStructureByPrepFileInformation();
            break;
        case gmml::TOP_CRD:
            break;
        case gmml::MULTIPLE:
            break;
        case gmml::UNKNOWN:
            break;
    }
}

void Assembly::BuildStructureByPDBFileInformation()
{
    try{
        cout << "Building structure by pdb file information ..." << endl;
        gmml::log(__LINE__, __FILE__, gmml::INF, "Building structure by pdb file information ...");
        PdbFile* pdb_file = new PdbFile(this->GetSourceFile());
        AtomVector all_atoms_of_assembly = this->GetAllAtomsOfAssembly();
        int i = 0;
        for(AtomVector::iterator it = all_atoms_of_assembly.begin(); it != all_atoms_of_assembly.end(); it++)
        {
            Atom* atom = (*it);
            AtomNode* atom_node = new AtomNode();
            atom_node->SetAtom(atom);
            atom_node->SetId(i);
            i++;
            PdbAtom* pdb_atom = pdb_file->GetAtomOfResidueByAtomKey(atom->GetId());
            if(pdb_atom != NULL)
            {
                int atom_serial_number = pdb_atom->GetAtomSerialNumber();
                PdbConnectCard* connectivities = pdb_file->GetConnectivities();
                PdbConnectCard::BondedAtomsSerialNumbersMap bonded_atoms_map = connectivities->GetBondedAtomsSerialNumbers();
                vector<int> bonded_atoms_serial_number = bonded_atoms_map[atom_serial_number];
                for(vector<int>::iterator it1 = bonded_atoms_serial_number.begin(); it1 != bonded_atoms_serial_number.end(); it1++)
                {
                    int bonded_atom_serial_number = *it1;
                    PdbAtom* pdb_bonded_atom = pdb_file->GetAtomBySerialNumber(bonded_atom_serial_number);
                    stringstream sss;
                    sss << pdb_bonded_atom->GetAtomName() << "_" << pdb_bonded_atom->GetAtomSerialNumber() << "_" << pdb_bonded_atom->GetAtomResidueName()
                        << "_" << pdb_bonded_atom->GetAtomChainId() << "_" << pdb_bonded_atom->GetAtomResidueSequenceNumber()
                        << "_" << pdb_bonded_atom->GetAtomInsertionCode() << "_" << pdb_bonded_atom->GetAtomAlternateLocation();
                    string pdb_bonded_atom_key = sss.str();
                    for(AtomVector::iterator it2 = all_atoms_of_assembly.begin(); it2 != all_atoms_of_assembly.end(); it2++)
                    {
                        if(it != it2)
                        {
                            Atom* assembly_atom = (*it2);
                            string assembly_atom_key = assembly_atom->GetId();
                            if(assembly_atom_key.compare(pdb_bonded_atom_key) == 0)
                            {
                                atom_node->AddNodeNeighbor(assembly_atom);
                                break;
                            }
                        }
                    }
                }
            }
            atom->SetNode(atom_node);
        }
    }
    catch(PdbFileSpace::PdbFileProcessingException &ex)
    {}
}

void Assembly::BuildStructureByTOPFileInformation()
{
    cout << "Building structure by topology file information ..." << endl;
    gmml::log(__LINE__, __FILE__, gmml::INF, "Building structure by topology file information ...");
    TopologyFile* topology_file = new TopologyFile(gmml::Split(this->GetSourceFile(), ";")[0]);
    AtomVector all_atoms_of_assembly = this->GetAllAtomsOfAssembly();
    int i = 0;
    for (AtomVector::iterator it = all_atoms_of_assembly.begin(); it != all_atoms_of_assembly.end(); it++)
    {
        Atom* atom_1 = (*it);
        AtomNode* atom_node = new AtomNode();
        atom_node->SetAtom(atom_1);
        atom_node->SetId(i);
        i++;
        for(AtomVector::iterator it1 = all_atoms_of_assembly.begin(); it1 != all_atoms_of_assembly.end(); it1++)
        {
            if(it != it1)
            {
                Atom* atom_2 = (*it1);
                stringstream ss;
                ss << gmml::Split(atom_1->GetId(), "_").at(2) << "(" << gmml::Split(atom_1->GetId(), "_").at(4) << ")"
                   << ":" << gmml::Split(atom_1->GetId(), "_").at(0) << "(" <<  gmml::Split(atom_1->GetId(), "_").at(1) << ")" << "-"
                   << gmml::Split(atom_2->GetId(), "_").at(2) << "(" << gmml::Split(atom_2->GetId(), "_").at(4) << ")"
                   << ":" << gmml::Split(atom_2->GetId(), "_").at(0) << "(" << gmml::Split(atom_2->GetId(), "_").at(1) << ")";
                string key = ss.str();
                TopologyFile::TopologyBondMap topology_bond = topology_file->GetBonds();
                for(TopologyFile::TopologyBondMap::iterator it2 = topology_bond.begin(); it2 != topology_bond.end(); it2++)
                {
                    TopologyBond* bond = (*it2).second;
                    stringstream sss;
                    sss << bond->GetResidueNames().at(0) << ":" << bond->GetBonds().at(0) << "-" << bond->GetResidueNames().at(1) << ":" << bond->GetBonds().at(1);
                    string topology_bond_key = sss.str();
                    if(key.compare(topology_bond_key) == 0)
                    {
                        atom_node->AddNodeNeighbor(atom_2);
                        break;
                    }
                    stringstream ssss;
                    ssss << bond->GetResidueNames().at(1) << ":" << bond->GetBonds().at(1) << "-" << bond->GetResidueNames().at(0) << ":" << bond->GetBonds().at(0);
                    topology_bond_key = ssss.str();
                    if(key.compare(topology_bond_key) == 0)
                    {
                        atom_node->AddNodeNeighbor(atom_2);
                        break;
                    }
                }
            }
        }
        atom_1->SetNode(atom_node);
    }
}

void Assembly::BuildStructureByLIBFileInformation()
{
    cout << "Building structure by library file information..." << endl;
    gmml::log(__LINE__, __FILE__, gmml::INF, "Building structure by library file information ...");
    LibraryFile* library_file = new LibraryFile(this->GetSourceFile());
    AtomVector all_atoms_of_assembly = this->GetAllAtomsOfAssembly();
    int i = 0;
    for(AtomVector::iterator it = all_atoms_of_assembly.begin(); it != all_atoms_of_assembly.end(); it++)
    {
        Atom* atom = (*it);
        AtomNode* atom_node = new AtomNode();
        atom_node->SetAtom(atom);
        atom_node->SetId(i);
        i++;
        Residue* assembly_residue = atom->GetResidue();
        vector<string> atom_id_tokens = Split(atom->GetId(), "_");
        LibraryFileResidue* library_residue = library_file->GetLibraryResidueByResidueName(assembly_residue->GetName());
        if(library_residue != NULL)
        {
            LibraryFileAtom* library_atom = library_residue->GetAtomByOrder(ConvertString<int>(atom_id_tokens.at(1)));
            if(library_atom != NULL)
            {
                vector<int> library_bonded_atom_indices = library_atom->GetBondedAtomsIndices();
                for(vector<int>::iterator it1 = library_bonded_atom_indices.begin(); it1 != library_bonded_atom_indices.end(); it1++)
                {
                    int library_bonded_atom_index = (*it1);
                    LibraryFileAtom* library_atom = library_residue->GetAtomByOrder(library_bonded_atom_index);
                    for(AtomVector::iterator it2 = all_atoms_of_assembly.begin(); it2 != all_atoms_of_assembly.end(); it2++)
                    {
                        Atom* assembly_atom = (*it2);
                        string assembly_atom_id = assembly_atom->GetId();
                        stringstream ss;
                        ss << library_residue->GetName() << ":" << library_atom->GetName() << "(" << library_atom->GetAtomOrder() << ")";
                        string library_atom_id = ss.str();
                        vector<string> assembly_atom_id_tokens = gmml::Split(assembly_atom_id, "_");
                        stringstream sss;
                        sss << assembly_atom_id_tokens.at(2) << ":" << assembly_atom_id_tokens.at(0) << "(" << assembly_atom_id_tokens.at(1) << ")";
                        if(sss.str().compare(library_atom_id) == 0)
                        {
                            atom_node->AddNodeNeighbor(assembly_atom);
                            break;
                        }
                    }
                }
            }
        }
        atom->SetNode(atom_node);
    }
}

void Assembly::BuildStructureByPrepFileInformation()
{
    cout << "Building structure by prep file information ..." << endl;
    gmml::log(__LINE__, __FILE__, gmml::INF, "Building structure by prep file information ...");
    PrepFile* prep_file = new PrepFile(this->GetSourceFile());
    AtomVector all_atoms_of_assembly = this->GetAllAtomsOfAssembly();
    int i = 0;
    for(AtomVector::iterator it = all_atoms_of_assembly.begin(); it != all_atoms_of_assembly.end(); it++)
    {
        Atom* atom = (*it);
        AtomNode* atom_node = new AtomNode();
        atom_node->SetAtom(atom);
        atom_node->SetId(i);
        i++;
        Residue* assembly_residue = atom->GetResidue();
        //        cout << assembly_residue->GetName() << endl;
        PrepFileResidue* prep_residue = prep_file->GetResidues()[assembly_residue->GetName()];
        if(prep_residue != NULL)
        {
            PrepFileAtom* prep_atom = prep_residue->GetPrepAtomByName(atom->GetName());
            if(prep_atom != NULL)
            {
                vector<int> bonded_atoms_index = prep_residue->GetBondingsOfResidue()[prep_residue->GetAtomIndexByName(atom->GetName())];
                for(vector<int>::iterator it1 = bonded_atoms_index.begin(); it1 != bonded_atoms_index.end(); it1++)
                {
                    int bonded_atom_index = (*it1);
                    PrepFileAtom* bonded_atom = prep_residue->GetPrepAtomByName(prep_residue->GetAtomNameByIndex(bonded_atom_index));
                    stringstream ss;
                    ss << prep_residue->GetName() << ":" << bonded_atom->GetName();
                    for(AtomVector::iterator it2 = all_atoms_of_assembly.begin(); it2 != all_atoms_of_assembly.end(); it2++)
                    {
                        Atom* assembly_atom = (*it2);
                        vector<string> atom_id_tokens = gmml::Split(assembly_atom->GetId(), "_");
                        stringstream sss;
                        sss << atom_id_tokens.at(2) << ":" << atom_id_tokens.at(0);
                        if(sss.str().compare(ss.str()) == 0)
                        {
                            atom_node->AddNodeNeighbor(assembly_atom);
                            break;
                        }
                    }
                }
            }
        }
        atom->SetNode(atom_node);
    }
}

void Assembly::BuildStructureByDatabaseFilesBondingInformation(vector<gmml::InputFileType> types, vector<string> file_paths)
{
    cout << "Building structure by dataset files information ..." << endl;
    gmml::log(__LINE__, __FILE__, gmml::INF, "Building structure by dataset files information ...");
    AtomVector all_atoms_of_assembly = this->GetAllAtomsOfAssembly();
    int i = 0;
    for(AtomVector::iterator it = all_atoms_of_assembly.begin(); it != all_atoms_of_assembly.end(); it++)
    {
        Atom* atom = (*it);
        AtomNode* atom_node = new AtomNode();
        atom_node->SetAtom(atom);
        atom_node->SetId(i);
        i++;
        Residue* assembly_residue = atom->GetResidue();
        for(unsigned int i = 0; i < types.size(); i++)
        {
            if(types.at(i) == gmml::LIB)
            {
                string lib_path = file_paths.at(i);
                LibraryFile* library_file = new LibraryFile(lib_path);
                LibraryFileResidue* library_residue = library_file->GetLibraryResidueByResidueName(assembly_residue->GetName());
                if(library_residue != NULL)
                {
                    LibraryFileAtom* library_atom = library_residue->GetLibraryAtomByAtomName(atom->GetName());
                    if(library_atom != NULL)
                    {
                        vector<int> library_bonded_atom_indices = library_atom->GetBondedAtomsIndices();
                        for(vector<int>::iterator it1 = library_bonded_atom_indices.begin(); it1 != library_bonded_atom_indices.end(); it1++)
                        {
                            int library_bonded_atom_index = (*it1);
                            LibraryFileAtom* library_atom = library_residue->GetAtomByIndex(library_bonded_atom_index);
                            for(AtomVector::iterator it2 = all_atoms_of_assembly.begin(); it2 != all_atoms_of_assembly.end(); it2++)
                            {
                                Atom* assembly_atom = (*it2);
                                string assembly_atom_id = assembly_atom->GetId();
                                stringstream ss;
                                ss << library_residue->GetName() << ":" << library_atom->GetName();
                                string library_atom_id = ss.str();
                                if(assembly_atom_id.compare(library_atom_id) == 0)
                                {
                                    atom_node->AddNodeNeighbor(assembly_atom);
                                    break;
                                }
                            }
                        }
                    }
                }
            }
            if(types.at(i) == gmml::PREP)
            {
                string prep_path = file_paths.at(i);
                PrepFile* prep_file = new PrepFile(prep_path);
                PrepFileResidue* prep_residue = prep_file->GetResidues()[assembly_residue->GetName()];
                if(prep_residue != NULL)
                {
                    PrepFileAtom* prep_atom = prep_residue->GetPrepAtomByName(atom->GetName());
                    if(prep_atom != NULL)
                    {
                        vector<int> bonded_atoms_index = prep_residue->GetBondingsOfResidue()[prep_residue->GetAtomIndexByName(atom->GetName())];
                        for(vector<int>::iterator it1 = bonded_atoms_index.begin(); it1 != bonded_atoms_index.end(); it1++)
                        {
                            int bonded_atom_index = (*it1);
                            PrepFileAtom* bonded_atom = prep_residue->GetPrepAtomByName(prep_residue->GetAtomNameByIndex(bonded_atom_index));
                            stringstream ss;
                            ss << prep_residue->GetName() << ":" << bonded_atom->GetName();
                            for(AtomVector::iterator it2 = all_atoms_of_assembly.begin(); it2 != all_atoms_of_assembly.end(); it2++)
                            {
                                Atom* assembly_atom = (*it2);
                                if(assembly_atom->GetId().compare(ss.str()) == 0)
                                {
                                    atom_node->AddNodeNeighbor(assembly_atom);
                                    break;
                                }
                            }
                        }
                    }
                }

            }
            atom->SetNode(atom_node);
        }
    }
}

int Assembly::CountNumberOfAtoms()
{
    int counter = 0;
    for(AssemblyVector::iterator it = assemblies_.begin(); it != assemblies_.end(); it++)
    {
        Assembly* assembly = (*it);
        counter += assembly->CountNumberOfAtoms();
    }
    for(ResidueVector:: iterator it1 = residues_.begin(); it1 != residues_.end(); it1++)
    {
        Residue* residue = (*it1);
        Residue::AtomVector atoms = residue->GetAtoms();
        counter += atoms.size();
    }
    return counter;
}

int Assembly::CountNumberOfAtomTypes()
{
    vector<string> type_list = vector<string>();
    AtomVector atoms = GetAllAtomsOfAssembly();
    for(AtomVector::iterator it = atoms.begin(); it != atoms.end(); it++)
    {
        Atom* atom = (*it);
        string atom_type = atom->GetAtomType();
        if(find(type_list.begin(), type_list.end(), atom_type) != type_list.end())
        {}
        else
        {
            type_list.push_back(atom_type);
        }
    }
    return type_list.size();
}

int Assembly::CountNumberOfResidues()
{
    int counter = 0;
    for(AssemblyVector::iterator it = assemblies_.begin(); it != assemblies_.end(); it++)
    {
        Assembly* assembly = (*it);
        counter += assembly->CountNumberOfResidues();
    }
    counter += residues_.size();
    return counter;
}

int Assembly::CountNumberOfBondsIncludingHydrogen(string parameter_file_path)
{
    ParameterFile* parameter_file = new ParameterFile(parameter_file_path);
    ParameterFileSpace::ParameterFile::BondMap bonds = parameter_file->GetBonds();
    AtomVector atoms = GetAllAtomsOfAssembly();
    int counter = 0;
    for(AtomVector::iterator it = atoms.begin(); it != atoms.end(); it++)
    {
        Atom* atom = (*it);
        string atom_name = atom->GetName();
        AtomNode* atom_node = atom->GetNode();
        if(atom_node != NULL)
        {
            AtomVector node_neighbors = atom_node->GetNodeNeighbors();
            //            if((atom_name.substr(0,1).compare("H") == 0 ||
            //                (atom_name.substr(1,1).compare("H") == 0 && isdigit(ConvertString<char>(atom_name.substr(0,1))))))
            //            {
            //                counter += node_neighbors.size();
            //            }
            //            else
            //            {
            for(AtomVector::iterator it1 = node_neighbors.begin(); it1 != node_neighbors.end(); it1++)
            {
                Atom* node_neighbor = (*it1);
                string node_neighbor_name = node_neighbor->GetName();
                if((atom_name.substr(0,1).compare("H") == 0 || (atom_name.substr(1,1).compare("H") == 0 && isdigit(ConvertString<char>(atom_name.substr(0,1))))) ||
                        (node_neighbor_name.substr(0,1).compare("H") == 0 || (node_neighbor_name.substr(1,1).compare("H") == 0 && isdigit(ConvertString<char>(node_neighbor_name.substr(0,1))))))
                {
                    vector<string> atom_pair_type = vector<string>();
                    vector<string> reverse_atom_pair_type = vector<string>();
                    atom_pair_type.push_back(atom->GetAtomType());
                    atom_pair_type.push_back(node_neighbor->GetAtomType());
                    reverse_atom_pair_type.push_back(node_neighbor->GetAtomType());
                    reverse_atom_pair_type.push_back(atom->GetAtomType());

                    if(bonds.find(atom_pair_type) != bonds.end() || bonds.find(reverse_atom_pair_type) != bonds.end())
                        counter++;
                }
            }
            //            }
        }
    }
    return counter/2;
}

int Assembly::CountNumberOfBondsExcludingHydrogen(string parameter_file_path)
{
    ParameterFile* parameter_file = new ParameterFile(parameter_file_path);
    ParameterFileSpace::ParameterFile::BondMap bonds = parameter_file->GetBonds();
    AtomVector atoms = GetAllAtomsOfAssembly();
    int counter = 0;
    for(AtomVector::iterator it = atoms.begin(); it != atoms.end(); it++)
    {
        Atom* atom = (*it);
        string atom_name = atom->GetName();
        AtomNode* atom_node = atom->GetNode();
        if(atom_node != NULL)
        {
            AtomVector node_neighbors = atom_node->GetNodeNeighbors();

            //            if((atom_name.substr(0,1).compare("H") == 0 ||
            //                (atom_name.substr(1,1).compare("H") == 0 && isdigit(ConvertString<char>(atom_name.substr(0,1))))))
            //            {}
            //            else
            //            {
            for(AtomVector::iterator it1 = node_neighbors.begin(); it1 != node_neighbors.end(); it1++)
            {
                Atom* node_neighbor = (*it1);
                string node_neighbor_name = node_neighbor->GetName();
                if((atom_name.substr(0,1).compare("H") == 0 || (atom_name.substr(1,1).compare("H") == 0 && isdigit(ConvertString<char>(atom_name.substr(0,1))))) ||
                        (node_neighbor_name.substr(0,1).compare("H") == 0 || (node_neighbor_name.substr(1,1).compare("H") == 0 && isdigit(ConvertString<char>(node_neighbor_name.substr(0,1))))))
                {}
                else
                {
                    vector<string> atom_pair_type = vector<string>();
                    vector<string> reverse_atom_pair_type = vector<string>();
                    atom_pair_type.push_back(atom->GetAtomType());
                    atom_pair_type.push_back(node_neighbor->GetAtomType());
                    reverse_atom_pair_type.push_back(node_neighbor->GetAtomType());
                    reverse_atom_pair_type.push_back(atom->GetAtomType());

                    if(bonds.find(atom_pair_type) != bonds.end() || bonds.find(reverse_atom_pair_type) != bonds.end())
                        counter++;
                }
            }
            //            }
        }
    }
    return counter/2;
}

int Assembly::CountNumberOfBonds()
{
    AtomVector atoms = GetAllAtomsOfAssembly();
    int counter = 0;
    for(AtomVector::iterator it = atoms.begin(); it != atoms.end(); it++)
    {
        Atom* atom = (*it);
        AtomNode* atom_node = atom->GetNode();
        if(atom_node != NULL)
        {
            AtomVector node_neighbors = atom_node->GetNodeNeighbors();
            counter += node_neighbors.size();
        }
    }
    return counter/2;
}

int Assembly::CountNumberOfBondTypes(string parameter_file_path)
{
    AtomVector atoms = GetAllAtomsOfAssembly();
    vector<string> type_list = vector<string>();
    ParameterFile* parameter_file = new ParameterFile(parameter_file_path);
    ParameterFileSpace::ParameterFile::BondMap bonds = parameter_file->GetBonds();

    for(AtomVector::iterator it = atoms.begin(); it != atoms.end(); it++)
    {
        Atom* atom = (*it);
        string atom_bond_type = atom->GetAtomType();
        AtomNode* atom_node = atom->GetNode();
        if(atom_node != NULL)
        {
            AtomVector node_neighbors = atom_node->GetNodeNeighbors();
            for(AtomVector::iterator it1 = node_neighbors.begin(); it1 != node_neighbors.end(); it1++)
            {
                Atom* node_neighbor = (*it1);
                string node_neighbor_bond_type = node_neighbor->GetAtomType();
                stringstream ss;
                ss << atom_bond_type << "_" << node_neighbor_bond_type;
                string key = ss.str();
                stringstream ss1;
                ss1 << node_neighbor_bond_type << "_" << atom_bond_type;
                string key1 = ss1.str();
                if(find(type_list.begin(), type_list.end(), key ) != type_list.end() ||
                        find(type_list.begin(), type_list.end(), key1 ) != type_list.end())
                {}
                else
                {
                    vector<string> atom_pair_type = vector<string>();
                    vector<string> reverse_atom_pair_type = vector<string>();
                    atom_pair_type.push_back(atom->GetAtomType());
                    atom_pair_type.push_back(node_neighbor->GetAtomType());
                    reverse_atom_pair_type.push_back(node_neighbor->GetAtomType());
                    reverse_atom_pair_type.push_back(atom->GetAtomType());

                    if(bonds.find(atom_pair_type) != bonds.end() || bonds.find(reverse_atom_pair_type) != bonds.end())
                        type_list.push_back(key);
                }
            }
        }
    }
    return type_list.size();
}

int Assembly::CountNumberOfAnglesIncludingHydrogen(string parameter_file_path)
{
    ParameterFile* parameter_file = new ParameterFile(parameter_file_path);
    ParameterFileSpace::ParameterFile::AngleMap angles = parameter_file->GetAngles();
    AtomVector atoms = GetAllAtomsOfAssembly();
    int counter = 0;
    for(AtomVector::iterator it = atoms.begin(); it != atoms.end(); it++)
    {
        Atom* atom = (*it);
        stringstream ss;
        ss << atom->GetId();
        string atom_name = atom->GetName();
        AtomNode* atom_node = atom->GetNode();
        if(atom_node != NULL)
        {
            AtomVector neighbors = atom_node->GetNodeNeighbors();
            for(AtomVector::iterator it1 = neighbors.begin(); it1 != neighbors.end(); it1++)
            {
                Atom* neighbor = (*it1);
                string neighbor_name = neighbor->GetName();
                AtomNode* neighbor_atom_node = neighbor->GetNode();
                AtomVector neighbors_of_neighbor = neighbor_atom_node->GetNodeNeighbors();
                for(AtomVector::iterator it2 = neighbors_of_neighbor.begin(); it2 != neighbors_of_neighbor.end(); it2++)
                {
                    Atom* neighbor_of_neighbor = (*it2);
                    stringstream ss1;
                    ss1 << neighbor_of_neighbor->GetId();
                    if(ss.str().compare(ss1.str()) != 0)
                    {
                        string neighbor_of_neighbor_name = neighbor_of_neighbor->GetName();
                        if((atom_name.substr(0,1).compare("H") == 0 || (atom_name.substr(1,1).compare("H") == 0 && isdigit(ConvertString<char>(atom_name.substr(0,1))))) ||
                                (neighbor_name.substr(0,1).compare("H") == 0 || (neighbor_name.substr(1,1).compare("H") == 0 && isdigit(ConvertString<char>(neighbor_name.substr(0,1))))) ||
                                (neighbor_of_neighbor_name.substr(0,1).compare("H") == 0 || (neighbor_of_neighbor_name.substr(1,1).compare("H") == 0 && isdigit(ConvertString<char>(neighbor_of_neighbor_name.substr(0,1))))))
                        {
                            vector<string> angle_type = vector<string>();
                            vector<string> reverse_angle_type = vector<string>();
                            angle_type.push_back(atom->GetAtomType());
                            angle_type.push_back(neighbor->GetAtomType());
                            angle_type.push_back(neighbor_of_neighbor->GetAtomType());
                            reverse_angle_type.push_back(neighbor_of_neighbor->GetAtomType());
                            reverse_angle_type.push_back(neighbor->GetAtomType());
                            reverse_angle_type.push_back(atom->GetAtomType());
                            if(angles.find(angle_type) != angles.end() || angles.find(reverse_angle_type) != angles.end())
                                counter++;
                        }
                    }
                }
            }
        }
    }
    return counter/2;
}

int Assembly::CountNumberOfAnglesExcludingHydrogen(string parameter_file_path)
{
    ParameterFile* parameter_file = new ParameterFile(parameter_file_path);
    ParameterFileSpace::ParameterFile::AngleMap angles = parameter_file->GetAngles();
    AtomVector atoms = GetAllAtomsOfAssembly();
    int counter = 0;
    for(AtomVector::iterator it = atoms.begin(); it != atoms.end(); it++)
    {
        Atom* atom = (*it);
        stringstream ss;
        ss << atom->GetId();
        string atom_name = atom->GetName();
        AtomNode* atom_node = atom->GetNode();
        if(atom_node != NULL)
        {
            AtomVector neighbors = atom_node->GetNodeNeighbors();
            for(AtomVector::iterator it1 = neighbors.begin(); it1 != neighbors.end(); it1++)
            {
                Atom* neighbor = (*it1);
                string neighbor_name = neighbor->GetName();
                AtomNode* neighbor_atom_node = neighbor->GetNode();
                AtomVector neighbors_of_neighbor = neighbor_atom_node->GetNodeNeighbors();
                for(AtomVector::iterator it2 = neighbors_of_neighbor.begin(); it2 != neighbors_of_neighbor.end(); it2++)
                {
                    Atom* neighbor_of_neighbor = (*it2);
                    stringstream ss1;
                    ss1 << neighbor_of_neighbor->GetId();
                    if(ss.str().compare(ss1.str()) != 0)
                    {
                        string neighbor_of_neighbor_name = neighbor_of_neighbor->GetName();
                        if((atom_name.substr(0,1).compare("H") == 0 || (atom_name.substr(1,1).compare("H") == 0 && isdigit(ConvertString<char>(atom_name.substr(0,1))))) ||
                                (neighbor_name.substr(0,1).compare("H") == 0 || (neighbor_name.substr(1,1).compare("H") == 0 && isdigit(ConvertString<char>(neighbor_name.substr(0,1))))) ||
                                (neighbor_of_neighbor_name.substr(0,1).compare("H") == 0 || (neighbor_of_neighbor_name.substr(1,1).compare("H") == 0 && isdigit(ConvertString<char>(neighbor_of_neighbor_name.substr(0,1))))))
                        {}
                        else
                        {
                            vector<string> angle_type = vector<string>();
                            vector<string> reverse_angle_type = vector<string>();
                            angle_type.push_back(atom->GetAtomType());
                            angle_type.push_back(neighbor->GetAtomType());
                            angle_type.push_back(neighbor_of_neighbor->GetAtomType());
                            reverse_angle_type.push_back(neighbor_of_neighbor->GetAtomType());
                            reverse_angle_type.push_back(neighbor->GetAtomType());
                            reverse_angle_type.push_back(atom->GetAtomType());
                            if(angles.find(angle_type) != angles.end() || angles.find(reverse_angle_type) != angles.end())
                                counter++;
                        }
                    }
                }
            }
        }
    }
    return counter/2;
}

int Assembly::CountNumberOfAngles()
{
    AtomVector atoms = GetAllAtomsOfAssembly();
    int counter = 0;
    for(AtomVector::iterator it = atoms.begin(); it != atoms.end(); it++)
    {
        Atom* atom = (*it);
        stringstream ss;
        ss << atom->GetId();
        AtomNode* atom_node = atom->GetNode();
        if(atom_node != NULL)
        {
            AtomVector neighbors = atom_node->GetNodeNeighbors();
            for(AtomVector::iterator it1 = neighbors.begin(); it1 != neighbors.end(); it1++)
            {
                Atom* neighbor = (*it1);
                AtomNode* neighbor_atom_node = neighbor->GetNode();
                AtomVector neighbors_of_neighbor = neighbor_atom_node->GetNodeNeighbors();
                for(AtomVector::iterator it2 = neighbors_of_neighbor.begin(); it2 != neighbors_of_neighbor.end(); it2++)
                {
                    Atom* neighbor_of_neighbor = (*it2);
                    stringstream ss1;
                    ss1 << neighbor_of_neighbor->GetId();
                    if(ss.str().compare(ss1.str()) != 0)
                    {
                        counter++;
                    }
                }
            }
        }
    }
    return counter/2;
}

int Assembly::CountNumberOfAngleTypes(string parameter_file_path)
{
    ParameterFile* parameter_file = new ParameterFile(parameter_file_path);
    ParameterFileSpace::ParameterFile::AngleMap angles = parameter_file->GetAngles();
    AtomVector atoms = GetAllAtomsOfAssembly();

    vector<string> type_list = vector<string>();
    for(AtomVector::iterator it = atoms.begin(); it != atoms.end(); it++)
    {
        Atom* atom = (*it);
        string type1 = atom->GetAtomType();
        stringstream ss;
        ss << atom->GetId();
        AtomNode* atom_node = atom->GetNode();
        if(atom_node != NULL)
        {
            AtomVector neighbors = atom_node->GetNodeNeighbors();
            for(AtomVector::iterator it1 = neighbors.begin(); it1 != neighbors.end(); it1++)
            {
                Atom* neighbor = (*it1);
                string type2 = neighbor->GetAtomType();
                AtomNode* neighbor_atom_node = neighbor->GetNode();
                AtomVector neighbors_of_neighbor = neighbor_atom_node->GetNodeNeighbors();
                for(AtomVector::iterator it2 = neighbors_of_neighbor.begin(); it2 != neighbors_of_neighbor.end(); it2++)
                {
                    Atom* neighbor_of_neighbor = (*it2);
                    string type3 = neighbor_of_neighbor->GetAtomType();
                    stringstream ss1;
                    ss1 << neighbor_of_neighbor->GetId();
                    if(ss.str().compare(ss1.str()) != 0)
                    {
                        stringstream ss2;
                        ss2 << type1 << "_" << type2 << "_" << type3;
                        stringstream ss3;
                        ss3 << type3 << "_" << type2 << "_" << type1;
                        if(find(type_list.begin(), type_list.end(), ss2.str()) == type_list.end() &&
                                find(type_list.begin(), type_list.end(), ss3.str()) == type_list.end())
                        {
                            vector<string> angle_type = vector<string>();
                            vector<string> reverse_angle_type = vector<string>();
                            angle_type.push_back(atom->GetAtomType());
                            angle_type.push_back(neighbor->GetAtomType());
                            angle_type.push_back(neighbor_of_neighbor->GetAtomType());
                            reverse_angle_type.push_back(neighbor_of_neighbor->GetAtomType());
                            reverse_angle_type.push_back(neighbor->GetAtomType());
                            reverse_angle_type.push_back(atom->GetAtomType());
                            if(angles.find(angle_type) != angles.end() || angles.find(reverse_angle_type) != angles.end())
                            {
                                type_list.push_back(ss2.str());
                            }
                        }
                    }
                }
            }
        }
    }
    return type_list.size();
}

int Assembly::CountNumberOfDihedralsIncludingHydrogen(string parameter_file_path)
{
    ParameterFile* parameter_file = new ParameterFile(parameter_file_path);
    AtomVector atoms = GetAllAtomsOfAssembly();
    int counter = 0;
    int improper_counter = 0;
    //    int not_found_counter = 0;
    for(AtomVector::iterator it = atoms.begin(); it != atoms.end(); it++)
    {
        Atom* atom = (*it);
        stringstream ss;
        ss << atom->GetId();
        string atom_name = atom->GetName();
        AtomNode* atom_node = atom->GetNode();
        if(atom_node != NULL)
        {
            AtomVector neighbors = atom_node->GetNodeNeighbors();
            for(AtomVector::iterator it1 = neighbors.begin(); it1 != neighbors.end(); it1++)
            {
                Atom* neighbor = (*it1);
                stringstream ss1;
                ss1 << neighbor->GetId();
                string neighbor_name = neighbor->GetName();
                AtomNode* neighbor_node = neighbor->GetNode();
                AtomVector neighbors_of_neighbor = neighbor_node->GetNodeNeighbors();
                for(AtomVector::iterator it2 = neighbors_of_neighbor.begin(); it2 != neighbors_of_neighbor.end(); it2++)
                {
                    Atom* neighbor_of_neighbor = (*it2);
                    stringstream ss2;
                    ss2 << neighbor_of_neighbor->GetId();
                    if(ss.str().compare(ss2.str()) != 0)
                    {
                        string neighbor_of_neighbor_name = neighbor_of_neighbor->GetName();
                        AtomNode* neighbor_of_neighbor_node = neighbor_of_neighbor->GetNode();
                        AtomVector neighbors_of_neighbor_of_neighbor = neighbor_of_neighbor_node->GetNodeNeighbors();
                        for(AtomVector::iterator it3 = neighbors_of_neighbor_of_neighbor.begin(); it3 != neighbors_of_neighbor_of_neighbor.end(); it3++)
                        {
                            Atom* neighbor_of_neighbor_of_neighbor = (*it3);
                            stringstream ss3;
                            ss3 << neighbor_of_neighbor_of_neighbor->GetId();
                            if(ss1.str().compare(ss3.str()) != 0)
                            {
                                string neighbor_of_neighbor_of_neighbor_name = neighbor_of_neighbor_of_neighbor->GetName();
                                if((atom_name.substr(0,1).compare("H") == 0 || (atom_name.substr(1,1).compare("H") == 0 && isdigit(ConvertString<char>(atom_name.substr(0,1))))) ||
                                        (neighbor_name.substr(0,1).compare("H") == 0 || (neighbor_name.substr(1,1).compare("H") == 0 && isdigit(ConvertString<char>(neighbor_name.substr(0,1))))) ||
                                        (neighbor_of_neighbor_name.substr(0,1).compare("H") == 0 || (neighbor_of_neighbor_name.substr(1,1).compare("H") == 0 && isdigit(ConvertString<char>(neighbor_of_neighbor_name.substr(0,1))))) ||
                                        (neighbor_of_neighbor_of_neighbor_name.substr(0,1).compare("H") == 0 || (neighbor_of_neighbor_of_neighbor_name.substr(1,1).compare("H") == 0 && isdigit(ConvertString<char>(neighbor_of_neighbor_of_neighbor_name.substr(0,1))))))
                                {

                                    vector<vector<string> > all_atom_type_permutations = CreateAllAtomTypePermutationsforDihedralType(atom->GetAtomType(), neighbor->GetAtomType(),
                                                                                                                                      neighbor_of_neighbor->GetAtomType(), neighbor_of_neighbor_of_neighbor->GetAtomType());
                                    ParameterFile::DihedralMap dihedrals = parameter_file->GetDihedrals();
                                    for(vector<vector<string> >::iterator it4 = all_atom_type_permutations.begin(); it4 != all_atom_type_permutations.end(); it4++)
                                    {
                                        vector<string> atom_types = (*it4);
                                        if(dihedrals[atom_types] != NULL)
                                        {
                                            ParameterFileDihedral* parameter_file_dihedrals = dihedrals[atom_types];
                                            int terms_count = parameter_file_dihedrals->GetTerms().size();
                                            counter += terms_count;
                                            break;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            ///Improper Dihedrals
            if(neighbors.size() == 3)
            {
                Atom* neighbor1 = neighbors.at(0);
                Atom* neighbor2 = neighbors.at(1);
                Atom* neighbor3 = neighbors.at(2);
                string neighbor1_name = neighbor1->GetName();
                string neighbor2_name = neighbor2->GetName();
                string neighbor3_name = neighbor3->GetName();
                vector<vector<string> > all_improper_dihedrals_atom_type_permutations = CreateAllAtomTypePermutationsforImproperDihedralType(neighbor1->GetAtomType(), neighbor2->GetAtomType(),
                                                                                                                                             neighbor3->GetAtomType(), atom->GetAtomType());

                ParameterFile::DihedralMap dihedrals = parameter_file->GetDihedrals();
                for(vector<vector<string> >::iterator it1 = all_improper_dihedrals_atom_type_permutations.begin(); it1 != all_improper_dihedrals_atom_type_permutations.end(); it1++)
                {
                    vector<string> improper_dihedral_permutation = (*it1);
                    if((atom_name.substr(0,1).compare("H") == 0 || (atom_name.substr(1,1).compare("H") == 0 && isdigit(ConvertString<char>(atom_name.substr(0,1))))) ||
                            (neighbor1_name.substr(0,1).compare("H") == 0 || (neighbor1_name.substr(1,1).compare("H") == 0 && isdigit(ConvertString<char>(neighbor1_name.substr(0,1))))) ||
                            (neighbor2_name.substr(0,1).compare("H") == 0 || (neighbor2_name.substr(1,1).compare("H") == 0 && isdigit(ConvertString<char>(neighbor2_name.substr(0,1))))) ||
                            (neighbor3_name.substr(0,1).compare("H") == 0 || (neighbor3_name.substr(1,1).compare("H") == 0 && isdigit(ConvertString<char>(neighbor3_name.substr(0,1))))))
                    {
                        if(dihedrals[improper_dihedral_permutation] != NULL)
                        {
                            ParameterFileDihedral* parameter_file_dihedrals = dihedrals[improper_dihedral_permutation];
                            int terms_count = parameter_file_dihedrals->GetTerms().size();
                            improper_counter += terms_count;
                            break;
                        }
                    }
                }
            }
        }
    }
    //    cout << not_found_counter/2 << " dihedrals not found in parameter file" << endl;
    return counter/2 + improper_counter;
}

int Assembly::CountNumberOfDihedralsExcludingHydrogen(string parameter_file_path)
{
    ParameterFile* parameter_file = new ParameterFile(parameter_file_path);
    AtomVector atoms = GetAllAtomsOfAssembly();
    int counter = 0;
    int improper_counter = 0;
    //    int not_found_counter = 0;
    for(AtomVector::iterator it = atoms.begin(); it != atoms.end(); it++)
    {
        Atom* atom = (*it);
        stringstream ss;
        ss << atom->GetId();
        string atom_name = atom->GetName();
        AtomNode* atom_node = atom->GetNode();
        if(atom_node != NULL)
        {
            AtomVector neighbors = atom_node->GetNodeNeighbors();
            for(AtomVector::iterator it1 = neighbors.begin(); it1 != neighbors.end(); it1++)
            {
                Atom* neighbor = (*it1);
                stringstream ss1;
                ss1 << neighbor->GetId();
                string neighbor_name = neighbor->GetName();
                AtomNode* neighbor_node = neighbor->GetNode();
                AtomVector neighbors_of_neighbor = neighbor_node->GetNodeNeighbors();
                for(AtomVector::iterator it2 = neighbors_of_neighbor.begin(); it2 != neighbors_of_neighbor.end(); it2++)
                {
                    Atom* neighbor_of_neighbor = (*it2);
                    stringstream ss2;
                    ss2 << neighbor_of_neighbor->GetId();
                    if(ss.str().compare(ss2.str()) != 0)
                    {
                        string neighbor_of_neighbor_name = neighbor_of_neighbor->GetName();
                        AtomNode* neighbor_of_neighbor_node = neighbor_of_neighbor->GetNode();
                        AtomVector neighbors_of_neighbor_of_neighbor = neighbor_of_neighbor_node->GetNodeNeighbors();
                        for(AtomVector::iterator it3 = neighbors_of_neighbor_of_neighbor.begin(); it3 != neighbors_of_neighbor_of_neighbor.end(); it3++)
                        {
                            Atom* neighbor_of_neighbor_of_neighbor = (*it3);
                            stringstream ss3;
                            ss3 << neighbor_of_neighbor_of_neighbor->GetId();
                            if(ss1.str().compare(ss3.str()) != 0)
                            {


                                string neighbor_of_neighbor_of_neighbor_name = neighbor_of_neighbor_of_neighbor->GetName();
                                if((atom_name.substr(0,1).compare("H") == 0 || (atom_name.substr(1,1).compare("H") == 0 && isdigit(ConvertString<char>(atom_name.substr(0,1))))) ||
                                        (neighbor_name.substr(0,1).compare("H") == 0 || (neighbor_name.substr(1,1).compare("H") == 0 && isdigit(ConvertString<char>(neighbor_name.substr(0,1))))) ||
                                        (neighbor_of_neighbor_name.substr(0,1).compare("H") == 0 || (neighbor_of_neighbor_name.substr(1,1).compare("H") == 0 && isdigit(ConvertString<char>(neighbor_of_neighbor_name.substr(0,1))))) ||
                                        (neighbor_of_neighbor_of_neighbor_name.substr(0,1).compare("H") == 0 || (neighbor_of_neighbor_of_neighbor_name.substr(1,1).compare("H") == 0 && isdigit(ConvertString<char>(neighbor_of_neighbor_of_neighbor_name.substr(0,1))))))
                                {}
                                else
                                {

                                    vector<vector<string> > all_atom_type_permutations = CreateAllAtomTypePermutationsforDihedralType(atom->GetAtomType(), neighbor->GetAtomType(),
                                                                                                                                      neighbor_of_neighbor->GetAtomType(), neighbor_of_neighbor_of_neighbor->GetAtomType());
                                    ParameterFile::DihedralMap dihedrals = parameter_file->GetDihedrals();
                                    for(vector<vector<string> >::iterator it4 = all_atom_type_permutations.begin(); it4 != all_atom_type_permutations.end(); it4++)
                                    {
                                        vector<string> atom_types = (*it4);
                                        if(dihedrals[atom_types] != NULL)
                                        {
                                            ParameterFileDihedral* parameter_file_dihedrals = dihedrals[atom_types];
                                            int terms_count = parameter_file_dihedrals->GetTerms().size();
                                            counter += terms_count;
                                            break;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            ///Improper Dihedrals
            if(neighbors.size() == 3)
            {
                Atom* neighbor1 = neighbors.at(0);
                Atom* neighbor2 = neighbors.at(1);
                Atom* neighbor3 = neighbors.at(2);
                string neighbor1_name = neighbor1->GetName();
                string neighbor2_name = neighbor2->GetName();
                string neighbor3_name = neighbor3->GetName();
                vector<vector<string> > all_improper_dihedrals_atom_type_permutations = CreateAllAtomTypePermutationsforImproperDihedralType(neighbor1->GetAtomType(), neighbor2->GetAtomType(),
                                                                                                                                             neighbor3->GetAtomType(), atom->GetAtomType());

                ParameterFile::DihedralMap dihedrals = parameter_file->GetDihedrals();
                for(vector<vector<string> >::iterator it1 = all_improper_dihedrals_atom_type_permutations.begin(); it1 != all_improper_dihedrals_atom_type_permutations.end(); it1++)
                {
                    vector<string> improper_dihedral_permutation = (*it1);
                    if((atom_name.substr(0,1).compare("H") == 0 || (atom_name.substr(1,1).compare("H") == 0 && isdigit(ConvertString<char>(atom_name.substr(0,1))))) ||
                            (neighbor1_name.substr(0,1).compare("H") == 0 || (neighbor1_name.substr(1,1).compare("H") == 0 && isdigit(ConvertString<char>(neighbor1_name.substr(0,1))))) ||
                            (neighbor2_name.substr(0,1).compare("H") == 0 || (neighbor2_name.substr(1,1).compare("H") == 0 && isdigit(ConvertString<char>(neighbor2_name.substr(0,1))))) ||
                            (neighbor3_name.substr(0,1).compare("H") == 0 || (neighbor3_name.substr(1,1).compare("H") == 0 && isdigit(ConvertString<char>(neighbor3_name.substr(0,1))))))
                    {}
                    else
                    {
                        if(dihedrals[improper_dihedral_permutation] != NULL)
                        {
                            ParameterFileDihedral* parameter_file_dihedrals = dihedrals[improper_dihedral_permutation];
                            int terms_count = parameter_file_dihedrals->GetTerms().size();
                            improper_counter += terms_count;
                            break;
                        }
                    }
                }
            }
        }
    }
    //    cout << not_found_counter/2 << " dihedrals not found in parameter file" << endl;
    return counter/2 + improper_counter;
}

int Assembly::CountNumberOfDihedrals(string parameter_file_path)
{
    ParameterFile* parameter_file = new ParameterFile(parameter_file_path);
    AtomVector atoms = GetAllAtomsOfAssembly();
    int counter = 0;
    int improper_counter = 0;
    //    int not_found_counter = 0;
    for(AtomVector::iterator it = atoms.begin(); it != atoms.end(); it++)
    {
        Atom* atom = (*it);
        stringstream ss;
        ss << atom->GetId();
        AtomNode* atom_node = atom->GetNode();
        if(atom_node != NULL)
        {
            AtomVector neighbors = atom_node->GetNodeNeighbors();
            for(AtomVector::iterator it1 = neighbors.begin(); it1 != neighbors.end(); it1++)
            {
                Atom* neighbor = (*it1);
                stringstream ss1;
                ss1 << neighbor->GetId();
                AtomNode* neighbor_node = neighbor->GetNode();
                AtomVector neighbors_of_neighbor = neighbor_node->GetNodeNeighbors();
                for(AtomVector::iterator it2 = neighbors_of_neighbor.begin(); it2 != neighbors_of_neighbor.end(); it2++)
                {
                    Atom* neighbor_of_neighbor = (*it2);
                    stringstream ss2;
                    ss2 << neighbor_of_neighbor->GetId();
                    if(ss.str().compare(ss2.str()) != 0)
                    {
                        AtomNode* neighbor_of_neighbor_node = neighbor_of_neighbor->GetNode();
                        AtomVector neighbors_of_neighbor_of_neighbor = neighbor_of_neighbor_node->GetNodeNeighbors();
                        for(AtomVector::iterator it3 = neighbors_of_neighbor_of_neighbor.begin(); it3 != neighbors_of_neighbor_of_neighbor.end(); it3++)
                        {
                            Atom* neighbor_of_neighbor_of_neighbor = (*it3);
                            stringstream ss3;
                            ss3 << neighbor_of_neighbor_of_neighbor->GetId();
                            if(ss1.str().compare(ss3.str()) != 0)
                            {
                                vector<vector<string> > all_atom_type_permutations = CreateAllAtomTypePermutationsforDihedralType(atom->GetAtomType(), neighbor->GetAtomType(),
                                                                                                                                  neighbor_of_neighbor->GetAtomType(), neighbor_of_neighbor_of_neighbor->GetAtomType());
                                ParameterFile::DihedralMap dihedrals = parameter_file->GetDihedrals();
                                for(vector<vector<string> >::iterator it4 = all_atom_type_permutations.begin(); it4 != all_atom_type_permutations.end(); it4++)
                                {
                                    vector<string> atom_types = (*it4);
                                    //                                cout << atom_types.at(0) << atom_types.at(1) << atom_types.at(2) << atom_types.at(3) << endl;
                                    if(dihedrals[atom_types] != NULL)
                                    {
                                        ParameterFileDihedral* parameter_file_dihedrals = dihedrals[atom_types];
                                        int terms_count = parameter_file_dihedrals->GetTerms().size();
                                        counter += terms_count;
                                        break;
                                    }

                                }
                            }
                        }
                    }
                }
            }
            ///Improper Dihedrals
            if(neighbors.size() == 3)
            {

                Atom* neighbor1 = neighbors.at(0);
                Atom* neighbor2 = neighbors.at(1);
                Atom* neighbor3 = neighbors.at(2);

                vector<vector<string> > all_improper_dihedrals_atom_type_permutations = CreateAllAtomTypePermutationsforImproperDihedralType(neighbor1->GetAtomType(), neighbor2->GetAtomType(),
                                                                                                                                             neighbor3->GetAtomType(), atom->GetAtomType());
                ParameterFile::DihedralMap dihedrals = parameter_file->GetDihedrals();
                for(vector<vector<string> >::iterator it1 = all_improper_dihedrals_atom_type_permutations.begin(); it1 != all_improper_dihedrals_atom_type_permutations.end(); it1++)
                {
                    vector<string> improper_dihedral_permutation = (*it1);
                    if(dihedrals[improper_dihedral_permutation] != NULL)
                    {
                        ParameterFileDihedral* parameter_file_dihedrals = dihedrals[improper_dihedral_permutation];
                        int terms_count = parameter_file_dihedrals->GetTerms().size();
                        improper_counter += terms_count;
                        break;
                    }
                }
            }
        }
    }
    //    cout << not_found_counter/2 << " dihedrals not found in parameter file" << endl;
    return counter/2 + improper_counter;
}

int Assembly::CountNumberOfDihedralTypes(string parameter_file_path)
{
    vector<string> type_list = vector<string>();
    ParameterFile* parameter_file = new ParameterFile(parameter_file_path);
    AtomVector atoms = GetAllAtomsOfAssembly();
    int counter = 0;
    //    int not_found_counter = 0;
    for(AtomVector::iterator it = atoms.begin(); it != atoms.end(); it++)
    {
        Atom* atom = (*it);
        stringstream ss;
        ss << atom->GetId();
        AtomNode* atom_node = atom->GetNode();
        if(atom_node != NULL)
        {
            AtomVector neighbors = atom_node->GetNodeNeighbors();
            for(AtomVector::iterator it1 = neighbors.begin(); it1 != neighbors.end(); it1++)
            {
                Atom* neighbor = (*it1);
                stringstream ss1;
                ss1 << neighbor->GetId();
                AtomNode* neighbor_node = neighbor->GetNode();
                AtomVector neighbors_of_neighbor = neighbor_node->GetNodeNeighbors();
                for(AtomVector::iterator it2 = neighbors_of_neighbor.begin(); it2 != neighbors_of_neighbor.end(); it2++)
                {
                    Atom* neighbor_of_neighbor = (*it2);
                    stringstream ss2;
                    ss2 << neighbor_of_neighbor->GetId();
                    if(ss.str().compare(ss2.str()) != 0)
                    {
                        AtomNode* neighbor_of_neighbor_node = neighbor_of_neighbor->GetNode();
                        AtomVector neighbors_of_neighbor_of_neighbor = neighbor_of_neighbor_node->GetNodeNeighbors();
                        for(AtomVector::iterator it3 = neighbors_of_neighbor_of_neighbor.begin(); it3 != neighbors_of_neighbor_of_neighbor.end(); it3++)
                        {
                            Atom* neighbor_of_neighbor_of_neighbor = (*it3);
                            stringstream ss3;
                            ss3 << neighbor_of_neighbor_of_neighbor->GetId();
                            if(ss1.str().compare(ss3.str()) != 0)
                            {

                                vector<vector<string> > all_atom_type_permutations = CreateAllAtomTypePermutationsforDihedralType(atom->GetAtomType(), neighbor->GetAtomType(),
                                                                                                                                  neighbor_of_neighbor->GetAtomType(), neighbor_of_neighbor_of_neighbor->GetAtomType());
                                ParameterFile::DihedralMap dihedrals = parameter_file->GetDihedrals();
                                for(vector<vector<string> >::iterator it4 = all_atom_type_permutations.begin(); it4 != all_atom_type_permutations.end(); it4++)
                                {
                                    vector<string> atom_types = (*it4);


                                    //                            vector<string> atom_types = vector<string>();
                                    //                            atom_types.push_back(atom->GetAtomType());
                                    //                            atom_types.push_back(neighbor->GetAtomType());
                                    //                            atom_types.push_back(neighbor_of_neighbor->GetAtomType());
                                    //                            atom_types.push_back(neighbor_of_neighbor_of_neighbor->GetAtomType());
                                    if(dihedrals[atom_types] != NULL)
                                    {
                                        stringstream ss4;
                                        ss4 << atom_types.at(0) << "_" << atom_types.at(1) << "_" << atom_types.at(2) << "_" << atom_types.at(3);
                                        stringstream ss5;
                                        ss5 << atom_types.at(3) << "_" << atom_types.at(2) << "_" << atom_types.at(1) << "_" << atom_types.at(0);
                                        if(find(type_list.begin(), type_list.end(), ss4.str()) == type_list.end() &&
                                                find(type_list.begin(), type_list.end(), ss5.str()) == type_list.end())
                                        {
                                            type_list.push_back(ss4.str());
                                            ParameterFileDihedral* parameter_file_dihedrals = dihedrals[atom_types];
                                            int terms_count = parameter_file_dihedrals->GetTerms().size();
                                            counter += terms_count;
                                            break;
                                        }
                                    }
                                    //                                else
                                    //                                {
                                    //                                    atom_types[0] = neighbor_of_neighbor_of_neighbor->GetAtomType();
                                    //                                    atom_types[1] = neighbor_of_neighbor->GetAtomType();
                                    //                                    atom_types[2] = neighbor->GetAtomType();
                                    //                                    atom_types[3] = atom->GetAtomType();
                                    //                                    if(dihedrals[atom_types] != NULL)
                                    //                                    {
                                    //                                        stringstream ss6;
                                    //                                        ss6 << atom->GetAtomType() << "_" << neighbor->GetAtomType() << "_" << neighbor_of_neighbor->GetAtomType() << "_" << neighbor_of_neighbor_of_neighbor->GetAtomType();
                                    //                                        stringstream ss7;
                                    //                                        ss7 << neighbor_of_neighbor_of_neighbor->GetAtomType() << "_" << neighbor_of_neighbor->GetAtomType() << "_" << neighbor->GetAtomType() << "_" << atom->GetAtomType();
                                    //                                        if(find(type_list.begin(), type_list.end(), ss6.str()) == type_list.end() &&
                                    //                                                find(type_list.begin(), type_list.end(), ss7.str()) == type_list.end() )
                                    //                                        {
                                    //                                            type_list.push_back(ss7.str());
                                    //                                            ParameterFileDihedral* parameter_file_dihedrals = dihedrals[atom_types];
                                    //                                            int terms_count = parameter_file_dihedrals->GetTerms().size();
                                    //                                            counter += terms_count;
                                    //                                        }
                                    //                                    }
                                    //                                    else
                                    //                                    {
                                    //                                        not_found_counter++;
                                    //                                    }
                                    //                                }
                                }
                            }
                        }
                    }
                }
            }
            ///Improper Dihedrals
            if(neighbors.size() == 3)
            {
                Atom* neighbor1 = neighbors.at(0);
                Atom* neighbor2 = neighbors.at(1);
                Atom* neighbor3 = neighbors.at(2);
                vector<vector<string> > all_improper_dihedrals_atom_type_permutations = CreateAllAtomTypePermutationsforImproperDihedralType(neighbor1->GetAtomType(), neighbor2->GetAtomType(),
                                                                                                                                             neighbor3->GetAtomType(), atom->GetAtomType());

                for(vector<vector<string> >::iterator it1 = all_improper_dihedrals_atom_type_permutations.begin(); it1 != all_improper_dihedrals_atom_type_permutations.end(); it1++)
                {
                    vector<string> improper_dihedral_permutation = (*it1);
                    ParameterFile::DihedralMap dihedrals = parameter_file->GetDihedrals();
                    if(dihedrals[improper_dihedral_permutation] != NULL)
                    {
                        stringstream ss;
                        ss << improper_dihedral_permutation.at(0) << "_" << improper_dihedral_permutation.at(1) << "_" << improper_dihedral_permutation.at(2) << "_" << improper_dihedral_permutation.at(3);
                        stringstream ss1;
                        ss1 << improper_dihedral_permutation.at(3) << "_" << improper_dihedral_permutation.at(2) << "_" << improper_dihedral_permutation.at(1) << "_" << improper_dihedral_permutation.at(0);
                        if(find(type_list.begin(), type_list.end(), ss.str()) == type_list.end() &&
                                find(type_list.begin(), type_list.end(), ss1.str()) == type_list.end())
                        {
                            type_list.push_back(ss.str());
                            ParameterFileDihedral* parameter_file_dihedrals = dihedrals[improper_dihedral_permutation];
                            int terms_count = parameter_file_dihedrals->GetTerms().size();
                            counter += terms_count;
                            break;
                        }
                    }
                }
            }
        }
    }
    //    cout << not_found_counter << " dihedrals not found in parameter file" << endl;
    //    cout << type_list.size() << endl;
    return counter;
}

vector<vector<string> > Assembly::CreateAllAtomTypePermutationsforDihedralType(string atom_type1, string atom_type2, string atom_type3, string atom_type4)
{
    vector<vector<string> > all_permutations = vector<vector<string> >();
    vector<string> normal_order = vector<string>();
    vector<string> reverse_order = vector<string>();
    normal_order.push_back(atom_type1);
    normal_order.push_back(atom_type2);
    normal_order.push_back(atom_type3);
    normal_order.push_back(atom_type4);
    reverse_order.push_back(atom_type4);
    reverse_order.push_back(atom_type3);
    reverse_order.push_back(atom_type2);
    reverse_order.push_back(atom_type1);
    all_permutations.push_back(normal_order);
    all_permutations.push_back(reverse_order);

    ///Permutations w/ one X
    for(int i = 0; i < 4; i++)
    {
        normal_order.clear();
        normal_order.push_back(atom_type1);
        normal_order.push_back(atom_type2);
        normal_order.push_back(atom_type3);
        normal_order.push_back(atom_type4);
        normal_order.at(i) = "X";
        reverse_order.clear();
        reverse_order.push_back(atom_type4);
        reverse_order.push_back(atom_type3);
        reverse_order.push_back(atom_type2);
        reverse_order.push_back(atom_type1);
        reverse_order.at(i) = "X";
        all_permutations.push_back(normal_order);
        all_permutations.push_back(reverse_order);
    }
    ///Permutations w/ two X
    for(int i = 0; i < 3; i++)
        for(int j = i+1; j < 4; j++)
        {
            normal_order.clear();
            normal_order.push_back(atom_type1);
            normal_order.push_back(atom_type2);
            normal_order.push_back(atom_type3);
            normal_order.push_back(atom_type4);
            normal_order.at(i) = "X";
            normal_order.at(j) = "X";
            reverse_order.clear();
            reverse_order.push_back(atom_type4);
            reverse_order.push_back(atom_type3);
            reverse_order.push_back(atom_type2);
            reverse_order.push_back(atom_type1);
            reverse_order.at(i) = "X";
            reverse_order.at(j) = "X";
            all_permutations.push_back(normal_order);
            all_permutations.push_back(reverse_order);
        }
    ///Permutations w/ three X
    for(int i = 0; i < 2; i++)
        for(int j = i+1; j < 3; j++)
            for(int k = j+1; k < 4; k++)
            {
                normal_order.clear();
                normal_order.push_back(atom_type1);
                normal_order.push_back(atom_type2);
                normal_order.push_back(atom_type3);
                normal_order.push_back(atom_type4);
                normal_order.at(i) = "X";
                normal_order.at(j) = "X";
                normal_order.at(k) = "X";
                reverse_order.clear();
                reverse_order.push_back(atom_type4);
                reverse_order.push_back(atom_type3);
                reverse_order.push_back(atom_type2);
                reverse_order.push_back(atom_type1);
                reverse_order.at(i) = "X";
                reverse_order.at(j) = "X";
                reverse_order.at(k) = "X";
                all_permutations.push_back(normal_order);
                all_permutations.push_back(reverse_order);
            }
    return all_permutations;
}

vector<vector<string> > Assembly::CreateAllAtomTypePermutationsforImproperDihedralType(string neighbor1_type, string neighbor2_type, string neighbor3_type, string atom_type)
{
    vector<vector<string> > all_permutations = vector<vector<string> >();
    vector<string> order1 = vector<string>();
    vector<string> order2 = vector<string>();
    vector<string> order3 = vector<string>();
    vector<string> reverse_order1 = vector<string>();
    vector<string> reverse_order2 = vector<string>();
    vector<string> reverse_order3 = vector<string>();

    order1.push_back(neighbor1_type);
    order1.push_back(neighbor2_type);
    order1.push_back(atom_type);
    order1.push_back(neighbor3_type);
    reverse_order1.push_back(neighbor3_type);
    reverse_order1.push_back(atom_type);
    reverse_order1.push_back(neighbor2_type);
    reverse_order1.push_back(neighbor1_type);

    order2.push_back(neighbor1_type);
    order2.push_back(atom_type);
    order2.push_back(neighbor3_type);
    order2.push_back(neighbor2_type);
    reverse_order2.push_back(neighbor2_type);
    reverse_order2.push_back(neighbor3_type);
    reverse_order2.push_back(atom_type);
    reverse_order2.push_back(neighbor1_type);

    order3.push_back(neighbor1_type);
    order3.push_back(neighbor3_type);
    order3.push_back(atom_type);
    order3.push_back(neighbor2_type);
    reverse_order3.push_back(neighbor2_type);
    reverse_order3.push_back(atom_type);
    reverse_order3.push_back(neighbor3_type);
    reverse_order3.push_back(neighbor1_type);

    all_permutations.push_back(order1);
    all_permutations.push_back(order2);
    all_permutations.push_back(order3);
    all_permutations.push_back(reverse_order1);
    all_permutations.push_back(reverse_order2);
    all_permutations.push_back(reverse_order3);

    ///Permutations w/ one X
    for(int i = 0; i < 4; i++)
    {
        order1.clear();
        order2.clear();
        order3.clear();
        reverse_order1.clear();
        reverse_order2.clear();
        reverse_order3.clear();

        order1.push_back(neighbor1_type);
        order1.push_back(neighbor2_type);
        order1.push_back(atom_type);
        order1.push_back(neighbor3_type);
        order1.at(i) = "X";
        reverse_order1.push_back(neighbor3_type);
        reverse_order1.push_back(atom_type);
        reverse_order1.push_back(neighbor2_type);
        reverse_order1.push_back(neighbor1_type);
        reverse_order1.at(i) = "X";
        order2.push_back(neighbor1_type);
        order2.push_back(atom_type);
        order2.push_back(neighbor3_type);
        order2.push_back(neighbor2_type);
        order2.at(i) = "X";
        reverse_order2.push_back(neighbor2_type);
        reverse_order2.push_back(neighbor3_type);
        reverse_order2.push_back(atom_type);
        reverse_order2.push_back(neighbor1_type);
        reverse_order2.at(i) = "X";
        order3.push_back(neighbor1_type);
        order3.push_back(neighbor3_type);
        order3.push_back(atom_type);
        order3.push_back(neighbor2_type);
        order3.at(i) = "X";
        reverse_order3.push_back(neighbor2_type);
        reverse_order3.push_back(atom_type);
        reverse_order3.push_back(neighbor3_type);
        reverse_order3.push_back(neighbor1_type);
        reverse_order3.at(i) = "X";

        all_permutations.push_back(order1);
        all_permutations.push_back(reverse_order1);
        all_permutations.push_back(order2);
        all_permutations.push_back(reverse_order2);
        all_permutations.push_back(order3);
        all_permutations.push_back(reverse_order3);
    }
    ///Permutations w/ two X
    for(int i = 0; i < 3; i++)
        for(int j = i+1; j < 4; j++)
        {
            order1.clear();
            order2.clear();
            order3.clear();
            reverse_order1.clear();
            reverse_order2.clear();
            reverse_order3.clear();

            order1.push_back(neighbor1_type);
            order1.push_back(neighbor2_type);
            order1.push_back(atom_type);
            order1.push_back(neighbor3_type);
            order1.at(i) = "X";
            order1.at(j) = "X";
            reverse_order1.push_back(neighbor3_type);
            reverse_order1.push_back(atom_type);
            reverse_order1.push_back(neighbor2_type);
            reverse_order1.push_back(neighbor1_type);
            reverse_order1.at(i) = "X";
            order2.push_back(neighbor1_type);
            order2.push_back(atom_type);
            order2.push_back(neighbor3_type);
            order2.push_back(neighbor2_type);
            order2.at(i) = "X";
            order2.at(j) = "X";
            reverse_order2.push_back(neighbor2_type);
            reverse_order2.push_back(neighbor3_type);
            reverse_order2.push_back(atom_type);
            reverse_order2.push_back(neighbor1_type);
            reverse_order2.at(i) = "X";
            reverse_order2.at(j) = "X";
            order3.push_back(neighbor1_type);
            order3.push_back(neighbor3_type);
            order3.push_back(atom_type);
            order3.push_back(neighbor2_type);
            order3.at(i) = "X";
            order3.at(j) = "X";
            reverse_order3.push_back(neighbor2_type);
            reverse_order3.push_back(atom_type);
            reverse_order3.push_back(neighbor3_type);
            reverse_order3.push_back(neighbor1_type);
            reverse_order3.at(i) = "X";
            reverse_order3.at(j) = "X";

            all_permutations.push_back(order1);
            all_permutations.push_back(reverse_order1);
            all_permutations.push_back(order2);
            all_permutations.push_back(reverse_order2);
            all_permutations.push_back(order3);
            all_permutations.push_back(reverse_order3);;
        }
    ///Permutations w/ three X
    for(int i = 0; i < 2; i++)
        for(int j = i+1; j < 3; j++)
            for(int k = j+1; k < 4; k++)
            {
                order1.clear();
                order2.clear();
                order3.clear();
                reverse_order1.clear();
                reverse_order2.clear();
                reverse_order3.clear();

                order1.push_back(neighbor1_type);
                order1.push_back(neighbor2_type);
                order1.push_back(atom_type);
                order1.push_back(neighbor3_type);
                order1.at(i) = "X";
                order1.at(j) = "X";
                order1.at(k) = "X";
                reverse_order1.push_back(neighbor3_type);
                reverse_order1.push_back(atom_type);
                reverse_order1.push_back(neighbor2_type);
                reverse_order1.push_back(neighbor1_type);
                reverse_order1.at(i) = "X";
                order2.push_back(neighbor1_type);
                order2.push_back(atom_type);
                order2.push_back(neighbor3_type);
                order2.push_back(neighbor2_type);
                order2.at(i) = "X";
                order2.at(j) = "X";
                order2.at(k) = "X";
                reverse_order2.push_back(neighbor2_type);
                reverse_order2.push_back(neighbor3_type);
                reverse_order2.push_back(atom_type);
                reverse_order2.push_back(neighbor1_type);
                reverse_order2.at(i) = "X";
                reverse_order2.at(j) = "X";
                reverse_order2.at(k) = "X";
                order3.push_back(neighbor1_type);
                order3.push_back(neighbor3_type);
                order3.push_back(atom_type);
                order3.push_back(neighbor2_type);
                order3.at(i) = "X";
                order3.at(j) = "X";
                order3.at(k) = "X";
                reverse_order3.push_back(neighbor2_type);
                reverse_order3.push_back(atom_type);
                reverse_order3.push_back(neighbor3_type);
                reverse_order3.push_back(neighbor1_type);
                reverse_order3.at(i) = "X";
                reverse_order3.at(j) = "X";
                reverse_order3.at(k) = "X";

                all_permutations.push_back(order1);
                all_permutations.push_back(reverse_order1);
                all_permutations.push_back(order2);
                all_permutations.push_back(reverse_order2);
                all_permutations.push_back(order3);
                all_permutations.push_back(reverse_order3);;
            }
    return all_permutations;
}

Assembly::AtomVector Assembly::GetAllAtomsOfAssemblyWithAtLeastThreeNeighbors()
{
    AtomVector all_atoms = GetAllAtomsOfAssembly();
    AtomVector atoms_with_at_least_three_neighbors = AtomVector();
    for(AtomVector::iterator it = all_atoms.begin(); it != all_atoms.end(); it++)
    {
        Atom* atom = (*it);
        AtomNode* node = atom->GetNode();
        AtomVector neighbors = node->GetNodeNeighbors();
        if(neighbors.size() > 2)
            atoms_with_at_least_three_neighbors.push_back(atom);
    }
    return atoms_with_at_least_three_neighbors;
}

int Assembly::CountNumberOfExcludedAtoms()
{
    AtomVector atoms = GetAllAtomsOfAssembly();
    vector<string> excluded_atom_list = vector<string>();
    map<string, vector<string> > excluded_atom_list_map = map<string, vector<string> >();
    for(AtomVector::iterator it = atoms.begin(); it != atoms.end(); it++)
    {
        Atom* atom = (*it);
        stringstream ss;
        ss << atom->GetId();
        AtomNode* node = atom->GetNode();
        excluded_atom_list_map[atom->GetId()] = vector<string>();
        if(node != NULL)
        {
            AtomVector neighbors = node->GetNodeNeighbors();
            for(AtomVector::iterator it1 = neighbors.begin(); it1 != neighbors.end(); it1++)
            {
                Atom* neighbor = (*it1);
                stringstream ss1;
                ss1 << neighbor->GetId();
                stringstream first_order_interaction;
                stringstream reverse_first_order_interaction;
                first_order_interaction << ss.str() << "-" << ss1.str();
                reverse_first_order_interaction << ss1.str() << "-" << ss.str();
                if(find(excluded_atom_list.begin(), excluded_atom_list.end(), first_order_interaction.str()) == excluded_atom_list.end() &&
                        find(excluded_atom_list.begin(), excluded_atom_list.end(), reverse_first_order_interaction.str()) == excluded_atom_list.end())
                {
                    excluded_atom_list.push_back(first_order_interaction.str());
                    excluded_atom_list_map[atom->GetId()].push_back(neighbor->GetId());
                }
                AtomNode* neighbor_node = neighbor->GetNode();
                AtomVector neighbor_of_neighbors = neighbor_node->GetNodeNeighbors();
                for(AtomVector::iterator it2 = neighbor_of_neighbors.begin(); it2 != neighbor_of_neighbors.end(); it2++)
                {
                    Atom* neighbor_of_neighbor = (*it2);
                    stringstream ss2;
                    ss2 << neighbor_of_neighbor->GetId();
                    if(ss.str().compare(ss2.str()) != 0)
                    {
                        stringstream second_order_interaction;
                        stringstream reverse_second_order_interaction;
                        second_order_interaction << ss.str() << "-" << ss2.str();
                        reverse_second_order_interaction << ss2.str() << "-" << ss.str();
                        if(find(excluded_atom_list.begin(), excluded_atom_list.end(), second_order_interaction.str()) == excluded_atom_list.end() &&
                                find(excluded_atom_list.begin(), excluded_atom_list.end(), reverse_second_order_interaction.str()) == excluded_atom_list.end())
                        {
                            excluded_atom_list.push_back(second_order_interaction.str());
                            excluded_atom_list_map[atom->GetId()].push_back(neighbor_of_neighbor->GetId());
                        }
                        AtomNode* neighbor_of_neighbor_node = neighbor_of_neighbor->GetNode();
                        AtomVector neighbor_of_neighbor_of_neighbors = neighbor_of_neighbor_node->GetNodeNeighbors();
                        for(AtomVector::iterator it3 = neighbor_of_neighbor_of_neighbors.begin(); it3 != neighbor_of_neighbor_of_neighbors.end(); it3++)
                        {
                            Atom* neighbor_of_neighbor_of_neighbor = (*it3);
                            stringstream ss3;
                            ss3 << neighbor_of_neighbor_of_neighbor->GetId();
                            if(ss1.str().compare(ss3.str()) != 0)
                            {
                                stringstream third_order_interaction;
                                stringstream reverse_third_order_interaction;
                                third_order_interaction << ss.str() << "-" << ss3.str();
                                reverse_third_order_interaction << ss3.str() << "-" << ss.str();
                                if(find(excluded_atom_list.begin(), excluded_atom_list.end(), third_order_interaction.str()) == excluded_atom_list.end() &&
                                        find(excluded_atom_list.begin(), excluded_atom_list.end(), reverse_third_order_interaction.str()) == excluded_atom_list.end())
                                {
                                    excluded_atom_list.push_back(third_order_interaction.str());
                                    excluded_atom_list_map[atom->GetId()].push_back(neighbor_of_neighbor_of_neighbor->GetId());
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    int number_of_excluded_atoms = 0;
    for(map<string, vector<string> >::iterator it = excluded_atom_list_map.begin(); it != excluded_atom_list_map.end(); it++)
        (*it).second.size() == 0 ? number_of_excluded_atoms++ : number_of_excluded_atoms += (*it).second.size();
    return number_of_excluded_atoms;
}

int Assembly::CountMaxNumberOfAtomsInLargestResidue()
{
    int max = 0;
    for(AssemblyVector::iterator it = assemblies_.begin(); it != assemblies_.end(); it++)
    {
        Assembly* assembly = (*it);
        if(max <= assembly->CountNumberOfAtoms())
            max = assembly->CountNumberOfAtoms();
    }
    for(ResidueVector:: iterator it1 = residues_.begin(); it1 != residues_.end(); it1++)
    {
        Residue* residue = (*it1);
        Residue::AtomVector atoms = residue->GetAtoms();
        if(max <= atoms.size())
            max = atoms.size();
    }
    return max;
}

Assembly::AtomVector Assembly::Select(string pattern)
{
    AtomVector selection = AtomVector();

    HierarchicalContainmentMap hierarchical_map;
    stringstream index;
    index << this->GetSequenceNumber();
    this->GetHierarchicalMapOfAssembly(hierarchical_map, index);

    for(HierarchicalContainmentMap::iterator it = hierarchical_map.begin(); it != hierarchical_map.end(); it++)
        cout << (*it).first << " " << (*it).second.size() << endl;

    SelectPatternMap select_pattern_map = ParsePatternString(pattern);

    //    for(SelectPatternMap::iterator it = select_pattern_map.begin(); it != select_pattern_map.end(); it++)
    //    {
    //        cout << "================================= " << (*it).first << " ================================" << endl;
    //        map<string, vector<string> > value = (*it).second;
    //        for(map<string, vector<string> >::iterator it1 = value.begin(); it1 != value.end(); it1++)
    //        {
    //            cout << "-------------------------------- " << (*it1).first << " --------------------------" << endl;
    //            vector<string> atom_names = (*it1).second;
    //            for(vector<string>::iterator it3 = atom_names.begin(); it3 != atom_names.end(); it3++)
    //            {
    //                cout << (*it3) << endl;
    //            }
    //        }
    //    }

    for(SelectPatternMap::iterator it = select_pattern_map.begin(); it != select_pattern_map.end(); it++)
    {
        string key = (*it).first;
        if(key.find("*") == string::npos)
        {
            ResidueVector residues_of_assembly = hierarchical_map.at(key);
            map<string, vector<string> > value = (*it).second;
            for(map<string, vector<string> >::iterator it1 = value.begin(); it1 != value.end(); it1++)
            {
                string residue_name = (*it1).first;
                int residue_name_search_type = 0;
                if(residue_name.at(0) == '^')
                {
                    residue_name = residue_name.substr(1);
                    residue_name_search_type = -1;
                }
                else if(residue_name.at(residue_name.size() - 1) == '$')
                {
                    residue_name = residue_name.substr(0, residue_name.size() - 2);
                    residue_name_search_type = 1;
                }
                else if(residue_name.at(0) == '#')
                {
                    residue_name = residue_name.substr(1);
                    residue_name_search_type = -2;
                }
                else
                    residue_name_search_type = 0;
                if(residue_name.find("*") == string::npos)
                {
                    for(ResidueVector::iterator it2 = residues_of_assembly.begin(); it2 != residues_of_assembly.end(); it2++)
                    {
                        Residue* residue = *it2;
                        switch(residue_name_search_type)
                        {
                            case 0:  /// Search in residue set by matching the whole name of the residue
                            {
                                if(residue->GetName().compare(residue_name) == 0)
                                {
                                    vector<string> atom_names = (*it1).second;
                                    if(find(atom_names.begin(), atom_names.end(), "*") == atom_names.end())
                                    {
                                        AtomVector atoms = residue->GetAtoms();
                                        for(vector<string>::iterator it3 = atom_names.begin(); it3 != atom_names.end(); it3++)
                                        {
                                            string atom_name = *it3;
                                            int atom_name_search_type = 0;
                                            if(atom_name.at(0) == '^')
                                            {
                                                atom_name = atom_name.substr(1);
                                                atom_name_search_type = -1;
                                            }
                                            else if(atom_name.at(atom_name.size() - 1) == '$')
                                            {
                                                atom_name = atom_name.substr(0, atom_name.size() - 2);
                                                atom_name_search_type = 1;
                                            }
                                            else if(atom_name.at(0) == '#')
                                            {
                                                atom_name = atom_name.substr(1);
                                                atom_name_search_type = -2;
                                            }
                                            else
                                                atom_name_search_type = 0;
                                            for(AtomVector::iterator it4 = atoms.begin(); it4 != atoms.end(); it4++)
                                            {
                                                Atom* atom = *it4;
                                                switch(atom_name_search_type)
                                                {
                                                    case 0:
                                                    {
                                                        if(atom->GetName().compare(atom_name) == 0)
                                                            selection.push_back(atom);
                                                        break;
                                                    }
                                                    case 1:
                                                    {
                                                        if(atom->GetName().find(atom_name) != string::npos &&
                                                                atom->GetName().find(atom_name) == atom->GetName().size() - atom_name.size())
                                                            selection.push_back(atom);
                                                        break;
                                                    }
                                                    case -1:
                                                    {
                                                        if(atom->GetName().find(atom_name) != string::npos &&
                                                                atom->GetName().find(atom_name) == 0)
                                                            selection.push_back(atom);
                                                        break;
                                                    }
                                                    case -2:
                                                    {
                                                        string atom_serial_number = Split(atom->GetId(), "_").at(1);
                                                        int range_selection = 0;
                                                        if(atom_name.find("-") != string::npos)
                                                        {
                                                            range_selection = 1;

                                                        }
                                                        else
                                                            range_selection = 0;
                                                        switch(range_selection)
                                                        {
                                                            case 0:
                                                            {
                                                                if(atom_serial_number.find(atom_name) != string::npos &&
                                                                        atom_serial_number.find(atom_name) == 0)
                                                                    selection.push_back(atom);
                                                                break;
                                                            }
                                                            case 1:
                                                            {
                                                                string start_index = Split(atom_name, "-").at(0);
                                                                string end_index = Split(atom_name, "-").at(1);
                                                                int atom_serial_number_int = ConvertString<int>(atom_serial_number);
                                                                int start_index_int = ConvertString<int>(start_index);
                                                                int end_index_int = ConvertString<int>(end_index);
                                                                if(atom_serial_number_int >= start_index_int && atom_serial_number_int <= end_index_int)
                                                                    selection.push_back(atom);
                                                                break;
                                                            }
                                                        }
                                                        break;
                                                    }
                                                }

                                            }
                                        }
                                    }
                                    else
                                    {
                                        AtomVector atoms = residue->GetAtoms();
                                        for(AtomVector::iterator it4 = atoms.begin(); it4 != atoms.end(); it4++)
                                        {
                                            Atom* atom = *it4;
                                            selection.push_back(atom);
                                        }
                                    }
                                }
                                break;
                            }
                            case 1:   /// Searching the residue set by maching the end of the residue names
                            {
                                if(residue->GetName().find(residue_name) != string::npos &&
                                        residue->GetName().find(residue_name) == residue->GetName().size() - residue_name.size())
                                {
                                    vector<string> atom_names = (*it1).second;
                                    if(find(atom_names.begin(), atom_names.end(), "*") == atom_names.end())
                                    {
                                        AtomVector atoms = residue->GetAtoms();
                                        for(vector<string>::iterator it3 = atom_names.begin(); it3 != atom_names.end(); it3++)
                                        {
                                            string atom_name = *it3;
                                            int atom_name_search_type = 0;
                                            if(atom_name.at(0) == '^')
                                            {
                                                atom_name = atom_name.substr(1);
                                                atom_name_search_type = -1;
                                            }
                                            else if(atom_name.at(atom_name.size() - 1) == '$')
                                            {
                                                atom_name = atom_name.substr(0, atom_name.size() - 2);
                                                atom_name_search_type = 1;
                                            }
                                            else if(atom_name.at(0) == '#')
                                            {
                                                atom_name = atom_name.substr(1);
                                                atom_name_search_type = -2;
                                            }
                                            else
                                                atom_name_search_type = 0;
                                            for(AtomVector::iterator it4 = atoms.begin(); it4 != atoms.end(); it4++)
                                            {
                                                Atom* atom = *it4;
                                                switch(atom_name_search_type)
                                                {
                                                    case 0:
                                                    {
                                                        if(atom->GetName().compare(atom_name) == 0)
                                                            selection.push_back(atom);
                                                        break;
                                                    }
                                                    case 1:
                                                    {
                                                        if(atom->GetName().find(atom_name) != string::npos &&
                                                                atom->GetName().find(atom_name) == atom->GetName().size() - atom_name.size())
                                                            selection.push_back(atom);
                                                        break;
                                                    }
                                                    case -1:
                                                    {
                                                        if(atom->GetName().find(atom_name) != string::npos &&
                                                                atom->GetName().find(atom_name) == 0)
                                                            selection.push_back(atom);
                                                        break;
                                                    }
                                                    case -2:
                                                    {
                                                        string atom_serial_number = Split(atom->GetId(), "_").at(1);
                                                        int range_selection = 0;
                                                        if(atom_name.find("-") != string::npos)
                                                        {
                                                            range_selection = 1;

                                                        }
                                                        else
                                                            range_selection = 0;
                                                        switch(range_selection)
                                                        {
                                                            case 0:
                                                            {
                                                                if(atom_serial_number.find(atom_name) != string::npos &&
                                                                        atom_serial_number.find(atom_name) == 0)
                                                                    selection.push_back(atom);
                                                                break;
                                                            }
                                                            case 1:
                                                            {
                                                                string start_index = Split(atom_name, "-").at(0);
                                                                string end_index = Split(atom_name, "-").at(1);
                                                                int atom_serial_number_int = ConvertString<int>(atom_serial_number);
                                                                int start_index_int = ConvertString<int>(start_index);
                                                                int end_index_int = ConvertString<int>(end_index);
                                                                if(atom_serial_number_int >= start_index_int && atom_serial_number_int <= end_index_int)
                                                                    selection.push_back(atom);
                                                                break;
                                                            }
                                                        }
                                                        break;
                                                    }
                                                }

                                            }
                                        }
                                    }
                                    else
                                    {
                                        AtomVector atoms = residue->GetAtoms();
                                        for(AtomVector::iterator it4 = atoms.begin(); it4 != atoms.end(); it4++)
                                        {
                                            Atom* atom = *it4;
                                            selection.push_back(atom);
                                        }
                                    }
                                }
                                break;
                            }
                            case -1:  /// Searching the residue set by matching the begining of the residue names
                            {
                                if(residue->GetName().find(residue_name) != string::npos &&
                                        residue->GetName().find(residue_name) == 0)
                                {
                                    vector<string> atom_names = (*it1).second;
                                    if(find(atom_names.begin(), atom_names.end(), "*") == atom_names.end())
                                    {
                                        AtomVector atoms = residue->GetAtoms();
                                        for(vector<string>::iterator it3 = atom_names.begin(); it3 != atom_names.end(); it3++)
                                        {
                                            string atom_name = *it3;
                                            int atom_name_search_type = 0;
                                            if(atom_name.at(0) == '^')
                                            {
                                                atom_name = atom_name.substr(1);
                                                atom_name_search_type = -1;
                                            }
                                            else if(atom_name.at(atom_name.size() - 1) == '$')
                                            {
                                                atom_name = atom_name.substr(0, atom_name.size() - 2);
                                                atom_name_search_type = 1;
                                            }
                                            else if(atom_name.at(0) == '#')
                                            {
                                                atom_name = atom_name.substr(1);
                                                atom_name_search_type = -2;
                                            }
                                            else
                                                atom_name_search_type = 0;
                                            for(AtomVector::iterator it4 = atoms.begin(); it4 != atoms.end(); it4++)
                                            {
                                                Atom* atom = *it4;
                                                switch(atom_name_search_type)
                                                {
                                                    case 0:
                                                    {
                                                        if(atom->GetName().compare(atom_name) == 0)
                                                            selection.push_back(atom);
                                                        break;
                                                    }
                                                    case 1:
                                                    {
                                                        if(atom->GetName().find(atom_name) != string::npos &&
                                                                atom->GetName().find(atom_name) == atom->GetName().size() - atom_name.size())
                                                            selection.push_back(atom);
                                                        break;
                                                    }
                                                    case -1:
                                                    {
                                                        if(atom->GetName().find(atom_name) != string::npos &&
                                                                atom->GetName().find(atom_name) == 0)
                                                            selection.push_back(atom);
                                                        break;
                                                    }
                                                    case -2:
                                                    {
                                                        string atom_serial_number = Split(atom->GetId(), "_").at(1);
                                                        int range_selection = 0;
                                                        if(atom_name.find("-") != string::npos)
                                                        {
                                                            range_selection = 1;

                                                        }
                                                        else
                                                            range_selection = 0;
                                                        switch(range_selection)
                                                        {
                                                            case 0:
                                                            {
                                                                if(atom_serial_number.find(atom_name) != string::npos &&
                                                                        atom_serial_number.find(atom_name) == 0)
                                                                    selection.push_back(atom);
                                                                break;
                                                            }
                                                            case 1:
                                                            {
                                                                string start_index = Split(atom_name, "-").at(0);
                                                                string end_index = Split(atom_name, "-").at(1);
                                                                int atom_serial_number_int = ConvertString<int>(atom_serial_number);
                                                                int start_index_int = ConvertString<int>(start_index);
                                                                int end_index_int = ConvertString<int>(end_index);
                                                                if(atom_serial_number_int >= start_index_int && atom_serial_number_int <= end_index_int)
                                                                    selection.push_back(atom);
                                                                break;
                                                            }
                                                        }
                                                        break;
                                                    }
                                                }

                                            }
                                        }
                                    }
                                    else
                                    {
                                        AtomVector atoms = residue->GetAtoms();
                                        for(AtomVector::iterator it4 = atoms.begin(); it4 != atoms.end(); it4++)
                                        {
                                            Atom* atom = *it4;
                                            selection.push_back(atom);
                                        }
                                    }
                                }
                                break;
                            }
                            case -2:  /// Searching the residue set by matching the residue sequence number
                            {
                                string residue_sequence_number = Split(residue->GetId(), "_").at(2);
                                int range_residue_selection = 0;
                                if(residue_sequence_number.find("-") != string::npos)
                                {
                                    range_residue_selection = 1;
                                }
                                else
                                    range_residue_selection = 0;
                                switch(range_residue_selection)
                                {
                                    case 0:
                                    {
                                        if(residue_sequence_number.find(residue_name) != string::npos &&
                                                residue_sequence_number.find(residue_name) == 0)
                                        {
                                            vector<string> atom_names = (*it1).second;
                                            if(find(atom_names.begin(), atom_names.end(), "*") == atom_names.end())
                                            {
                                                AtomVector atoms = residue->GetAtoms();
                                                for(vector<string>::iterator it3 = atom_names.begin(); it3 != atom_names.end(); it3++)
                                                {
                                                    string atom_name = *it3;
                                                    int atom_name_search_type = 0;
                                                    if(atom_name.at(0) == '^')
                                                    {
                                                        atom_name = atom_name.substr(1);
                                                        atom_name_search_type = -1;
                                                    }
                                                    else if(atom_name.at(atom_name.size() - 1) == '$')
                                                    {
                                                        atom_name = atom_name.substr(0, atom_name.size() - 2);
                                                        atom_name_search_type = 1;
                                                    }
                                                    else if(atom_name.at(0) == '#')
                                                    {
                                                        atom_name = atom_name.substr(1);
                                                        atom_name_search_type = -2;
                                                    }
                                                    else
                                                        atom_name_search_type = 0;
                                                    for(AtomVector::iterator it4 = atoms.begin(); it4 != atoms.end(); it4++)
                                                    {
                                                        Atom* atom = *it4;
                                                        switch(atom_name_search_type)
                                                        {
                                                            case 0:
                                                            {
                                                                if(atom->GetName().compare(atom_name) == 0)
                                                                    selection.push_back(atom);
                                                                break;
                                                            }
                                                            case 1:
                                                            {
                                                                if(atom->GetName().find(atom_name) != string::npos &&
                                                                        atom->GetName().find(atom_name) == atom->GetName().size() - atom_name.size())
                                                                    selection.push_back(atom);
                                                                break;
                                                            }
                                                            case -1:
                                                            {
                                                                if(atom->GetName().find(atom_name) != string::npos &&
                                                                        atom->GetName().find(atom_name) == 0)
                                                                    selection.push_back(atom);
                                                                break;
                                                            }
                                                            case -2:
                                                            {
                                                                string atom_serial_number = Split(atom->GetId(), "_").at(1);
                                                                int range_selection = 0;
                                                                if(atom_name.find("-") != string::npos)
                                                                {
                                                                    range_selection = 1;

                                                                }
                                                                else
                                                                    range_selection = 0;
                                                                switch(range_selection)
                                                                {
                                                                    case 0:
                                                                    {
                                                                        if(atom_serial_number.find(atom_name) != string::npos &&
                                                                                atom_serial_number.find(atom_name) == 0)
                                                                            selection.push_back(atom);
                                                                        break;
                                                                    }
                                                                    case 1:
                                                                    {
                                                                        string start_index = Split(atom_name, "-").at(0);
                                                                        string end_index = Split(atom_name, "-").at(1);
                                                                        int atom_serial_number_int = ConvertString<int>(atom_serial_number);
                                                                        int start_index_int = ConvertString<int>(start_index);
                                                                        int end_index_int = ConvertString<int>(end_index);
                                                                        if(atom_serial_number_int >= start_index_int && atom_serial_number_int <= end_index_int)
                                                                            selection.push_back(atom);
                                                                        break;
                                                                    }
                                                                }
                                                                break;
                                                            }
                                                        }

                                                    }
                                                }
                                            }
                                            else
                                            {
                                                AtomVector atoms = residue->GetAtoms();
                                                for(AtomVector::iterator it4 = atoms.begin(); it4 != atoms.end(); it4++)
                                                {
                                                    Atom* atom = *it4;
                                                    selection.push_back(atom);
                                                }
                                            }
                                        }
                                        break;
                                    }
                                    case 1:
                                    {
                                        string residue_start_index = Split(residue_name, "-").at(0);
                                        string residue_end_index = Split(residue_name, "-").at(1);
                                        int residue_sequence_number_int = ConvertString<int>(residue_sequence_number);
                                        int residue_start_index_int = ConvertString<int>(residue_start_index);
                                        int residue_end_index_int = ConvertString<int>(residue_end_index);
                                        if(residue_sequence_number_int >= residue_start_index_int && residue_sequence_number_int <= residue_end_index_int)
                                        {
                                            vector<string> atom_names = (*it1).second;
                                            if(find(atom_names.begin(), atom_names.end(), "*") == atom_names.end())
                                            {
                                                AtomVector atoms = residue->GetAtoms();
                                                for(vector<string>::iterator it3 = atom_names.begin(); it3 != atom_names.end(); it3++)
                                                {
                                                    string atom_name = *it3;
                                                    int atom_name_search_type = 0;
                                                    if(atom_name.at(0) == '^')
                                                    {
                                                        atom_name = atom_name.substr(1);
                                                        atom_name_search_type = -1;
                                                    }
                                                    else if(atom_name.at(atom_name.size() - 1) == '$')
                                                    {
                                                        atom_name = atom_name.substr(0, atom_name.size() - 2);
                                                        atom_name_search_type = 1;
                                                    }
                                                    else if(atom_name.at(0) == '#')
                                                    {
                                                        atom_name = atom_name.substr(1);
                                                        atom_name_search_type = -2;
                                                    }
                                                    else
                                                        atom_name_search_type = 0;
                                                    for(AtomVector::iterator it4 = atoms.begin(); it4 != atoms.end(); it4++)
                                                    {
                                                        Atom* atom = *it4;
                                                        switch(atom_name_search_type)
                                                        {
                                                            case 0:
                                                            {
                                                                if(atom->GetName().compare(atom_name) == 0)
                                                                    selection.push_back(atom);
                                                                break;
                                                            }
                                                            case 1:
                                                            {
                                                                if(atom->GetName().find(atom_name) != string::npos &&
                                                                        atom->GetName().find(atom_name) == atom->GetName().size() - atom_name.size())
                                                                    selection.push_back(atom);
                                                                break;
                                                            }
                                                            case -1:
                                                            {
                                                                if(atom->GetName().find(atom_name) != string::npos &&
                                                                        atom->GetName().find(atom_name) == 0)
                                                                    selection.push_back(atom);
                                                                break;
                                                            }
                                                            case -2:
                                                            {
                                                                string atom_serial_number = Split(atom->GetId(), "_").at(1);
                                                                int range_selection = 0;
                                                                if(atom_name.find("-") != string::npos)
                                                                {
                                                                    range_selection = 1;

                                                                }
                                                                else
                                                                    range_selection = 0;
                                                                switch(range_selection)
                                                                {
                                                                    case 0:
                                                                    {
                                                                        if(atom_serial_number.find(atom_name) != string::npos &&
                                                                                atom_serial_number.find(atom_name) == 0)
                                                                            selection.push_back(atom);
                                                                        break;
                                                                    }
                                                                    case 1:
                                                                    {
                                                                        string start_index = Split(atom_name, "-").at(0);
                                                                        string end_index = Split(atom_name, "-").at(1);
                                                                        int atom_serial_number_int = ConvertString<int>(atom_serial_number);
                                                                        int start_index_int = ConvertString<int>(start_index);
                                                                        int end_index_int = ConvertString<int>(end_index);
                                                                        if(atom_serial_number_int >= start_index_int && atom_serial_number_int <= end_index_int)
                                                                            selection.push_back(atom);
                                                                        break;
                                                                    }
                                                                }
                                                                break;
                                                            }
                                                        }

                                                    }
                                                }
                                            }
                                            else
                                            {
                                                AtomVector atoms = residue->GetAtoms();
                                                for(AtomVector::iterator it4 = atoms.begin(); it4 != atoms.end(); it4++)
                                                {
                                                    Atom* atom = *it4;
                                                    selection.push_back(atom);
                                                }
                                            }
                                            break;
                                        }
                                    }
                                }
                                break;
                            }
                        }
                    }
                }
                else
                {
                    for(ResidueVector::iterator it2 = residues_of_assembly.begin(); it2 != residues_of_assembly.end(); it2++)
                    {
                        Residue* residue = *it2;
                        vector<string> atom_names = (*it1).second;
                        if(find(atom_names.begin(), atom_names.end(), "*") == atom_names.end())
                        {
                            AtomVector atoms = residue->GetAtoms();
                            for(vector<string>::iterator it3 = atom_names.begin(); it3 != atom_names.end(); it3++)
                            {
                                string atom_name = *it3;
                                int atom_name_search_type = 0;
                                if(atom_name.at(0) == '^')
                                {
                                    atom_name = atom_name.substr(1);
                                    atom_name_search_type = -1;
                                }
                                else if(atom_name.at(atom_name.size() - 1) == '$')
                                {
                                    atom_name = atom_name.substr(0, atom_name.size() - 2);
                                    atom_name_search_type = 1;
                                }
                                else if(atom_name.at(0) == '#')
                                {
                                    atom_name = atom_name.substr(1);
                                    atom_name_search_type = -2;
                                }
                                else
                                    atom_name_search_type = 0;
                                for(AtomVector::iterator it4 = atoms.begin(); it4 != atoms.end(); it4++)
                                {
                                    Atom* atom = *it4;
                                    switch(atom_name_search_type)
                                    {
                                        case 0:
                                        {
                                            if(atom->GetName().compare(atom_name) == 0)
                                                selection.push_back(atom);
                                            break;
                                        }
                                        case 1:
                                        {
                                            if(atom->GetName().find(atom_name) != string::npos &&
                                                    atom->GetName().find(atom_name) == atom->GetName().size() - atom_name.size())
                                                selection.push_back(atom);
                                            break;
                                        }
                                        case -1:
                                        {
                                            if(atom->GetName().find(atom_name) != string::npos &&
                                                    atom->GetName().find(atom_name) == 0)
                                                selection.push_back(atom);
                                            break;
                                        }
                                        case -2:
                                        {
                                            string atom_serial_number = Split(atom->GetId(), "_").at(1);
                                            int range_selection = 0;
                                            if(atom_name.find("-") != string::npos)
                                            {
                                                range_selection = 1;

                                            }
                                            else
                                                range_selection = 0;
                                            switch(range_selection)
                                            {
                                                case 0:
                                                {
                                                    if(atom_serial_number.find(atom_name) != string::npos &&
                                                            atom_serial_number.find(atom_name) == 0)
                                                        selection.push_back(atom);
                                                    break;
                                                }
                                                case 1:
                                                {
                                                    string start_index = Split(atom_name, "-").at(0);
                                                    string end_index = Split(atom_name, "-").at(1);
                                                    int atom_serial_number_int = ConvertString<int>(atom_serial_number);
                                                    int start_index_int = ConvertString<int>(start_index);
                                                    int end_index_int = ConvertString<int>(end_index);
                                                    if(atom_serial_number_int >= start_index_int && atom_serial_number_int <= end_index_int)
                                                        selection.push_back(atom);
                                                    break;
                                                }
                                            }
                                            break;
                                        }
                                    }

                                }
                            }
                        }
                        else
                        {
                            AtomVector atoms = residue->GetAtoms();
                            for(AtomVector::iterator it4 = atoms.begin(); it4 != atoms.end(); it4++)
                            {
                                Atom* atom = *it4;
                                selection.push_back(atom);
                            }
                        }
                    }
                }

            }

        }
        else
        {
            int start_index = key.find("*");
            if(start_index == 0)
            {
                for(HierarchicalContainmentMap::iterator it5 = hierarchical_map.begin(); it5 != hierarchical_map.end(); it5++)
                {
                    ResidueVector residues_of_assembly = (*it5).second;
                    map<string, vector<string> > value = (*it).second;
                    for(map<string, vector<string> >::iterator it1 = value.begin(); it1 != value.end(); it1++)
                    {
                        string residue_name = (*it1).first;
                        int residue_name_search_type = 0;
                        if(residue_name.at(0) == '^')
                        {
                            residue_name = residue_name.substr(1);
                            residue_name_search_type = -1;
                        }
                        else if(residue_name.at(residue_name.size() - 1) == '$')
                        {
                            residue_name = residue_name.substr(0, residue_name.size() - 2);
                            residue_name_search_type = 1;
                        }
                        else if(residue_name.at(0) == '#')
                        {
                            residue_name = residue_name.substr(1);
                            residue_name_search_type = -2;
                        }
                        else
                            residue_name_search_type = 0;
                        if(residue_name.find("*") == string::npos)
                        {
                            for(ResidueVector::iterator it2 = residues_of_assembly.begin(); it2 != residues_of_assembly.end(); it2++)
                            {
                                Residue* residue = *it2;
                                switch(residue_name_search_type)
                                {
                                    case 0:  /// Search in residue set by matching the whole name of the residue
                                    {
                                        if(residue->GetName().compare(residue_name) == 0)
                                        {
                                            vector<string> atom_names = (*it1).second;
                                            if(find(atom_names.begin(), atom_names.end(), "*") == atom_names.end())
                                            {
                                                AtomVector atoms = residue->GetAtoms();
                                                for(vector<string>::iterator it3 = atom_names.begin(); it3 != atom_names.end(); it3++)
                                                {
                                                    string atom_name = *it3;
                                                    int atom_name_search_type = 0;
                                                    if(atom_name.at(0) == '^')
                                                    {
                                                        atom_name = atom_name.substr(1);
                                                        atom_name_search_type = -1;
                                                    }
                                                    else if(atom_name.at(atom_name.size() - 1) == '$')
                                                    {
                                                        atom_name = atom_name.substr(0, atom_name.size() - 2);
                                                        atom_name_search_type = 1;
                                                    }
                                                    else if(atom_name.at(0) == '#')
                                                    {
                                                        atom_name = atom_name.substr(1);
                                                        atom_name_search_type = -2;
                                                    }
                                                    else
                                                        atom_name_search_type = 0;
                                                    for(AtomVector::iterator it4 = atoms.begin(); it4 != atoms.end(); it4++)
                                                    {
                                                        Atom* atom = *it4;
                                                        switch(atom_name_search_type)
                                                        {
                                                            case 0:
                                                            {
                                                                if(atom->GetName().compare(atom_name) == 0)
                                                                    selection.push_back(atom);
                                                                break;
                                                            }
                                                            case 1:
                                                            {
                                                                if(atom->GetName().find(atom_name) != string::npos &&
                                                                        atom->GetName().find(atom_name) == atom->GetName().size() - atom_name.size())
                                                                    selection.push_back(atom);
                                                                break;
                                                            }
                                                            case -1:
                                                            {
                                                                if(atom->GetName().find(atom_name) != string::npos &&
                                                                        atom->GetName().find(atom_name) == 0)
                                                                    selection.push_back(atom);
                                                                break;
                                                            }
                                                            case -2:
                                                            {
                                                                string atom_serial_number = Split(atom->GetId(), "_").at(1);
                                                                int range_selection = 0;
                                                                if(atom_name.find("-") != string::npos)
                                                                {
                                                                    range_selection = 1;

                                                                }
                                                                else
                                                                    range_selection = 0;
                                                                switch(range_selection)
                                                                {
                                                                    case 0:
                                                                    {
                                                                        if(atom_serial_number.find(atom_name) != string::npos &&
                                                                                atom_serial_number.find(atom_name) == 0)
                                                                            selection.push_back(atom);
                                                                        break;
                                                                    }
                                                                    case 1:
                                                                    {
                                                                        string start_index = Split(atom_name, "-").at(0);
                                                                        string end_index = Split(atom_name, "-").at(1);
                                                                        int atom_serial_number_int = ConvertString<int>(atom_serial_number);
                                                                        int start_index_int = ConvertString<int>(start_index);
                                                                        int end_index_int = ConvertString<int>(end_index);
                                                                        if(atom_serial_number_int >= start_index_int && atom_serial_number_int <= end_index_int)
                                                                            selection.push_back(atom);
                                                                        break;
                                                                    }
                                                                }
                                                                break;
                                                            }
                                                        }

                                                    }
                                                }
                                            }
                                            else
                                            {
                                                AtomVector atoms = residue->GetAtoms();
                                                for(AtomVector::iterator it4 = atoms.begin(); it4 != atoms.end(); it4++)
                                                {
                                                    Atom* atom = *it4;
                                                    selection.push_back(atom);
                                                }
                                            }
                                        }
                                        break;
                                    }
                                    case 1:   /// Searching the residue set by maching the end of the residue names
                                    {
                                        if(residue->GetName().find(residue_name) != string::npos &&
                                                residue->GetName().find(residue_name) == residue->GetName().size() - residue_name.size())
                                        {
                                            vector<string> atom_names = (*it1).second;
                                            if(find(atom_names.begin(), atom_names.end(), "*") == atom_names.end())
                                            {
                                                AtomVector atoms = residue->GetAtoms();
                                                for(vector<string>::iterator it3 = atom_names.begin(); it3 != atom_names.end(); it3++)
                                                {
                                                    string atom_name = *it3;
                                                    int atom_name_search_type = 0;
                                                    if(atom_name.at(0) == '^')
                                                    {
                                                        atom_name = atom_name.substr(1);
                                                        atom_name_search_type = -1;
                                                    }
                                                    else if(atom_name.at(atom_name.size() - 1) == '$')
                                                    {
                                                        atom_name = atom_name.substr(0, atom_name.size() - 2);
                                                        atom_name_search_type = 1;
                                                    }
                                                    else if(atom_name.at(0) == '#')
                                                    {
                                                        atom_name = atom_name.substr(1);
                                                        atom_name_search_type = -2;
                                                    }
                                                    else
                                                        atom_name_search_type = 0;
                                                    for(AtomVector::iterator it4 = atoms.begin(); it4 != atoms.end(); it4++)
                                                    {
                                                        Atom* atom = *it4;
                                                        switch(atom_name_search_type)
                                                        {
                                                            case 0:
                                                            {
                                                                if(atom->GetName().compare(atom_name) == 0)
                                                                    selection.push_back(atom);
                                                                break;
                                                            }
                                                            case 1:
                                                            {
                                                                if(atom->GetName().find(atom_name) != string::npos &&
                                                                        atom->GetName().find(atom_name) == atom->GetName().size() - atom_name.size())
                                                                    selection.push_back(atom);
                                                                break;
                                                            }
                                                            case -1:
                                                            {
                                                                if(atom->GetName().find(atom_name) != string::npos &&
                                                                        atom->GetName().find(atom_name) == 0)
                                                                    selection.push_back(atom);
                                                                break;
                                                            }
                                                            case -2:
                                                            {
                                                                string atom_serial_number = Split(atom->GetId(), "_").at(1);
                                                                int range_selection = 0;
                                                                if(atom_name.find("-") != string::npos)
                                                                {
                                                                    range_selection = 1;

                                                                }
                                                                else
                                                                    range_selection = 0;
                                                                switch(range_selection)
                                                                {
                                                                    case 0:
                                                                    {
                                                                        if(atom_serial_number.find(atom_name) != string::npos &&
                                                                                atom_serial_number.find(atom_name) == 0)
                                                                            selection.push_back(atom);
                                                                        break;
                                                                    }
                                                                    case 1:
                                                                    {
                                                                        string start_index = Split(atom_name, "-").at(0);
                                                                        string end_index = Split(atom_name, "-").at(1);
                                                                        int atom_serial_number_int = ConvertString<int>(atom_serial_number);
                                                                        int start_index_int = ConvertString<int>(start_index);
                                                                        int end_index_int = ConvertString<int>(end_index);
                                                                        if(atom_serial_number_int >= start_index_int && atom_serial_number_int <= end_index_int)
                                                                            selection.push_back(atom);
                                                                        break;
                                                                    }
                                                                }
                                                                break;
                                                            }
                                                        }

                                                    }
                                                }
                                            }
                                            else
                                            {
                                                AtomVector atoms = residue->GetAtoms();
                                                for(AtomVector::iterator it4 = atoms.begin(); it4 != atoms.end(); it4++)
                                                {
                                                    Atom* atom = *it4;
                                                    selection.push_back(atom);
                                                }
                                            }
                                        }
                                        break;
                                    }
                                    case -1:  /// Searching the residue set by matching the begining of the residue names
                                    {
                                        if(residue->GetName().find(residue_name) != string::npos &&
                                                residue->GetName().find(residue_name) == 0)
                                        {
                                            vector<string> atom_names = (*it1).second;
                                            if(find(atom_names.begin(), atom_names.end(), "*") == atom_names.end())
                                            {
                                                AtomVector atoms = residue->GetAtoms();
                                                for(vector<string>::iterator it3 = atom_names.begin(); it3 != atom_names.end(); it3++)
                                                {
                                                    string atom_name = *it3;
                                                    int atom_name_search_type = 0;
                                                    if(atom_name.at(0) == '^')
                                                    {
                                                        atom_name = atom_name.substr(1);
                                                        atom_name_search_type = -1;
                                                    }
                                                    else if(atom_name.at(atom_name.size() - 1) == '$')
                                                    {
                                                        atom_name = atom_name.substr(0, atom_name.size() - 2);
                                                        atom_name_search_type = 1;
                                                    }
                                                    else if(atom_name.at(0) == '#')
                                                    {
                                                        atom_name = atom_name.substr(1);
                                                        atom_name_search_type = -2;
                                                    }
                                                    else
                                                        atom_name_search_type = 0;
                                                    for(AtomVector::iterator it4 = atoms.begin(); it4 != atoms.end(); it4++)
                                                    {
                                                        Atom* atom = *it4;
                                                        switch(atom_name_search_type)
                                                        {
                                                            case 0:
                                                            {
                                                                if(atom->GetName().compare(atom_name) == 0)
                                                                    selection.push_back(atom);
                                                                break;
                                                            }
                                                            case 1:
                                                            {
                                                                if(atom->GetName().find(atom_name) != string::npos &&
                                                                        atom->GetName().find(atom_name) == atom->GetName().size() - atom_name.size())
                                                                    selection.push_back(atom);
                                                                break;
                                                            }
                                                            case -1:
                                                            {
                                                                if(atom->GetName().find(atom_name) != string::npos &&
                                                                        atom->GetName().find(atom_name) == 0)
                                                                    selection.push_back(atom);
                                                                break;
                                                            }
                                                            case -2:
                                                            {
                                                                string atom_serial_number = Split(atom->GetId(), "_").at(1);
                                                                int range_selection = 0;
                                                                if(atom_name.find("-") != string::npos)
                                                                {
                                                                    range_selection = 1;

                                                                }
                                                                else
                                                                    range_selection = 0;
                                                                switch(range_selection)
                                                                {
                                                                    case 0:
                                                                    {
                                                                        if(atom_serial_number.find(atom_name) != string::npos &&
                                                                                atom_serial_number.find(atom_name) == 0)
                                                                            selection.push_back(atom);
                                                                        break;
                                                                    }
                                                                    case 1:
                                                                    {
                                                                        string start_index = Split(atom_name, "-").at(0);
                                                                        string end_index = Split(atom_name, "-").at(1);
                                                                        int atom_serial_number_int = ConvertString<int>(atom_serial_number);
                                                                        int start_index_int = ConvertString<int>(start_index);
                                                                        int end_index_int = ConvertString<int>(end_index);
                                                                        if(atom_serial_number_int >= start_index_int && atom_serial_number_int <= end_index_int)
                                                                            selection.push_back(atom);
                                                                        break;
                                                                    }
                                                                }
                                                                break;
                                                            }
                                                        }

                                                    }
                                                }
                                            }
                                            else
                                            {
                                                AtomVector atoms = residue->GetAtoms();
                                                for(AtomVector::iterator it4 = atoms.begin(); it4 != atoms.end(); it4++)
                                                {
                                                    Atom* atom = *it4;
                                                    selection.push_back(atom);
                                                }
                                            }
                                        }
                                        break;
                                    }
                                    case -2:  /// Searching the residue set by matching the residue sequence number
                                    {
                                        string residue_sequence_number = Split(residue->GetId(), "_").at(2);
                                        int range_residue_selection = 0;
                                        if(residue_sequence_number.find("-") != string::npos)
                                        {
                                            range_residue_selection = 1;
                                        }
                                        else
                                            range_residue_selection = 0;
                                        switch(range_residue_selection)
                                        {
                                            case 0:
                                            {
                                                if(residue_sequence_number.find(residue_name) != string::npos &&
                                                        residue_sequence_number.find(residue_name) == 0)
                                                {
                                                    vector<string> atom_names = (*it1).second;
                                                    if(find(atom_names.begin(), atom_names.end(), "*") == atom_names.end())
                                                    {
                                                        AtomVector atoms = residue->GetAtoms();
                                                        for(vector<string>::iterator it3 = atom_names.begin(); it3 != atom_names.end(); it3++)
                                                        {
                                                            string atom_name = *it3;
                                                            int atom_name_search_type = 0;
                                                            if(atom_name.at(0) == '^')
                                                            {
                                                                atom_name = atom_name.substr(1);
                                                                atom_name_search_type = -1;
                                                            }
                                                            else if(atom_name.at(atom_name.size() - 1) == '$')
                                                            {
                                                                atom_name = atom_name.substr(0, atom_name.size() - 2);
                                                                atom_name_search_type = 1;
                                                            }
                                                            else if(atom_name.at(0) == '#')
                                                            {
                                                                atom_name = atom_name.substr(1);
                                                                atom_name_search_type = -2;
                                                            }
                                                            else
                                                                atom_name_search_type = 0;
                                                            for(AtomVector::iterator it4 = atoms.begin(); it4 != atoms.end(); it4++)
                                                            {
                                                                Atom* atom = *it4;
                                                                switch(atom_name_search_type)
                                                                {
                                                                    case 0:
                                                                    {
                                                                        if(atom->GetName().compare(atom_name) == 0)
                                                                            selection.push_back(atom);
                                                                        break;
                                                                    }
                                                                    case 1:
                                                                    {
                                                                        if(atom->GetName().find(atom_name) != string::npos &&
                                                                                atom->GetName().find(atom_name) == atom->GetName().size() - atom_name.size())
                                                                            selection.push_back(atom);
                                                                        break;
                                                                    }
                                                                    case -1:
                                                                    {
                                                                        if(atom->GetName().find(atom_name) != string::npos &&
                                                                                atom->GetName().find(atom_name) == 0)
                                                                            selection.push_back(atom);
                                                                        break;
                                                                    }
                                                                    case -2:
                                                                    {
                                                                        string atom_serial_number = Split(atom->GetId(), "_").at(1);
                                                                        int range_selection = 0;
                                                                        if(atom_name.find("-") != string::npos)
                                                                        {
                                                                            range_selection = 1;

                                                                        }
                                                                        else
                                                                            range_selection = 0;
                                                                        switch(range_selection)
                                                                        {
                                                                            case 0:
                                                                            {
                                                                                if(atom_serial_number.find(atom_name) != string::npos &&
                                                                                        atom_serial_number.find(atom_name) == 0)
                                                                                    selection.push_back(atom);
                                                                                break;
                                                                            }
                                                                            case 1:
                                                                            {
                                                                                string start_index = Split(atom_name, "-").at(0);
                                                                                string end_index = Split(atom_name, "-").at(1);
                                                                                int atom_serial_number_int = ConvertString<int>(atom_serial_number);
                                                                                int start_index_int = ConvertString<int>(start_index);
                                                                                int end_index_int = ConvertString<int>(end_index);
                                                                                if(atom_serial_number_int >= start_index_int && atom_serial_number_int <= end_index_int)
                                                                                    selection.push_back(atom);
                                                                                break;
                                                                            }
                                                                        }
                                                                        break;
                                                                    }
                                                                }

                                                            }
                                                        }
                                                    }
                                                    else
                                                    {
                                                        AtomVector atoms = residue->GetAtoms();
                                                        for(AtomVector::iterator it4 = atoms.begin(); it4 != atoms.end(); it4++)
                                                        {
                                                            Atom* atom = *it4;
                                                            selection.push_back(atom);
                                                        }
                                                    }
                                                }
                                                break;
                                            }
                                            case 1:
                                            {
                                                string residue_start_index = Split(residue_name, "-").at(0);
                                                string residue_end_index = Split(residue_name, "-").at(1);
                                                int residue_sequence_number_int = ConvertString<int>(residue_sequence_number);
                                                int residue_start_index_int = ConvertString<int>(residue_start_index);
                                                int residue_end_index_int = ConvertString<int>(residue_end_index);
                                                if(residue_sequence_number_int >= residue_start_index_int && residue_sequence_number_int <= residue_end_index_int)
                                                {
                                                    vector<string> atom_names = (*it1).second;
                                                    if(find(atom_names.begin(), atom_names.end(), "*") == atom_names.end())
                                                    {
                                                        AtomVector atoms = residue->GetAtoms();
                                                        for(vector<string>::iterator it3 = atom_names.begin(); it3 != atom_names.end(); it3++)
                                                        {
                                                            string atom_name = *it3;
                                                            int atom_name_search_type = 0;
                                                            if(atom_name.at(0) == '^')
                                                            {
                                                                atom_name = atom_name.substr(1);
                                                                atom_name_search_type = -1;
                                                            }
                                                            else if(atom_name.at(atom_name.size() - 1) == '$')
                                                            {
                                                                atom_name = atom_name.substr(0, atom_name.size() - 2);
                                                                atom_name_search_type = 1;
                                                            }
                                                            else if(atom_name.at(0) == '#')
                                                            {
                                                                atom_name = atom_name.substr(1);
                                                                atom_name_search_type = -2;
                                                            }
                                                            else
                                                                atom_name_search_type = 0;
                                                            for(AtomVector::iterator it4 = atoms.begin(); it4 != atoms.end(); it4++)
                                                            {
                                                                Atom* atom = *it4;
                                                                switch(atom_name_search_type)
                                                                {
                                                                    case 0:
                                                                    {
                                                                        if(atom->GetName().compare(atom_name) == 0)
                                                                            selection.push_back(atom);
                                                                        break;
                                                                    }
                                                                    case 1:
                                                                    {
                                                                        if(atom->GetName().find(atom_name) != string::npos &&
                                                                                atom->GetName().find(atom_name) == atom->GetName().size() - atom_name.size())
                                                                            selection.push_back(atom);
                                                                        break;
                                                                    }
                                                                    case -1:
                                                                    {
                                                                        if(atom->GetName().find(atom_name) != string::npos &&
                                                                                atom->GetName().find(atom_name) == 0)
                                                                            selection.push_back(atom);
                                                                        break;
                                                                    }
                                                                    case -2:
                                                                    {
                                                                        string atom_serial_number = Split(atom->GetId(), "_").at(1);
                                                                        int range_selection = 0;
                                                                        if(atom_name.find("-") != string::npos)
                                                                        {
                                                                            range_selection = 1;

                                                                        }
                                                                        else
                                                                            range_selection = 0;
                                                                        switch(range_selection)
                                                                        {
                                                                            case 0:
                                                                            {
                                                                                if(atom_serial_number.find(atom_name) != string::npos &&
                                                                                        atom_serial_number.find(atom_name) == 0)
                                                                                    selection.push_back(atom);
                                                                                break;
                                                                            }
                                                                            case 1:
                                                                            {
                                                                                string start_index = Split(atom_name, "-").at(0);
                                                                                string end_index = Split(atom_name, "-").at(1);
                                                                                int atom_serial_number_int = ConvertString<int>(atom_serial_number);
                                                                                int start_index_int = ConvertString<int>(start_index);
                                                                                int end_index_int = ConvertString<int>(end_index);
                                                                                if(atom_serial_number_int >= start_index_int && atom_serial_number_int <= end_index_int)
                                                                                    selection.push_back(atom);
                                                                                break;
                                                                            }
                                                                        }
                                                                        break;
                                                                    }
                                                                }

                                                            }
                                                        }
                                                    }
                                                    else
                                                    {
                                                        AtomVector atoms = residue->GetAtoms();
                                                        for(AtomVector::iterator it4 = atoms.begin(); it4 != atoms.end(); it4++)
                                                        {
                                                            Atom* atom = *it4;
                                                            selection.push_back(atom);
                                                        }
                                                    }
                                                    break;
                                                }
                                            }
                                        }
                                        break;
                                    }
                                }
                            }
                        }
                        else
                        {
                            for(ResidueVector::iterator it2 = residues_of_assembly.begin(); it2 != residues_of_assembly.end(); it2++)
                            {
                                Residue* residue = *it2;
                                vector<string> atom_names = (*it1).second;
                                if(find(atom_names.begin(), atom_names.end(), "*") == atom_names.end())
                                {
                                    AtomVector atoms = residue->GetAtoms();
                                    for(vector<string>::iterator it3 = atom_names.begin(); it3 != atom_names.end(); it3++)
                                    {
                                        string atom_name = *it3;
                                        int atom_name_search_type = 0;
                                        if(atom_name.at(0) == '^')
                                        {
                                            atom_name = atom_name.substr(1);
                                            atom_name_search_type = -1;
                                        }
                                        else if(atom_name.at(atom_name.size() - 1) == '$')
                                        {
                                            atom_name = atom_name.substr(0, atom_name.size() - 2);
                                            atom_name_search_type = 1;
                                        }
                                        else if(atom_name.at(0) == '#')
                                        {
                                            atom_name = atom_name.substr(1);
                                            atom_name_search_type = -2;
                                        }
                                        else
                                            atom_name_search_type = 0;
                                        for(AtomVector::iterator it4 = atoms.begin(); it4 != atoms.end(); it4++)
                                        {
                                            Atom* atom = *it4;
                                            switch(atom_name_search_type)
                                            {
                                                case 0:
                                                {
                                                    if(atom->GetName().compare(atom_name) == 0)
                                                        selection.push_back(atom);
                                                    break;
                                                }
                                                case 1:
                                                {
                                                    if(atom->GetName().find(atom_name) != string::npos &&
                                                            atom->GetName().find(atom_name) == atom->GetName().size() - atom_name.size())
                                                        selection.push_back(atom);
                                                    break;
                                                }
                                                case -1:
                                                {
                                                    if(atom->GetName().find(atom_name) != string::npos &&
                                                            atom->GetName().find(atom_name) == 0)
                                                        selection.push_back(atom);
                                                    break;
                                                }
                                                case -2:
                                                {
                                                    string atom_serial_number = Split(atom->GetId(), "_").at(1);
                                                    int range_selection = 0;
                                                    if(atom_name.find("-") != string::npos)
                                                    {
                                                        range_selection = 1;

                                                    }
                                                    else
                                                        range_selection = 0;
                                                    switch(range_selection)
                                                    {
                                                        case 0:
                                                        {
                                                            if(atom_serial_number.find(atom_name) != string::npos &&
                                                                    atom_serial_number.find(atom_name) == 0)
                                                                selection.push_back(atom);
                                                            break;
                                                        }
                                                        case 1:
                                                        {
                                                            string start_index = Split(atom_name, "-").at(0);
                                                            string end_index = Split(atom_name, "-").at(1);
                                                            int atom_serial_number_int = ConvertString<int>(atom_serial_number);
                                                            int start_index_int = ConvertString<int>(start_index);
                                                            int end_index_int = ConvertString<int>(end_index);
                                                            if(atom_serial_number_int >= start_index_int && atom_serial_number_int <= end_index_int)
                                                                selection.push_back(atom);
                                                            break;
                                                        }
                                                    }
                                                    break;
                                                }
                                            }

                                        }
                                    }
                                }
                                else
                                {
                                    AtomVector atoms = residue->GetAtoms();
                                    for(AtomVector::iterator it4 = atoms.begin(); it4 != atoms.end(); it4++)
                                    {
                                        Atom* atom = *it4;
                                        selection.push_back(atom);
                                    }
                                }
                            }
                        }

                    }

                }
            }
            else
            {
                string partial_key = key.substr(0, key.find(".*"));
                for(HierarchicalContainmentMap::iterator it5 = hierarchical_map.begin(); it5 != hierarchical_map.end(); it5++)
                {
                    if((*it5).first.find(partial_key) == 0)
                    {
                        ResidueVector residues_of_assembly = (*it5).second;
                        map<string, vector<string> > value = (*it).second;
                        for(map<string, vector<string> >::iterator it1 = value.begin(); it1 != value.end(); it1++)
                        {
                            string residue_name = (*it1).first;
                            int residue_name_search_type = 0;
                            if(residue_name.at(0) == '^')
                            {
                                residue_name = residue_name.substr(1);
                                residue_name_search_type = -1;
                            }
                            else if(residue_name.at(residue_name.size() - 1) == '$')
                            {
                                residue_name = residue_name.substr(0, residue_name.size() - 2);
                                residue_name_search_type = 1;
                            }
                            else if(residue_name.at(0) == '#')
                            {
                                residue_name = residue_name.substr(1);
                                residue_name_search_type = -2;
                            }
                            else
                                residue_name_search_type = 0;
                            if(residue_name.find("*") == string::npos)
                            {
                                for(ResidueVector::iterator it2 = residues_of_assembly.begin(); it2 != residues_of_assembly.end(); it2++)
                                {
                                    Residue* residue = *it2;
                                    switch(residue_name_search_type)
                                    {
                                        case 0:  /// Search in residue set by matching the whole name of the residue
                                        {
                                            if(residue->GetName().compare(residue_name) == 0)
                                            {
                                                vector<string> atom_names = (*it1).second;
                                                if(find(atom_names.begin(), atom_names.end(), "*") == atom_names.end())
                                                {
                                                    AtomVector atoms = residue->GetAtoms();
                                                    for(vector<string>::iterator it3 = atom_names.begin(); it3 != atom_names.end(); it3++)
                                                    {
                                                        string atom_name = *it3;
                                                        int atom_name_search_type = 0;
                                                        if(atom_name.at(0) == '^')
                                                        {
                                                            atom_name = atom_name.substr(1);
                                                            atom_name_search_type = -1;
                                                        }
                                                        else if(atom_name.at(atom_name.size() - 1) == '$')
                                                        {
                                                            atom_name = atom_name.substr(0, atom_name.size() - 2);
                                                            atom_name_search_type = 1;
                                                        }
                                                        else if(atom_name.at(0) == '#')
                                                        {
                                                            atom_name = atom_name.substr(1);
                                                            atom_name_search_type = -2;
                                                        }
                                                        else
                                                            atom_name_search_type = 0;
                                                        for(AtomVector::iterator it4 = atoms.begin(); it4 != atoms.end(); it4++)
                                                        {
                                                            Atom* atom = *it4;
                                                            switch(atom_name_search_type)
                                                            {
                                                                case 0:
                                                                {
                                                                    if(atom->GetName().compare(atom_name) == 0)
                                                                        selection.push_back(atom);
                                                                    break;
                                                                }
                                                                case 1:
                                                                {
                                                                    if(atom->GetName().find(atom_name) != string::npos &&
                                                                            atom->GetName().find(atom_name) == atom->GetName().size() - atom_name.size())
                                                                        selection.push_back(atom);
                                                                    break;
                                                                }
                                                                case -1:
                                                                {
                                                                    if(atom->GetName().find(atom_name) != string::npos &&
                                                                            atom->GetName().find(atom_name) == 0)
                                                                        selection.push_back(atom);
                                                                    break;
                                                                }
                                                                case -2:
                                                                {
                                                                    string atom_serial_number = Split(atom->GetId(), "_").at(1);
                                                                    int range_selection = 0;
                                                                    if(atom_name.find("-") != string::npos)
                                                                    {
                                                                        range_selection = 1;

                                                                    }
                                                                    else
                                                                        range_selection = 0;
                                                                    switch(range_selection)
                                                                    {
                                                                        case 0:
                                                                        {
                                                                            if(atom_serial_number.find(atom_name) != string::npos &&
                                                                                    atom_serial_number.find(atom_name) == 0)
                                                                                selection.push_back(atom);
                                                                            break;
                                                                        }
                                                                        case 1:
                                                                        {
                                                                            string start_index = Split(atom_name, "-").at(0);
                                                                            string end_index = Split(atom_name, "-").at(1);
                                                                            int atom_serial_number_int = ConvertString<int>(atom_serial_number);
                                                                            int start_index_int = ConvertString<int>(start_index);
                                                                            int end_index_int = ConvertString<int>(end_index);
                                                                            if(atom_serial_number_int >= start_index_int && atom_serial_number_int <= end_index_int)
                                                                                selection.push_back(atom);
                                                                            break;
                                                                        }
                                                                    }
                                                                    break;
                                                                }
                                                            }

                                                        }
                                                    }
                                                }
                                                else
                                                {
                                                    AtomVector atoms = residue->GetAtoms();
                                                    for(AtomVector::iterator it4 = atoms.begin(); it4 != atoms.end(); it4++)
                                                    {
                                                        Atom* atom = *it4;
                                                        selection.push_back(atom);
                                                    }
                                                }
                                            }
                                            break;
                                        }
                                        case 1:   /// Searching the residue set by maching the end of the residue names
                                        {
                                            if(residue->GetName().find(residue_name) != string::npos &&
                                                    residue->GetName().find(residue_name) == residue->GetName().size() - residue_name.size())
                                            {
                                                vector<string> atom_names = (*it1).second;
                                                if(find(atom_names.begin(), atom_names.end(), "*") == atom_names.end())
                                                {
                                                    AtomVector atoms = residue->GetAtoms();
                                                    for(vector<string>::iterator it3 = atom_names.begin(); it3 != atom_names.end(); it3++)
                                                    {
                                                        string atom_name = *it3;
                                                        int atom_name_search_type = 0;
                                                        if(atom_name.at(0) == '^')
                                                        {
                                                            atom_name = atom_name.substr(1);
                                                            atom_name_search_type = -1;
                                                        }
                                                        else if(atom_name.at(atom_name.size() - 1) == '$')
                                                        {
                                                            atom_name = atom_name.substr(0, atom_name.size() - 2);
                                                            atom_name_search_type = 1;
                                                        }
                                                        else if(atom_name.at(0) == '#')
                                                        {
                                                            atom_name = atom_name.substr(1);
                                                            atom_name_search_type = -2;
                                                        }
                                                        else
                                                            atom_name_search_type = 0;
                                                        for(AtomVector::iterator it4 = atoms.begin(); it4 != atoms.end(); it4++)
                                                        {
                                                            Atom* atom = *it4;
                                                            switch(atom_name_search_type)
                                                            {
                                                                case 0:
                                                                {
                                                                    if(atom->GetName().compare(atom_name) == 0)
                                                                        selection.push_back(atom);
                                                                    break;
                                                                }
                                                                case 1:
                                                                {
                                                                    if(atom->GetName().find(atom_name) != string::npos &&
                                                                            atom->GetName().find(atom_name) == atom->GetName().size() - atom_name.size())
                                                                        selection.push_back(atom);
                                                                    break;
                                                                }
                                                                case -1:
                                                                {
                                                                    if(atom->GetName().find(atom_name) != string::npos &&
                                                                            atom->GetName().find(atom_name) == 0)
                                                                        selection.push_back(atom);
                                                                    break;
                                                                }
                                                                case -2:
                                                                {
                                                                    string atom_serial_number = Split(atom->GetId(), "_").at(1);
                                                                    int range_selection = 0;
                                                                    if(atom_name.find("-") != string::npos)
                                                                    {
                                                                        range_selection = 1;

                                                                    }
                                                                    else
                                                                        range_selection = 0;
                                                                    switch(range_selection)
                                                                    {
                                                                        case 0:
                                                                        {
                                                                            if(atom_serial_number.find(atom_name) != string::npos &&
                                                                                    atom_serial_number.find(atom_name) == 0)
                                                                                selection.push_back(atom);
                                                                            break;
                                                                        }
                                                                        case 1:
                                                                        {
                                                                            string start_index = Split(atom_name, "-").at(0);
                                                                            string end_index = Split(atom_name, "-").at(1);
                                                                            int atom_serial_number_int = ConvertString<int>(atom_serial_number);
                                                                            int start_index_int = ConvertString<int>(start_index);
                                                                            int end_index_int = ConvertString<int>(end_index);
                                                                            if(atom_serial_number_int >= start_index_int && atom_serial_number_int <= end_index_int)
                                                                                selection.push_back(atom);
                                                                            break;
                                                                        }
                                                                    }
                                                                    break;
                                                                }
                                                            }

                                                        }
                                                    }
                                                }
                                                else
                                                {
                                                    AtomVector atoms = residue->GetAtoms();
                                                    for(AtomVector::iterator it4 = atoms.begin(); it4 != atoms.end(); it4++)
                                                    {
                                                        Atom* atom = *it4;
                                                        selection.push_back(atom);
                                                    }
                                                }
                                            }
                                            break;
                                        }
                                        case -1:  /// Searching the residue set by matching the begining of the residue names
                                        {
                                            if(residue->GetName().find(residue_name) != string::npos &&
                                                    residue->GetName().find(residue_name) == 0)
                                            {
                                                vector<string> atom_names = (*it1).second;
                                                if(find(atom_names.begin(), atom_names.end(), "*") == atom_names.end())
                                                {
                                                    AtomVector atoms = residue->GetAtoms();
                                                    for(vector<string>::iterator it3 = atom_names.begin(); it3 != atom_names.end(); it3++)
                                                    {
                                                        string atom_name = *it3;
                                                        int atom_name_search_type = 0;
                                                        if(atom_name.at(0) == '^')
                                                        {
                                                            atom_name = atom_name.substr(1);
                                                            atom_name_search_type = -1;
                                                        }
                                                        else if(atom_name.at(atom_name.size() - 1) == '$')
                                                        {
                                                            atom_name = atom_name.substr(0, atom_name.size() - 2);
                                                            atom_name_search_type = 1;
                                                        }
                                                        else if(atom_name.at(0) == '#')
                                                        {
                                                            atom_name = atom_name.substr(1);
                                                            atom_name_search_type = -2;
                                                        }
                                                        else
                                                            atom_name_search_type = 0;
                                                        for(AtomVector::iterator it4 = atoms.begin(); it4 != atoms.end(); it4++)
                                                        {
                                                            Atom* atom = *it4;
                                                            switch(atom_name_search_type)
                                                            {
                                                                case 0:
                                                                {
                                                                    if(atom->GetName().compare(atom_name) == 0)
                                                                        selection.push_back(atom);
                                                                    break;
                                                                }
                                                                case 1:
                                                                {
                                                                    if(atom->GetName().find(atom_name) != string::npos &&
                                                                            atom->GetName().find(atom_name) == atom->GetName().size() - atom_name.size())
                                                                        selection.push_back(atom);
                                                                    break;
                                                                }
                                                                case -1:
                                                                {
                                                                    if(atom->GetName().find(atom_name) != string::npos &&
                                                                            atom->GetName().find(atom_name) == 0)
                                                                        selection.push_back(atom);
                                                                    break;
                                                                }
                                                                case -2:
                                                                {
                                                                    string atom_serial_number = Split(atom->GetId(), "_").at(1);
                                                                    int range_selection = 0;
                                                                    if(atom_name.find("-") != string::npos)
                                                                    {
                                                                        range_selection = 1;

                                                                    }
                                                                    else
                                                                        range_selection = 0;
                                                                    switch(range_selection)
                                                                    {
                                                                        case 0:
                                                                        {
                                                                            if(atom_serial_number.find(atom_name) != string::npos &&
                                                                                    atom_serial_number.find(atom_name) == 0)
                                                                                selection.push_back(atom);
                                                                            break;
                                                                        }
                                                                        case 1:
                                                                        {
                                                                            string start_index = Split(atom_name, "-").at(0);
                                                                            string end_index = Split(atom_name, "-").at(1);
                                                                            int atom_serial_number_int = ConvertString<int>(atom_serial_number);
                                                                            int start_index_int = ConvertString<int>(start_index);
                                                                            int end_index_int = ConvertString<int>(end_index);
                                                                            if(atom_serial_number_int >= start_index_int && atom_serial_number_int <= end_index_int)
                                                                                selection.push_back(atom);
                                                                            break;
                                                                        }
                                                                    }
                                                                    break;
                                                                }
                                                            }

                                                        }
                                                    }
                                                }
                                                else
                                                {
                                                    AtomVector atoms = residue->GetAtoms();
                                                    for(AtomVector::iterator it4 = atoms.begin(); it4 != atoms.end(); it4++)
                                                    {
                                                        Atom* atom = *it4;
                                                        selection.push_back(atom);
                                                    }
                                                }
                                            }
                                            break;
                                        }
                                        case -2:  /// Searching the residue set by matching the residue sequence number
                                        {
                                            string residue_sequence_number = Split(residue->GetId(), "_").at(2);
                                            int range_residue_selection = 0;
                                            if(residue_sequence_number.find("-") != string::npos)
                                            {
                                                range_residue_selection = 1;
                                            }
                                            else
                                                range_residue_selection = 0;
                                            switch(range_residue_selection)
                                            {
                                                case 0:
                                                {
                                                    if(residue_sequence_number.find(residue_name) != string::npos &&
                                                            residue_sequence_number.find(residue_name) == 0)
                                                    {
                                                        vector<string> atom_names = (*it1).second;
                                                        if(find(atom_names.begin(), atom_names.end(), "*") == atom_names.end())
                                                        {
                                                            AtomVector atoms = residue->GetAtoms();
                                                            for(vector<string>::iterator it3 = atom_names.begin(); it3 != atom_names.end(); it3++)
                                                            {
                                                                string atom_name = *it3;
                                                                int atom_name_search_type = 0;
                                                                if(atom_name.at(0) == '^')
                                                                {
                                                                    atom_name = atom_name.substr(1);
                                                                    atom_name_search_type = -1;
                                                                }
                                                                else if(atom_name.at(atom_name.size() - 1) == '$')
                                                                {
                                                                    atom_name = atom_name.substr(0, atom_name.size() - 2);
                                                                    atom_name_search_type = 1;
                                                                }
                                                                else if(atom_name.at(0) == '#')
                                                                {
                                                                    atom_name = atom_name.substr(1);
                                                                    atom_name_search_type = -2;
                                                                }
                                                                else
                                                                    atom_name_search_type = 0;
                                                                for(AtomVector::iterator it4 = atoms.begin(); it4 != atoms.end(); it4++)
                                                                {
                                                                    Atom* atom = *it4;
                                                                    switch(atom_name_search_type)
                                                                    {
                                                                        case 0:
                                                                        {
                                                                            if(atom->GetName().compare(atom_name) == 0)
                                                                                selection.push_back(atom);
                                                                            break;
                                                                        }
                                                                        case 1:
                                                                        {
                                                                            if(atom->GetName().find(atom_name) != string::npos &&
                                                                                    atom->GetName().find(atom_name) == atom->GetName().size() - atom_name.size())
                                                                                selection.push_back(atom);
                                                                            break;
                                                                        }
                                                                        case -1:
                                                                        {
                                                                            if(atom->GetName().find(atom_name) != string::npos &&
                                                                                    atom->GetName().find(atom_name) == 0)
                                                                                selection.push_back(atom);
                                                                            break;
                                                                        }
                                                                        case -2:
                                                                        {
                                                                            string atom_serial_number = Split(atom->GetId(), "_").at(1);
                                                                            int range_selection = 0;
                                                                            if(atom_name.find("-") != string::npos)
                                                                            {
                                                                                range_selection = 1;

                                                                            }
                                                                            else
                                                                                range_selection = 0;
                                                                            switch(range_selection)
                                                                            {
                                                                                case 0:
                                                                                {
                                                                                    if(atom_serial_number.find(atom_name) != string::npos &&
                                                                                            atom_serial_number.find(atom_name) == 0)
                                                                                        selection.push_back(atom);
                                                                                    break;
                                                                                }
                                                                                case 1:
                                                                                {
                                                                                    string start_index = Split(atom_name, "-").at(0);
                                                                                    string end_index = Split(atom_name, "-").at(1);
                                                                                    int atom_serial_number_int = ConvertString<int>(atom_serial_number);
                                                                                    int start_index_int = ConvertString<int>(start_index);
                                                                                    int end_index_int = ConvertString<int>(end_index);
                                                                                    if(atom_serial_number_int >= start_index_int && atom_serial_number_int <= end_index_int)
                                                                                        selection.push_back(atom);
                                                                                    break;
                                                                                }
                                                                            }
                                                                            break;
                                                                        }
                                                                    }

                                                                }
                                                            }
                                                        }
                                                        else
                                                        {
                                                            AtomVector atoms = residue->GetAtoms();
                                                            for(AtomVector::iterator it4 = atoms.begin(); it4 != atoms.end(); it4++)
                                                            {
                                                                Atom* atom = *it4;
                                                                selection.push_back(atom);
                                                            }
                                                        }
                                                    }
                                                    break;
                                                }
                                                case 1:
                                                {
                                                    string residue_start_index = Split(residue_name, "-").at(0);
                                                    string residue_end_index = Split(residue_name, "-").at(1);
                                                    int residue_sequence_number_int = ConvertString<int>(residue_sequence_number);
                                                    int residue_start_index_int = ConvertString<int>(residue_start_index);
                                                    int residue_end_index_int = ConvertString<int>(residue_end_index);
                                                    if(residue_sequence_number_int >= residue_start_index_int && residue_sequence_number_int <= residue_end_index_int)
                                                    {
                                                        vector<string> atom_names = (*it1).second;
                                                        if(find(atom_names.begin(), atom_names.end(), "*") == atom_names.end())
                                                        {
                                                            AtomVector atoms = residue->GetAtoms();
                                                            for(vector<string>::iterator it3 = atom_names.begin(); it3 != atom_names.end(); it3++)
                                                            {
                                                                string atom_name = *it3;
                                                                int atom_name_search_type = 0;
                                                                if(atom_name.at(0) == '^')
                                                                {
                                                                    atom_name = atom_name.substr(1);
                                                                    atom_name_search_type = -1;
                                                                }
                                                                else if(atom_name.at(atom_name.size() - 1) == '$')
                                                                {
                                                                    atom_name = atom_name.substr(0, atom_name.size() - 2);
                                                                    atom_name_search_type = 1;
                                                                }
                                                                else if(atom_name.at(0) == '#')
                                                                {
                                                                    atom_name = atom_name.substr(1);
                                                                    atom_name_search_type = -2;
                                                                }
                                                                else
                                                                    atom_name_search_type = 0;
                                                                for(AtomVector::iterator it4 = atoms.begin(); it4 != atoms.end(); it4++)
                                                                {
                                                                    Atom* atom = *it4;
                                                                    switch(atom_name_search_type)
                                                                    {
                                                                        case 0:
                                                                        {
                                                                            if(atom->GetName().compare(atom_name) == 0)
                                                                                selection.push_back(atom);
                                                                            break;
                                                                        }
                                                                        case 1:
                                                                        {
                                                                            if(atom->GetName().find(atom_name) != string::npos &&
                                                                                    atom->GetName().find(atom_name) == atom->GetName().size() - atom_name.size())
                                                                                selection.push_back(atom);
                                                                            break;
                                                                        }
                                                                        case -1:
                                                                        {
                                                                            if(atom->GetName().find(atom_name) != string::npos &&
                                                                                    atom->GetName().find(atom_name) == 0)
                                                                                selection.push_back(atom);
                                                                            break;
                                                                        }
                                                                        case -2:
                                                                        {
                                                                            string atom_serial_number = Split(atom->GetId(), "_").at(1);
                                                                            int range_selection = 0;
                                                                            if(atom_name.find("-") != string::npos)
                                                                            {
                                                                                range_selection = 1;

                                                                            }
                                                                            else
                                                                                range_selection = 0;
                                                                            switch(range_selection)
                                                                            {
                                                                                case 0:
                                                                                {
                                                                                    if(atom_serial_number.find(atom_name) != string::npos &&
                                                                                            atom_serial_number.find(atom_name) == 0)
                                                                                        selection.push_back(atom);
                                                                                    break;
                                                                                }
                                                                                case 1:
                                                                                {
                                                                                    string start_index = Split(atom_name, "-").at(0);
                                                                                    string end_index = Split(atom_name, "-").at(1);
                                                                                    int atom_serial_number_int = ConvertString<int>(atom_serial_number);
                                                                                    int start_index_int = ConvertString<int>(start_index);
                                                                                    int end_index_int = ConvertString<int>(end_index);
                                                                                    if(atom_serial_number_int >= start_index_int && atom_serial_number_int <= end_index_int)
                                                                                        selection.push_back(atom);
                                                                                    break;
                                                                                }
                                                                            }
                                                                            break;
                                                                        }
                                                                    }

                                                                }
                                                            }
                                                        }
                                                        else
                                                        {
                                                            AtomVector atoms = residue->GetAtoms();
                                                            for(AtomVector::iterator it4 = atoms.begin(); it4 != atoms.end(); it4++)
                                                            {
                                                                Atom* atom = *it4;
                                                                selection.push_back(atom);
                                                            }
                                                        }
                                                        break;
                                                    }
                                                }
                                            }
                                            break;
                                        }
                                    }
                                }
                            }
                            else
                            {
                                for(ResidueVector::iterator it2 = residues_of_assembly.begin(); it2 != residues_of_assembly.end(); it2++)
                                {
                                    Residue* residue = *it2;
                                    vector<string> atom_names = (*it1).second;
                                    if(find(atom_names.begin(), atom_names.end(), "*") == atom_names.end())
                                    {
                                        AtomVector atoms = residue->GetAtoms();
                                        for(vector<string>::iterator it3 = atom_names.begin(); it3 != atom_names.end(); it3++)
                                        {
                                            string atom_name = *it3;
                                            int atom_name_search_type = 0;
                                            if(atom_name.at(0) == '^')
                                            {
                                                atom_name = atom_name.substr(1);
                                                atom_name_search_type = -1;
                                            }
                                            else if(atom_name.at(atom_name.size() - 1) == '$')
                                            {
                                                atom_name = atom_name.substr(0, atom_name.size() - 2);
                                                atom_name_search_type = 1;
                                            }
                                            else if(atom_name.at(0) == '#')
                                            {
                                                atom_name = atom_name.substr(1);
                                                atom_name_search_type = -2;
                                            }
                                            else
                                                atom_name_search_type = 0;
                                            for(AtomVector::iterator it4 = atoms.begin(); it4 != atoms.end(); it4++)
                                            {
                                                Atom* atom = *it4;
                                                switch(atom_name_search_type)
                                                {
                                                    case 0:
                                                    {
                                                        if(atom->GetName().compare(atom_name) == 0)
                                                            selection.push_back(atom);
                                                        break;
                                                    }
                                                    case 1:
                                                    {
                                                        if(atom->GetName().find(atom_name) != string::npos &&
                                                                atom->GetName().find(atom_name) == atom->GetName().size() - atom_name.size())
                                                            selection.push_back(atom);
                                                        break;
                                                    }
                                                    case -1:
                                                    {
                                                        if(atom->GetName().find(atom_name) != string::npos &&
                                                                atom->GetName().find(atom_name) == 0)
                                                            selection.push_back(atom);
                                                        break;
                                                    }
                                                    case -2:
                                                    {
                                                        string atom_serial_number = Split(atom->GetId(), "_").at(1);
                                                        int range_selection = 0;
                                                        if(atom_name.find("-") != string::npos)
                                                        {
                                                            range_selection = 1;

                                                        }
                                                        else
                                                            range_selection = 0;
                                                        switch(range_selection)
                                                        {
                                                            case 0:
                                                            {
                                                                if(atom_serial_number.find(atom_name) != string::npos &&
                                                                        atom_serial_number.find(atom_name) == 0)
                                                                    selection.push_back(atom);
                                                                break;
                                                            }
                                                            case 1:
                                                            {
                                                                string start_index = Split(atom_name, "-").at(0);
                                                                string end_index = Split(atom_name, "-").at(1);
                                                                int atom_serial_number_int = ConvertString<int>(atom_serial_number);
                                                                int start_index_int = ConvertString<int>(start_index);
                                                                int end_index_int = ConvertString<int>(end_index);
                                                                if(atom_serial_number_int >= start_index_int && atom_serial_number_int <= end_index_int)
                                                                    selection.push_back(atom);
                                                                break;
                                                            }
                                                        }
                                                        break;
                                                    }
                                                }

                                            }
                                        }
                                    }
                                    else
                                    {
                                        AtomVector atoms = residue->GetAtoms();
                                        for(AtomVector::iterator it4 = atoms.begin(); it4 != atoms.end(); it4++)
                                        {
                                            Atom* atom = *it4;
                                            selection.push_back(atom);
                                        }
                                    }
                                }
                            }

                        }

                    }
                }

            }
        }
    }

    return selection;
}

Assembly::SelectPatternMap Assembly::ParsePatternString(string pattern)
{
    SelectPatternMap select_pattern_map = SelectPatternMap();
    vector<string> tokens = Split(pattern, ";");
    for(int i = 0; i < tokens.size(); i++)
    {
        string token = tokens.at(i);
        vector<string> assembly_tokens = Split(token, ":");
        string assembly_id_token = assembly_tokens.at(0);
        vector<string> assembly_ids = Split(assembly_id_token, ",");
        map<string, vector<string> > residues = map<string, vector<string> >();
        for(int j = 1; j < assembly_tokens.size(); j++)
        {
            string assembly_token = assembly_tokens.at(j);
            vector<string> assembly_residue_tokens = Split(assembly_token, "@");
            string residue_token = assembly_residue_tokens.at(0);
            string atom_token = assembly_residue_tokens.at(1);
            vector<string> residue_tokens = Split(residue_token, ",");
            vector<string> atom_tokens = Split(atom_token, ",");
            for(int k = 0; k < residue_tokens.size(); k++)
            {
                for(int l = 0; l < atom_tokens.size(); l++)
                {
                    if(find(residues[residue_tokens.at(k)].begin(), residues[residue_tokens.at(k)].end(), atom_tokens.at(l)) == residues[residue_tokens.at(k)].end())
                        residues[residue_tokens.at(k)].push_back(atom_tokens.at(l));
                }
            }
        }
        for(int j = 0; j < assembly_ids.size(); j++)
            select_pattern_map[assembly_ids.at(j)] = residues;
    }
    return select_pattern_map;
}

void Assembly::GetHierarchicalMapOfAssembly(HierarchicalContainmentMap &hierarchical_map, stringstream &index)
{
    hierarchical_map[index.str()] = this->residues_;
    if(this->assemblies_.size() == 0)
        return;
    for(AssemblyVector::iterator it = this->assemblies_.begin(); it != this->assemblies_.end(); it++)
    {
        stringstream i;
        i << index.str() << "." << (*it)->GetSequenceNumber();
        Assembly* assembly = (*it);
        assembly->GetHierarchicalMapOfAssembly(hierarchical_map, i);
    }
}

void Assembly::ClearAssembly()
{
    this->residues_.clear();
    this->assemblies_.clear();
    //    this->source_file_ = "";
    //    this->source_file_type_ = UNKNOWN;
    //    this->chemical_type_ = "";
    //    this->description_ = "";
    //    this->sequence_number_ = 1;
    //    this->total_mass_ = dNotSet;
    //    this->center_of_geometry_ = Coordinate();
    //    this->center_of_mass_ = Coordinate();
    //    this->model_index_ = 0;
}

LibraryFileSpace::LibraryFile::ResidueMap Assembly::GetAllResiduesFromMultipleLibFilesMap(vector<string> lib_files)
{
    LibraryFileSpace::LibraryFile::ResidueMap all_residues;
    LibraryFileSpace::LibraryFile::ResidueMap residues;
    for(vector<string>::iterator it = lib_files.begin(); it != lib_files.end(); it++)
    {
        LibraryFileSpace::LibraryFile* lib_file = new LibraryFileSpace::LibraryFile(*it);
        residues = lib_file->GetResidues();
        for(LibraryFileSpace::LibraryFile::ResidueMap::iterator it1 = residues.begin(); it1 != residues.end(); it1++)
        {
            string lib_residue_name = (*it1).first;
            all_residues[lib_residue_name] = (*it1).second;
        }
    }
    return all_residues;
}

PrepFileSpace::PrepFile::ResidueMap Assembly::GetAllResiduesFromMultiplePrepFilesMap(vector<string> prep_files)
{
    PrepFileSpace::PrepFile::ResidueMap all_residues;
    PrepFileSpace::PrepFile::ResidueMap residues;
    for(vector<string>::iterator it = prep_files.begin(); it != prep_files.end(); it++)
    {
        PrepFileSpace::PrepFile* prep_file = new PrepFileSpace::PrepFile(*it);
        residues = prep_file->GetResidues();
        for(PrepFileSpace::PrepFile::ResidueMap::iterator it1 = residues.begin(); it1 != residues.end(); it1++)
        {
            string prep_residue_name = (*it1).first;
            all_residues[prep_residue_name] = (*it1).second;
        }
    }
    return all_residues;
}

ResidueNameMap Assembly::GetAllResidueNamesFromMultipleLibFilesMap(vector<string> lib_files)
{
    ResidueNameMap all_residue_names;
    vector<string> residue_names;
    for(vector<string>::iterator it = lib_files.begin(); it != lib_files.end(); it++)
    {
        LibraryFileSpace::LibraryFile* lib_file = new LibraryFileSpace::LibraryFile(*it);
        residue_names = lib_file->GetAllResidueNames();
        for(vector<string>::iterator it1 = residue_names.begin(); it1 != residue_names.end(); it1++)
        {
            string residue_name = (*it1);
            all_residue_names[residue_name] = (residue_name);
        }
    }
    return all_residue_names;
}

GlycamResidueNamingMap Assembly::ExtractResidueGlycamNamingMap(vector<Oligosaccharide*> oligosaccharides)
{
    //TODO: Done
    //Update the map
    //Key: string -> residue_id or atom_id
    //Value: vector<string> -> list of possible glycam naming
    GlycamResidueNamingMap pdb_glycam_residue_map = GlycamResidueNamingMap();
    //Iterates on all oligosaccharides and update the naming map
    for(vector<Oligosaccharide*>::iterator it = oligosaccharides.begin(); it != oligosaccharides.end(); it++)
    {
        int index = 0;
        Oligosaccharide* oligo = *it;
        string oligo_name = oligo->oligosaccharide_name_;
        //In case that there is no terminal attached to the reducing end adds a temporary terminal residue to make the sequence parser able to parse the sequence
        if(oligo->terminal_.compare("") == 0)
            oligo_name = oligo_name + "1-OH";
        CondensedSequence* condensed_sequence = new CondensedSequence(oligo_name);
        //Gets the three letter code of all carbohydrates involved in current oligosaccharide
        CondensedSequence::CondensedSequenceAmberPrepResidueTree condensed_sequence_amber_residue_tree = condensed_sequence->GetCondensedSequenceAmberPrepResidueTree();
        if(oligo->terminal_.compare("") != 0)
        {
            Atom* anomeric_o = NULL;
            Atom* anomeric_c = NULL;
            if(oligo->root_->cycle_atoms_.at(0) != NULL)
                anomeric_c = oligo->root_->cycle_atoms_.at(0);
            if(oligo->root_->side_atoms_.at(0).at(1) != NULL)
                anomeric_o = oligo->root_->side_atoms_.at(0).at(1);
            if(anomeric_o != NULL && anomeric_c != NULL)
            {
                AtomVector terminal_atoms = AtomVector();
                string terminal = CheckTerminals(anomeric_o, terminal_atoms);
                //Separating terminal residue in glycam naming
                if(terminal.compare("") != 0)
                {
                    //TODO:
                    //Add the residue mismatch into a structure for the Ontology usage
                    for(AtomVector::iterator it1 = terminal_atoms.begin(); it1 != terminal_atoms.end(); it1++)
                    {
                        Atom* terminal_atom = *it1;
                        string terminal_atom_id = terminal_atom->GetId();
                        string terminal_residue_id = terminal_atom->GetResidue()->GetId();
                        if(pdb_glycam_residue_map.find(terminal_atom_id) == pdb_glycam_residue_map.end())
                            pdb_glycam_residue_map[terminal_atom_id] == vector<string>();
                        pdb_glycam_residue_map[terminal_atom_id].push_back(condensed_sequence_amber_residue_tree.at(index)->GetName());
                        if(pdb_glycam_residue_map.find(terminal_residue_id) == pdb_glycam_residue_map.end())
                            pdb_glycam_residue_map[terminal_residue_id] = vector<string>();
                        pdb_glycam_residue_map[terminal_residue_id].push_back(condensed_sequence_amber_residue_tree.at(index)->GetName());
                    }
                }
            }
        }

        index++;
        this->ExtractOligosaccharideNamingMap(pdb_glycam_residue_map, oligo, condensed_sequence_amber_residue_tree, index);
    }

    return pdb_glycam_residue_map;
}

void Assembly::ExtractOligosaccharideNamingMap(GlycamResidueNamingMap& pdb_glycam_map, Oligosaccharide *oligosaccharide,
                                               CondensedSequence::CondensedSequenceAmberPrepResidueTree condensed_sequence_amber_residue_tree, int &index)
{
    string name = condensed_sequence_amber_residue_tree.at(index)->GetName();
    //TODO: Done
    //Update to hold all possible three letter names for the specific residue_id
    string residue_id = oligosaccharide->root_->cycle_atoms_.at(0)->GetResidue()->GetId();
    if(pdb_glycam_map.find(residue_id) == pdb_glycam_map.end())
        pdb_glycam_map[residue_id] = vector<string>();
    pdb_glycam_map[residue_id].push_back(name);
    index++;
    //Separating SO3 and PO3 residues in glycam naming
    while(index < condensed_sequence_amber_residue_tree.size() && condensed_sequence_amber_residue_tree.at(index)->GetIsDerivative())
    {
        int parent_index = condensed_sequence_amber_residue_tree.at(index)->GetParentId();
        int carbon_index = ConvertString<int>(condensed_sequence_amber_residue_tree.at(index)->GetAnomericCarbon().substr(1));
        Atom* carbon_atom = NULL;
        if(oligosaccharide->root_->derivatives_map_.find("-1") == oligosaccharide->root_->derivatives_map_.end())
        {
            if((carbon_index == 6 || carbon_index == 7 || carbon_index == 8) && carbon_atom == NULL)
                if(oligosaccharide->root_->sugar_name_.ring_type_.compare("P") == 0)
                    carbon_atom = oligosaccharide->root_->side_atoms_.at(oligosaccharide->root_->side_atoms_.size()-1).at(carbon_index-6);
            if((carbon_index == 5 || carbon_index == 6 || carbon_index == 7) && carbon_atom == NULL)
                if(oligosaccharide->root_->sugar_name_.ring_type_.compare("F") == 0)
                    carbon_atom = oligosaccharide->root_->side_atoms_.at(oligosaccharide->root_->side_atoms_.size()-1).at(carbon_index-5);
            if(carbon_atom == NULL)
                carbon_atom = oligosaccharide->root_->side_atoms_.at(carbon_index).at(0);
        }
        else
        {
            if((carbon_index == 7 || carbon_index == 8 || carbon_index == 9) && carbon_atom == NULL)
                if(oligosaccharide->root_->sugar_name_.ring_type_.compare("P") == 0)
                    carbon_atom = oligosaccharide->root_->side_atoms_.at(oligosaccharide->root_->side_atoms_.size()-1).at(carbon_index-7);
            if((carbon_index == 6 || carbon_index == 7 || carbon_index == 8) && carbon_atom == NULL)
                if(oligosaccharide->root_->sugar_name_.ring_type_.compare("F") == 0)
                    carbon_atom = oligosaccharide->root_->side_atoms_.at(oligosaccharide->root_->side_atoms_.size()-1).at(carbon_index-6);
            if(carbon_atom == NULL)
                carbon_atom = oligosaccharide->root_->side_atoms_.at(carbon_index-1).at(0);
        }
        if(carbon_atom != NULL)
        {
            AtomVector n_linkage_derivative_atoms = AtomVector();
            string n_linkage_derivative_string = CheckxC_NxO_SO3(carbon_atom, oligosaccharide->root_->cycle_atoms_str_, 'N', n_linkage_derivative_atoms);
            if(n_linkage_derivative_string.compare("xC-N-SO3") == 0)
            {
                for(AtomVector::iterator it = n_linkage_derivative_atoms.begin(); it != n_linkage_derivative_atoms.end(); it++)
                {
                    Atom* atom = *it;
                    string derivative_atom_id = atom->GetId();
                    if(pdb_glycam_map.find(derivative_atom_id) == pdb_glycam_map.end())
                        pdb_glycam_map[derivative_atom_id] = vector<string>();
                    pdb_glycam_map[derivative_atom_id].push_back("SO3");
                }
                string new_name = condensed_sequence_amber_residue_tree.at(parent_index)->GetName();
                string derivative_residue_id = oligosaccharide->root_->cycle_atoms_.at(0)->GetResidue()->GetId();
                if(pdb_glycam_map.find(derivative_residue_id) == pdb_glycam_map.end())
                    pdb_glycam_map[derivative_residue_id] = vector<string>();
                pdb_glycam_map[derivative_residue_id].push_back(ConvertT<int>(carbon_index) + new_name.substr(1));
            }
            AtomVector o_linkage_derivative_atoms = AtomVector();
            string o_linkage_derivative_string = CheckxC_NxO_SO3(carbon_atom, oligosaccharide->root_->cycle_atoms_str_, 'O', o_linkage_derivative_atoms);
            if(o_linkage_derivative_string.compare("xC-O-SO3") == 0)
            {
                for(AtomVector::iterator it = o_linkage_derivative_atoms.begin(); it != o_linkage_derivative_atoms.end(); it++)
                {
                    Atom* atom = *it;
                    string derivative_atom_id = atom->GetId();
                    if(pdb_glycam_map.find(derivative_atom_id) == pdb_glycam_map.end())
                        pdb_glycam_map[derivative_atom_id] = vector<string>();
                    pdb_glycam_map[derivative_atom_id].push_back("SO3");
                }
                string new_name = condensed_sequence_amber_residue_tree.at(parent_index)->GetName();
                string derivative_residue_id = oligosaccharide->root_->cycle_atoms_.at(0)->GetResidue()->GetId();
                if(pdb_glycam_map.find(derivative_residue_id) == pdb_glycam_map.end())
                    pdb_glycam_map[derivative_residue_id] = vector<string>();
                pdb_glycam_map[derivative_residue_id].push_back(ConvertT<int>(carbon_index) + new_name.substr(1));
            }
            AtomVector n_linkage_derivative_atoms_1 = AtomVector();
            string n_linkage_derivative_string_1 = CheckxC_NxO_PO3(carbon_atom, oligosaccharide->root_->cycle_atoms_str_, 'N', n_linkage_derivative_atoms_1);
            if(n_linkage_derivative_string_1.compare("xC-N-PO3") == 0)
            {
                for(AtomVector::iterator it = n_linkage_derivative_atoms_1.begin(); it != n_linkage_derivative_atoms_1.end(); it++)
                {
                    Atom* atom = *it;
                    string derivative_atom_id = atom->GetId();
                    if(pdb_glycam_map.find(derivative_atom_id) == pdb_glycam_map.end())
                        pdb_glycam_map[derivative_atom_id] = vector<string>();
                    pdb_glycam_map[derivative_atom_id].push_back("PO3");
                }
                string new_name = condensed_sequence_amber_residue_tree.at(parent_index)->GetName();
                string derivative_residue_id = oligosaccharide->root_->cycle_atoms_.at(0)->GetResidue()->GetId();
                if(pdb_glycam_map.find(derivative_residue_id) == pdb_glycam_map.end())
                    pdb_glycam_map[derivative_residue_id] = vector<string>();
                pdb_glycam_map[derivative_residue_id].push_back(ConvertT<int>(carbon_index) + new_name.substr(1));
            }
            AtomVector o_linkage_derivative_atoms_1 = AtomVector();
            string o_linkage_derivative_string_1 = CheckxC_NxO_PO3(carbon_atom, oligosaccharide->root_->cycle_atoms_str_, 'O', o_linkage_derivative_atoms_1);
            if(o_linkage_derivative_string_1.compare("xC-O-PO3") == 0)
            {
                for(AtomVector::iterator it = o_linkage_derivative_atoms_1.begin(); it != o_linkage_derivative_atoms_1.end(); it++)
                {
                    Atom* atom = *it;
                    string derivative_atom_id = atom->GetId();
                    if(pdb_glycam_map.find(derivative_atom_id) == pdb_glycam_map.end())
                        pdb_glycam_map[derivative_atom_id] = vector<string>();
                    pdb_glycam_map[derivative_atom_id].push_back("PO3");
                }
                string new_name = condensed_sequence_amber_residue_tree.at(parent_index)->GetName();
                string derivative_residue_id = oligosaccharide->root_->cycle_atoms_.at(0)->GetResidue()->GetId();
                if(pdb_glycam_map.find(derivative_residue_id) == pdb_glycam_map.end())
                    pdb_glycam_map[derivative_residue_id] = vector<string>();
                pdb_glycam_map[derivative_residue_id].push_back(ConvertT<int>(carbon_index) + new_name.substr(1));
            }
        }
        index++;
    }

    //Recursively assign glycam naming to the monosaccharides of an oligosaccharide
    for(unsigned int i = 0; i < oligosaccharide->child_oligos_.size(); i++)
    {
        this->ExtractOligosaccharideNamingMap(pdb_glycam_map, oligosaccharide->child_oligos_.at(i), condensed_sequence_amber_residue_tree, index);
    }
}

void Assembly::UpdateResidueName2GlycamName(GlycamResidueNamingMap residue_glycam_map, string prep_file)
{
    for(AssemblyVector::iterator it = this->GetAssemblies().begin(); it != this->GetAssemblies().end(); it++)
        (*it)->UpdateResidueName2GlycamName(residue_glycam_map, prep_file);

    PrepFile* prep = new PrepFile(prep_file);
    PrepFile::ResidueMap prep_residues = prep->GetResidues();
    ResidueVector residues = this->GetResidues();
    ResidueVector updated_residues = ResidueVector();
    int residue_sequence_number = 0;
    ResidueVector terminal_residues = ResidueVector();
    for(ResidueVector::iterator it2 = residues.begin(); it2 != residues.end(); it2++)
    {
        Residue* residue = *it2;
        string residue_name = residue->GetName();
        string residue_id = residue->GetId();
        if(residue_glycam_map.find(residue_id) != residue_glycam_map.end())
        {
            terminal_residues.push_back(new Residue(this, ""));
            residue_sequence_number--;
            bool terminal = false;
            int residue_name_size = residue_name.size();
            //TODO: Done
            //Update to match the atoms of the residue with all the three-letter code prep residues to set the updated residue name correspondingly
            //Create a map between the current atom names and the corresponding prep atom names when the match is found
            //Add the residue mismatch into a structure for the Ontology usage
            vector<string> glycam_names = residue_glycam_map[residue_id];
            GlycamAtomNameMap pdb_glycam_map = GlycamAtomNameMap();
            GlycamAtomNameMap glycam_pdb_map = GlycamAtomNameMap();
            string glycam_name = *glycam_names.begin();
            ResidueVector query_residues = ResidueVector();
            for(vector<string>::iterator name_it = glycam_names.begin(); name_it != glycam_names.end(); name_it++)
            {
                string glycam_residue_name = *name_it;
                PrepFile::ResidueMap customized_prep_residues = PrepFile::ResidueMap();
                customized_prep_residues[glycam_residue_name] = prep_residues[glycam_residue_name];
                PrepFile* temp_prep_file = new PrepFile();
                temp_prep_file->SetResidues(customized_prep_residues);
                temp_prep_file->Write(glycam_residue_name + ".prep");

                Assembly* prep_assembly = new Assembly();
                prep_assembly->BuildAssemblyFromPrepFile(temp_prep_file);
                prep_assembly->SetSourceFile(glycam_residue_name + ".prep");
                prep_assembly->BuildStructureByPrepFileInformation();
                remove((glycam_residue_name + ".prep").c_str());

                query_residues.push_back(prep_assembly->GetResidues().at(0));
            }
                //TODO:
                //The two residues atom names are matched but the geometry needs to be checked
                //if GeometryCheck returns true
//            PatternMatching(residue, query_residues, pdb_glycam_map, glycam_pdb_map);

            AtomVector atoms = residue->GetAtoms();
            map<string, Residue*> residue_set = map<string, Residue*>();
            for(AtomVector::iterator it1 = atoms.begin(); it1 != atoms.end(); it1++)
            {
                Atom* atom = *it1;
                string atom_id = atom->GetId();
                //Terminal residue naming
                if(residue_glycam_map.find(atom_id) != residue_glycam_map.end())
                {
                    terminal_residues.at(terminal_residues.size() - 1)->SetAssembly(this);
                    glycam_name = *residue_glycam_map[atom_id].begin();
                    terminal_residues.at(terminal_residues.size() - 1)->SetName(glycam_name);
                    terminal_residues.at(terminal_residues.size() - 1)->AddAtom(atom);
                    vector<string> residue_id_tokens = Split(residue_id, "_");
                    string terminal_residue_id = residue_id_tokens.at(0) + "_" + residue_id_tokens.at(1) + "_" + ConvertT<int>(residue_sequence_number) + "_"
                            + residue_id_tokens.at(3) + "_" + residue_id_tokens.at(4);
                    terminal_residues.at(terminal_residues.size() - 1)->SetId(terminal_residue_id);
                    atom->SetResidue(terminal_residues.at(terminal_residues.size() - 1));
                    if(!terminal)
                    {
                        updated_residues.push_back(terminal_residues.at(terminal_residues.size() - 1));
                        terminal = true;
                    }
                    int index = atom_id.find(residue_name);
                    if(index >= 0)
                        atom->SetId(atom_id.replace(index, residue_name_size, glycam_name));
                    vector<string> atom_id_tokens = Split(atom->GetId(), "_");
                    string atom_id_new = atom_id_tokens.at(0) + "_" + atom_id_tokens.at(1) + "_" + atom_id_tokens.at(2) + "_" + atom_id_tokens.at(3) + "_" +
                            ConvertT<int>(residue_sequence_number) + "_" + atom_id_tokens.at(5) + "_" + atom_id_tokens.at(6);
                    atom->SetId(atom_id_new);
                }
                //Non-terminal residues glycam naming
                else
                {
                    //TODO:
                    //Update to match the atoms with the corresponding prep residue to change the atom names with respect to prep file
                    //Use the map to update the atom naming
                    //Add the atom name mismatch into a structure for the Ontology usage
                    string prep_atom_id;
                    if(pdb_glycam_map.find(atom_id) != pdb_glycam_map.end())
                        prep_atom_id = pdb_glycam_map[atom_id];
                    else
                        prep_atom_id = atom_id;
                    glycam_name = Split(prep_atom_id,"_")[2];
                    if(residue_set.find(glycam_name) == residue_set.end())
                    {
                        residue_set[glycam_name] = (new Residue(this, glycam_name));
                        residue_sequence_number--;
                        vector<string> res_id_tokens = Split(residue_id, "_");
                        string res_id = glycam_name + "_" + res_id_tokens.at(1) + "_" + ConvertT<int>(residue_sequence_number) + "_"
                                + res_id_tokens.at(3) + "_" + res_id_tokens.at(4);
                        residue_set[glycam_name]->SetId(res_id);
                    }
                    string atom_name = atom->GetName();
                    string new_atom_name = Split(prep_atom_id,"_")[0];
                    string new_atom_id = atom_id;
                    int index = new_atom_id.find(residue_name);
                    if(index >= 0)
                        new_atom_id = new_atom_id.replace(index, residue_name_size, glycam_name);
                    index = new_atom_id.find(atom_name);
                    if(index >= 0)
                        new_atom_id = new_atom_id.replace(index, atom_name.size(), new_atom_name);
                    atom->SetId(new_atom_id);
                    atom->SetResidue(residue_set[glycam_name]);
                    residue_set[glycam_name]->AddAtom(atom);
                }
            }
            for(map<string, Residue*>::iterator it1 = residue_set.begin(); it1 != residue_set.end(); it1++)
            {
                updated_residues.push_back((*it1).second);
            }
        }
        else
            updated_residues.push_back(residue);
    }
    this->SetResidues(updated_residues);
}

void Assembly::DetectShape(AtomVector cycle, Monosaccharide* mono)
{
    ///Creating a new assembly only from the ring atoms for external detect shape program
    Assembly* detect_shape_assembly = new Assembly();
    detect_shape_assembly->AddResidue(cycle.at(0)->GetResidue());
    Residue* detect_shape_residue = detect_shape_assembly->GetResidues().at(0);
    for(int i = 0; i < cycle.size(); i++)
    {
        string name = cycle.at(i)->GetName();
        string id = cycle.at(i)->GetId();

        ///Preparing atom name and id for detect shape program. It doesn't work with atoms containing special characters: C', C* etc
        //        replace( id.begin(), id.end(), '?', 'n'); // replace all '?' with 'n'
        FindReplaceString(id, "\'", "");
        FindReplaceString(id, ",", "");
        FindReplaceString(name, "*", "");
        replace( id.begin(), id.end(), '*', 's'); // replace all '*' with ''

        //        replace( name.begin(), name.end(), '?', 'n'); // replace all '?' with 'n'
        FindReplaceString(name, "\'", "");
        FindReplaceString(name, ",", "");
        FindReplaceString(name, "*", "");

        cycle.at(i)->SetName(name);
        cycle.at(i)->SetId(id);
    }
    detect_shape_residue->SetAtoms(cycle);

    ///Write a new PDB file from the new assembly
    PdbFile* pdb = detect_shape_assembly->BuildPdbFileStructureFromAssembly();
    pdb->Write("temp_gmml_pdb.pdb");

    ///Converting the written PDB file to fomrat readable by detect_shape program
    string line = "";
    ifstream gmml_pdb ("temp_gmml_pdb.pdb");
    ofstream detect_shape_pdb ("temp_detect_shape_pdb.pdb");
    int n = 0;
    if (gmml_pdb.is_open())
    {
        while (!gmml_pdb.eof()) {
            getline(gmml_pdb, line);
            if(line.find("HETATM") != string::npos)
            {
                detect_shape_pdb << line << endl;
            }
            n++;
        }
        gmml_pdb.close();
        detect_shape_pdb.close();
    }
    else cout << "Unable to open temp_gmml_pdb.pdb file" << endl;

    ///Writing a configuration file for the second argument of the detect_sugar program
    ofstream detect_shape_configuration ("temp_config");
    detect_shape_configuration << "Atom" << endl;
    for(int i = 0; i < cycle.size(); i++)
    {
        detect_shape_configuration << cycle.at(i)->GetName() << endl;
    }
    detect_shape_configuration << "Residue" << endl;
    detect_shape_configuration << "1" << endl;
    detect_shape_configuration << "Path" << endl;
    detect_shape_configuration << "apps/BFMP/canonicals.txt" << endl;
    detect_shape_configuration.close();

    ///Calling detect_shape program
    system("apps/BFMP/detect_shape temp_detect_shape_pdb.pdb temp_config > /dev/null");

    ///Adding the BFMP ring conformation infomration gained from the detect_sugar program to the monosaccharide
    ifstream shape_detection_result ("ring_conformations.txt");
    line = "";
    if (shape_detection_result.is_open())
    {
        getline (shape_detection_result,line);
        getline (shape_detection_result,line);
        vector<string> line_tokens = Split(line, "\t");
        if(line_tokens.at(1).compare("-") == 0)
            mono->bfmp_ring_conformation_ = line_tokens.at(2);
        else
            mono->bfmp_ring_conformation_ = line_tokens.at(1);
        shape_detection_result.close();
    }
    else cout << "Unable to open ring_conformations.txt file from detect shape program" << endl;

    ///Deleting temporary files
    remove("temp_detect_shape_pdb.pdb");
    remove("temp_gmml_pdb.pdb");
    remove("temp_config");
    remove("ring_conformations.txt");
}

bool Assembly::PatternMatching(Residue *residue, ResidueVector query_residues, GlycamAtomNameMap &pdb_glycam_map, GlycamAtomNameMap& glycam_atom_map)
{
    CreatePrunedMatchingGraph(residue, query_residues);
    for(ResidueVector::iterator it = query_residues.begin(); it != query_residues.end(); it++)
    {
        Residue* query_residue = *it;
        AtomVector lowest_degree_atoms = query_residue->GetAtomsWithLowestIntraDegree();
        stack<BacktrackingElements*> backtracking_stack_1 = stack<BacktrackingElements*>();
        for(unsigned int i = lowest_degree_atoms.size() - 1; i >= 0; i--)
            backtracking_stack_1.push(new BacktrackingElements(pdb_glycam_map, glycam_atom_map, lowest_degree_atoms, i));
        if(!backtracking_stack_1.empty())
        {
            // backtracking point
            BacktrackingElements* backtracking_point_1 = backtracking_stack_1.top();
            backtracking_stack_1.pop();
            pdb_glycam_map = backtracking_point_1->pdb_glycam_map_;
            glycam_atom_map = backtracking_point_1->glycam_atom_map_;
            lowest_degree_atoms = backtracking_point_1->atoms_;
            Atom* source_atom = lowest_degree_atoms.at(backtracking_point_1->index_);
            AtomVector source_atom_intra_neighbors = source_atom->GetNode()->GetIntraNodeNeighbors();
            stack<BacktrackingElements*> backtracking_stack_2 = stack<BacktrackingElements*>();
            for(unsigned int i = source_atom_intra_neighbors.size() - 1; i >= 0; i--)
                backtracking_stack_2.push(new BacktrackingElements(pdb_glycam_map, glycam_atom_map, source_atom_intra_neighbors, i));
            if(!backtracking_stack_2.empty())
            {
                // backtracking point
                BacktrackingElements* backtracking_point_2 = backtracking_stack_2.top();
                backtracking_stack_2.pop();
                pdb_glycam_map = backtracking_point_2->pdb_glycam_map_;
                glycam_atom_map = backtracking_point_2->glycam_atom_map_;
                source_atom_intra_neighbors = backtracking_point_2->atoms_;
                Atom* mapped_atom = source_atom_intra_neighbors.at(backtracking_point_2->index_);
                pdb_glycam_map[mapped_atom->GetId()] = source_atom->GetId();
                glycam_atom_map[source_atom->GetId()] = mapped_atom->GetId();
                AtomVector source_atom_neighbors = source_atom->GetNode()->GetNodeNeighbors();
                queue<Atom*> to_visit = queue<Atom*>();
                queue<Atom*> parents = queue<Atom*>();
                for(unsigned int i = 0; i < source_atom_neighbors.size(); i++)
                {
                    if(glycam_atom_map.find(source_atom_neighbors.at(i)->GetId()) == glycam_atom_map.end())
                    {
                        to_visit.push(source_atom_neighbors.at(i));
                        parents.push(mapped_atom);
                    }
                }
                stack<BacktrackingElements*> backtracking_stack_3 = stack<BacktrackingElements*>();
                while(!to_visit.empty())
                {
                    Atom* atom = to_visit.front();
                    Atom* parent = parents.front();
                    to_visit.pop();
                    parents.pop();
                    AtomVector atom_intra_neighbors = atom->GetNode()->GetIntraNodeNeighbors();
                    AtomVector atom_intra_neighbors_filtered_with_parent = AtomVector();
                    for(unsigned int i = 0; i < atom_intra_neighbors.size(); i++)
                    {
                        AtomVector neighbors = atom_intra_neighbors.at(i)->GetNode()->GetNodeNeighbors();
                        if(find(neighbors.begin(), neighbors.end(), parent) != neighbors.end() &&
                                pdb_glycam_map.find(atom_intra_neighbors.at(i)->GetId()) == pdb_glycam_map.end())
                            atom_intra_neighbors_filtered_with_parent.push_back(atom_intra_neighbors.at(i));
                    }
                    for(unsigned int i = atom_intra_neighbors_filtered_with_parent.size() - 1; i >= 0; i--)
                        backtracking_stack_3.push(new BacktrackingElements(pdb_glycam_map, glycam_atom_map,
                                                                           atom_intra_neighbors_filtered_with_parent, i, to_visit));
                    if(!backtracking_stack_3.empty())
                    {
                        // backtracking point
                        BacktrackingElements* backtracking_point_3 = backtracking_stack_3.top();
                        backtracking_stack_3.pop();
                        pdb_glycam_map = backtracking_point_3->pdb_glycam_map_;
                        glycam_atom_map = backtracking_point_3->glycam_atom_map_;
                        atom_intra_neighbors_filtered_with_parent = backtracking_point_3->atoms_;
                        Atom* mapped = atom_intra_neighbors_filtered_with_parent.at(backtracking_point_3->index_);
                        pdb_glycam_map[mapped->GetId()] = atom->GetId();
                        glycam_atom_map[atom->GetId()] = mapped->GetId();
                        AtomVector atom_neighbors = atom->GetNode()->GetNodeNeighbors();
                        for(unsigned int i = 0; i < atom_neighbors.size(); i++)
                        {
                            if(glycam_atom_map.find(atom_neighbors.at(i)->GetId()) == glycam_atom_map.end())
                            {
                                to_visit.push(atom_neighbors.at(i));
                                parents.push(mapped);
                            }
                        }
                    }
                }
                // Check if the mapping is complete
                // continue to the next residue
                // else
                // backtrack
            }
            // Check if the mapping is complete
            // continue to the next residue
            // else
            // backtrack
        }
        // Check if the mapping is complete
        // continue to the next residue
        // else
        // backtrack
    }
    // Check if the mapping is complete
    // return true
    // else
    // return false
}

void Assembly::CreateLabelGraph(Residue *residue, Residue *query_residue)
{
    AtomVector query_atoms = query_residue->GetAtoms();
    AtomVector atoms = residue->GetAtoms();
    for(AtomVector::iterator it = query_atoms.begin(); it != query_atoms.end(); it++)
    {
        Atom* query_atom = *it;
        AtomNode* query_atom_node = query_atom->GetNode();
        AtomVector query_atom_intra_neighbors = AtomVector();
        if(query_atom_node == NULL)
            query_atom_node = new AtomNode();
        for(AtomVector::iterator it1 = atoms.begin(); it1 != atoms.end(); it1++)
        {
            Atom* atom = *it1;
            AtomNode* atom_node = atom->GetNode();
            if(atom_node == NULL)
                atom_node = new AtomNode();
            if(query_atom_node->GetElementLabel().compare(atom_node->GetElementLabel()) == 0)
                query_atom_intra_neighbors.push_back(atom);
        }
        query_atom_node->SetIntraNodeNeighbors(query_atom_intra_neighbors);
    }
}

void Assembly::PruneLabelGraphByNeighboringLabels(Residue *query_residue)
{
    AtomVector query_atoms = query_residue->GetAtoms();
    for(AtomVector::iterator it = query_atoms.begin(); it != query_atoms.end(); it++)
    {
        Atom* query_atom = *it;
        AtomNode* query_atom_node = query_atom->GetNode();
        if(query_atom_node != NULL)
        {
            string query_atom_neighborhood_label = query_atom_node->CreateNeighboringLabel();
            AtomVector query_atom_intra_neighbors = query_atom_node->GetIntraNodeNeighbors();
            AtomVector updated_query_atom_intra_neighbors = AtomVector();
            for(AtomVector::iterator it1 = query_atom_intra_neighbors.begin(); it1 != query_atom_intra_neighbors.end(); it1++)
            {
                Atom* query_atom_intra_neighbor = *it1;
                AtomNode* query_atom_intra_neighbor_node = query_atom_intra_neighbor->GetNode();
                if(query_atom_intra_neighbor_node != NULL)
                {
                    string query_atom_intra_neighbor_neighborhood_label = query_atom_intra_neighbor_node->CreateNeighboringLabel();
                    if(query_atom_neighborhood_label.compare(query_atom_intra_neighbor_neighborhood_label) == 0)
                        updated_query_atom_intra_neighbors.push_back(query_atom_intra_neighbor);
                }
            }
            query_atom_node->SetIntraNodeNeighbors(updated_query_atom_intra_neighbors);
        }
    }
}

void Assembly::CreatePrunedMatchingGraph(Residue *residue, ResidueVector query_residues)
{
    if(residue->GraphElementLabeling())
    {
        for(ResidueVector::iterator it = query_residues.begin(); it != query_residues.end(); it++)
        {
            Residue* query_residue = *it;
            if(query_residue->GraphElementLabeling())
            {
                this->CreateLabelGraph(residue, query_residue);
                this->PruneLabelGraphByNeighboringLabels(query_residue);
//                query_residue->Print();
            }
        }
    }
}

vector<Oligosaccharide*> Assembly::ExtractSugars(vector<string> amino_lib_files, bool glyprobity_report, bool populate_ontology)
{
    ResidueNameMap dataset_residue_names = GetAllResidueNamesFromMultipleLibFilesMap(amino_lib_files);

    ///CYCLE DETECTION
    CycleMap cycles = DetectCyclesByExhaustiveRingPerception();

    //    CycleMap cycles = DetectCyclesByDFS();

    ///PRINTING ALL DETECTED CYCLES
    cout << endl << "All detected cycles" << endl;
    for(CycleMap::iterator it = cycles.begin(); it != cycles.end(); it++)
    {
        string cycle_atoms_str = (*it).first;
        cout << cycle_atoms_str << endl;
    }

    ///FILTERING OUT FUSED CYCLES. aka Cycles that are sharing an edge
    RemoveFusedCycles(cycles);

    ///FILTERING OUT OXYGENLESS CYCLES
    FilterAllCarbonCycles(cycles);

    ///ANOMERIC CARBON DETECTION and SORTING
    cout << endl << "Cycles after discarding rings that are all-carbon" << endl;
    vector<string> anomeric_carbons_status = vector<string>();
    vector<Note*> anomeric_notes = vector<Note*>();
    CycleMap sorted_cycles = CycleMap();
    for(CycleMap::iterator it = cycles.begin(); it != cycles.end(); it++)
    {
        string cycle_atoms_str = (*it).first;
        AtomVector cycle_atoms = (*it).second;

        cout << cycle_atoms_str << endl; ///e.g. C1_3810_NAG_A_1521_?_?_1-O5_3821_NAG_A_1521_?_?_1-C5_3814_NAG_A_1521_?_?_1-C4_3813_NAG_A_1521_?_?_1-C3_3812_NAG_A_1521_?_?_1-C2_3811_NAG_A_1521_?_?_1

        Note* anomeric_note = new Note();
        Atom* anomeric = FindAnomericCarbon(anomeric_note, anomeric_carbons_status, cycle_atoms, cycle_atoms_str);
        anomeric_notes.push_back(anomeric_note);
        if(anomeric != NULL) ///Sorting the cycle atoms and adding it to the sorted_cycles map if an anomeric carbon identified, otherwise the structure can't be a sugar
        {
            AtomVector sorted_cycle_atoms = AtomVector();
            stringstream sorted_cycle_stream;
            sorted_cycle_atoms = SortCycle(cycle_atoms, anomeric, sorted_cycle_stream);
            sorted_cycles[sorted_cycle_stream.str()] = sorted_cycle_atoms;
        }
    }
    cycles = sorted_cycles;

    ///CREATING MONOSACCHARIDE STRUCTURE. Ring atoms, side atoms, chemical code (Glycode), modifications/derivatives, names
    cout << endl << "Detailed information of sorted cycles after discarding fused or oxygenless rings: " << endl;
    vector<Monosaccharide*> monos = vector<Monosaccharide*>();
    int mono_id = 0;
    for(CycleMap::iterator it = cycles.begin(); it != cycles.end(); it++)
    {
        string cycle_atoms_str = (*it).first;
        AtomVector cycle = (*it).second;

        Monosaccharide* mono = new Monosaccharide();
        int status_index = distance(cycles.begin(), it);
        mono->anomeric_status_ = anomeric_carbons_status.at(status_index);

        ///ASSIGNING RING ATOMS TO MONOSACCHARIDE OBJECT
        cout << "Ring atoms: " << cycle_atoms_str << endl;
        mono->cycle_atoms_str_ = cycle_atoms_str;
        mono->cycle_atoms_ = cycle;

        ///ASSIGNING SIDE ATOMS (EXCOCYCLIC ATOMS) TO MONOSACCHARIDE OBJECT
        vector<string> orientations = GetSideGroupOrientations(mono, cycle_atoms_str);

        ///PRINTING ASSIGNED SIDE ATOMS
        cout << "Side group atoms: " << endl;
        for(vector<AtomVector>::iterator it1 = mono->side_atoms_.begin(); it1 != mono->side_atoms_.end(); it1++)
        {
            AtomVector sides = (*it1);
            if(it1 == mono->side_atoms_.begin())///side atoms of anomeric carbon
            {
                if(sides.at(0) != NULL && sides.at(1) != NULL)
                    cout << "[1] -> " << sides.at(0)->GetId() << ", " << sides.at(1)->GetId() << endl;
                else if(sides.at(1) != NULL)
                    cout << "[1] -> " << sides.at(1)->GetId() << endl;
                else if(sides.at(0) != NULL)
                    cout << "[1] -> " << sides.at(0)->GetId() << endl;
            }
            else if(it1 == mono->side_atoms_.end() - 1)//side atoms of last carbon of the ring
            {
                cout << "[" << mono->cycle_atoms_.size() - 1 << "] -> ";
                if(sides.at(0) != NULL)
                    cout << sides.at(0)->GetId() << endl;
            }
            else if(sides.at(1) != NULL)
            {
                int cycle_atom_index = distance(mono->side_atoms_.begin(), it1);
                cout << "[" << cycle_atom_index + 1 << "] -> " << sides.at(1)->GetId() << endl;
            }
        }

        ///PRINTING ANOMERIC STATUS
        cout << mono->anomeric_status_ << mono->cycle_atoms_.at(0)->GetId() << endl;

        ///CREATING CHEMICAL CODE (Glycode) OBJECT
        ChemicalCode* code = BuildChemicalCode(orientations);
        if(code != NULL)
        {
            mono->chemical_code_ = code;
        }
        cout << endl << "Stereo chemistry chemical code:"  << endl;
        code->Print(cout);
        cout << endl;

        ///DETECT SHAPE USING BFMP EXTERNAL PROGRAM. Currently, the program does not work for furanoses
        if(cycle.size() > 5)
        {
            DetectShape(cycle, mono);
            if(mono->bfmp_ring_conformation_.compare("") != 0)
                cout << "BFMP ring conformation: " << mono->bfmp_ring_conformation_ << endl << endl; ///Part of Glyprobity report
        }

        ///CHECKING FOR +2 and +3 SIDE CARBONS
        AtomVector plus_sides = ExtractAdditionalSideAtoms(mono);

        ///FINDING CHEMICAL CODE IN NAME LOOKUP TABLE
        string code_str = code->toString();
        mono->sugar_name_ = SugarStereoChemistryNameLookup(code_str);

        ///DERIVATIVE/MODIFICATION PATTERN EXTRACTION
        ExtractDerivatives(mono);
        bool minus_one = false;
        if(mono->derivatives_map_.find("-1") != mono->derivatives_map_.end())
            minus_one = true;
        for(map<string, string>::iterator it1 = mono->derivatives_map_.begin(); it1 != mono->derivatives_map_.end(); it1++)
        {
            string key = (*it1).first;
            string value = (*it1).second;
            if(minus_one)
            {
                if(key.compare("-1") == 0)
                    cout << "Carbon at position 1 is attached to " << value << endl;
                else if(key.compare("a") == 0)
                    cout << "Carbon at position 2 is attached to " << value << endl;
                else if(key.compare("+1") == 0 || key.compare("+2") == 0 || key.compare("+3") == 0)
                    cout << "Carbon at position " << mono->cycle_atoms_.size() - 1 + ConvertString<int>(key) + 1 << " is attached to " << value << endl;
                else
                    cout << "Carbon at position " << (gmml::ConvertString<int>(key) + 1) << " is attached to " << value << endl;
            }
            else
            {
                if(key.compare("a") == 0)
                    cout << "Carbon at position 1 is attached to " << value << endl;
                else if(key.compare("+1") == 0 || key.compare("+2") == 0 || key.compare("+3") == 0)
                    cout << "Carbon at position " << mono->cycle_atoms_.size() - 1 + ConvertString<int>(key) << " is attached to " << value << endl;
                else
                    cout << "Carbon at position " << key << " is attached to " << value << endl;
            }
        }

        ///GENERATING COMPLETE NAME
        if(plus_sides.size() <= 1)
        {
            ///COMPLETE NAME GENERATION BASED ON DERIVATIVE MAP
            GenerateCompleteSugarName(mono);
        }
        else///UPDATING SIDE ATOMS
        {
            vector<string>::iterator index_it;
            if((index_it = find(mono->chemical_code_->right_up_.begin(), mono->chemical_code_->right_up_.end(), "+1")) != mono->chemical_code_->right_up_.end()){}
            else if((index_it = find(mono->chemical_code_->right_down_.begin(), mono->chemical_code_->right_down_.end(), "+1")) != mono->chemical_code_->right_down_.end()){}

            ///CHECKING R or S
            stringstream plus_one;
            string orientation = CalculateRSOrientations(mono->cycle_atoms_.at(mono->cycle_atoms_.size() - 2), plus_sides.at(0), plus_sides.at(1));
            plus_one << "+1" << orientation;
            (*index_it) = plus_one.str();
            if(plus_sides.size() == 3)
            {
                stringstream plus_two;
                orientation = CalculateRSOrientations(plus_sides.at(0), plus_sides.at(1), plus_sides.at(2));
                plus_two << "+2" << orientation;
                mono->chemical_code_->right_up_.push_back(plus_two.str());
                mono->chemical_code_->right_up_.push_back("+3");
            }

            ///UPDATING CHEMICAL CODE
            UpdateComplexSugarChemicalCode(mono);

            ///PRINTING SIDE ATOMS OF COMPLEX STRUCTURE
            cout << "Complex structure side group atoms: " << endl;
            for(vector<AtomVector>::iterator it1 = mono->side_atoms_.begin(); it1 != mono->side_atoms_.end(); it1++)
            {
                stringstream complex_structure_side;
                AtomVector sides = (*it1);
                if(it1 == mono->side_atoms_.begin())///side atoms of anomeric carbon
                {
                    if(sides.at(0) != NULL && sides.at(1) != NULL)
                        complex_structure_side << "[1] -> " << sides.at(0)->GetId() << ", " << sides.at(1)->GetId() << endl;
                    else if(sides.at(1) != NULL)
                        complex_structure_side << "[1] -> " << sides.at(1)->GetId() << endl;
                    else if(sides.at(0) != NULL)
                        complex_structure_side << "[1] -> " << sides.at(0)->GetId() << endl;
                }
                else if(it1 == mono->side_atoms_.end() - 1)///side atoms of last carbon of the ring
                {
                    complex_structure_side << "[" << mono->cycle_atoms_.size() - 1 << "]";
                    for(int i = 0; i < plus_sides.size() ; i++)
                        complex_structure_side << " -> " << sides.at(i)->GetId();
                    cout << complex_structure_side.str() << endl;
                }
                else if(sides.at(1) != NULL)
                {
                    int cycle_atom_index = distance(mono->side_atoms_.begin(), it1);
                    complex_structure_side << "[" << cycle_atom_index + 1 << "] -> " << sides.at(1)->GetId() << endl;
                }
            }

            ///PRINTING COMPLEX SUGAR CHEMICAL CODE
            cout << endl << "Complex sugar chemical code:" << endl;
            mono->chemical_code_->Print(cout);

            ///FINDING COMPLEX CHEMICAL CODE IN COMPLEX SUGAR NAME LOOKUP TABLE
            mono->sugar_name_ = ComplexSugarNameLookup(mono->chemical_code_->toString());

            if(plus_sides.size() == 2)
            {
                ///COMPLETE NAME GENERATION BASED ON DERIVATIVE MAP
                GenerateCompleteSugarName(mono);
            }
        }
        cout << endl;

        if(mono->sugar_name_.monosaccharide_stereochemistry_name_.compare("") == 0 && mono->sugar_name_.monosaccharide_name_.compare("") == 0)
        {
            ///FINDING CLOSEST MATCH FOR THE CHEMICAL CODE IN THE LOOKUP TABLE
            vector<Glycan::SugarName> closest_matches = vector<Glycan::SugarName>();
            mono->sugar_name_ = ClosestMatchSugarStereoChemistryNameLookup(mono->chemical_code_->toString(), closest_matches);
            if(mono->sugar_name_.monosaccharide_name_.compare("") == 0)
                mono->sugar_name_.monosaccharide_name_ = mono->sugar_name_.monosaccharide_stereochemistry_name_;
            if(mono->sugar_name_.monosaccharide_short_name_.compare("") == 0)
                mono->sugar_name_.monosaccharide_short_name_ = mono->sugar_name_.monosaccharide_stereochemistry_short_name_;

            ///ADDING NOTES/ISSUES OF MONOSACCHARIDE STRUCTURE
            Note* matching_note = new Note();
            matching_note->type_ = Glycan::COMMENT;
            matching_note->category_ = Glycan::MONOSACCHARIDE;
            stringstream ss;
            ss << "No exact match for " << mono->sugar_name_.monosaccharide_stereochemistry_short_name_ << ". close matches: ";
            for(vector<Glycan::SugarName>::iterator ite = closest_matches.begin(); ite != closest_matches.end(); ite++)
            {
                Glycan::SugarName sn = (*ite);
                if(ite == closest_matches.end() - 1)
                    ss << sn.monosaccharide_stereochemistry_short_name_;
                else
                    ss << sn.monosaccharide_stereochemistry_short_name_ << ", ";
            }
            matching_note->description_ = ss.str();
            this->AddNote(matching_note);


            if(mono->sugar_name_.monosaccharide_stereochemistry_name_.compare("") == 0)
            {
                mono->sugar_name_.monosaccharide_stereochemistry_name_ = "Unknown";
                mono->sugar_name_.monosaccharide_stereochemistry_short_name_ = "Unknown";
                mono->sugar_name_.monosaccharide_name_ = "Unknown";
                mono->sugar_name_.monosaccharide_short_name_ = "Unknown";
            }
            else
                cout << "No exact match found for the chemical code, the following information comes from one of the closest matches:" << endl;
        }

        ///GLYPROBITY REPORT (GEOMETRY OUTLIERS)
        if(glyprobity_report)
            CalculateGlyprobityGeometryOutliers(mono);

        ///ADDING NOTES/ISSUES OF ANOMERIC CONFIGURATION
        Note* anomeric_note = anomeric_notes.at(status_index);
        stringstream n;
        if(anomeric_note->description_.compare("") != 0)
        {
            n << mono->sugar_name_.monosaccharide_short_name_ << ": " << anomeric_note->description_;
            anomeric_note->description_ = n.str();
            this->AddNote(anomeric_note);
        }

        ///ADDING NOTES/ISSUES OF RESIDUE NAMING
        string original_residue = mono->cycle_atoms_.at(0)->GetResidue()->GetName();
        vector<string> pdb_codes = Split(mono->sugar_name_.pdb_code_, ",");
        if(pdb_codes.size() > 0)
        {
            string pdb_code = "";
            bool found_code = false;
            for(vector<string>::iterator codes_it = pdb_codes.begin(); codes_it != pdb_codes.end(); codes_it++)
            {
                pdb_code = (*codes_it);
                if(pdb_code.compare(original_residue) == 0)
                {
                    found_code = true;
                    break;
                }
            }
            if(!found_code)
            {
                Note* residue_naming_note = new Note();
                residue_naming_note->category_ = Glycan::RESIDUE_NAME;
                stringstream res_ss;
                if(mono->sugar_name_.pdb_code_.compare("") == 0)
                {
                    residue_naming_note->type_ = Glycan::WARNING;
                    res_ss << "PDB 3 letter code not found for " << mono->sugar_name_.monosaccharide_short_name_;
                }
                else
                {
                    residue_naming_note->type_ = Glycan::ERROR;
                    res_ss << "Residue name in input PDB file for " << mono->sugar_name_.monosaccharide_short_name_ << " does not match to PDB code: " << mono->sugar_name_.pdb_code_;
                }
                residue_naming_note->description_ = res_ss.str();
                this->AddNote(residue_naming_note);
            }
        }

        ///PRINTING NAMES OF MONOSACCHARIDE
        cout << "Stereochemistry name: " << mono->sugar_name_.monosaccharide_stereochemistry_name_ << endl;
        cout << "Stereochemistry short name: " << mono->sugar_name_.monosaccharide_stereochemistry_short_name_ << endl;
        cout << "Complete name: " << mono->sugar_name_.monosaccharide_name_ << endl;
        cout << "Short name: " << mono->sugar_name_.monosaccharide_short_name_ << endl;

        cout << "-------------------------------------------------------------------------------------------------------------------------------------------" << endl;

        mono_id++;
        mono->mono_id = mono_id;
        monos.push_back(mono);
    }

    ///CREATING TREE-LIKE STRUCTURE OF OLIGOSACCHARIDE
    cout << endl << "Oligosaccharides:" << endl;
    int number_of_covalent_links = 0;
    int number_of_probable_non_covalent_complexes = 0;
    vector<Oligosaccharide*> oligosaccharides = ExtractOligosaccharides(monos, dataset_residue_names, number_of_covalent_links, number_of_probable_non_covalent_complexes);

    ///BUILDING OLIGOSACCHARIDE SEQUENCE
    int number_of_oligosaccharides = 0;
    int number_of_monosaccharides = 0;
    for(vector<Oligosaccharide*>::iterator it = oligosaccharides.begin(); it != oligosaccharides.end(); it++)
    {
        Oligosaccharide* oligo = (*it);
        if(oligo->child_oligos_linkages_.size() > 0)
            number_of_oligosaccharides++;
        else
            number_of_monosaccharides++;
        oligo->Print(cout);
    }

    ///PRINTING NOTES AND ISSUES FOUND WITH THE INPUT FILE
    vector<Note*> notes = this->GetNotes();
    cout << "-------------------------------------------------------------------------------------------------------------------------------------------" << endl;
    cout << endl << "NOTES/ISSUES:" << endl;
    for(vector<Note*>::iterator note_it = notes.begin(); note_it != notes.end(); note_it++)
    {
        Note* note = (*note_it);
        cout << endl << "Category: " << note->ConvertGlycanNoteCat2String(note->category_) << endl;
        cout << "Type: " << note->ConvertGlycanNoteType2String(note->type_) << endl;
        cout << "Description: " << note->description_ << endl;
    }

    cout << "-------------------------------------------------------------------------------------------------------------------------------------------" << endl;

    ///PRINTING STATISTICAL REPORT OF GLYPROBITY
    if(glyprobity_report)
    {
        cout << endl << "GLYPROBITY REPORT" << endl;
        cout << "<-------Topology------>" << endl;
        cout << "Monosaccharide detected: " << monos.size() << endl;
        cout << "Residue Distribution " << endl;
        cout << " Monosaccharides: " << number_of_monosaccharides << endl;
        cout << " Oligosaccharides: " << number_of_oligosaccharides << endl;
        cout << "Carbohydrate Context " << endl;
        cout << " Covalently linked to protein: " << number_of_covalent_links << endl;
        cout << " Non-covalent complex: " << number_of_probable_non_covalent_complexes << endl;
        cout << "<--------------------->" << endl;
    }

    ///POPULATING GMMO ONTOLOGY
    if(populate_ontology)
    {
        if(oligosaccharides.size() > 0)
        {
            string gmmo = "gmmo.ttl";
            std::ofstream out_file;
            out_file.open(gmmo.c_str(), fstream::app);

            ifstream in("gmmo.ttl");///Checking if the file is empty
            size_t out_file_size = 0;
            in.seekg(0,ios_base::end);
            out_file_size = in.tellg();
            in.close();
            if(out_file_size == 0) ///If the file is empty add the prefixes first
                out_file << Ontology::TTL_FILE_PREFIX << endl;

            this->PopulateOntology(out_file, oligosaccharides);
            out_file.close();
        }
    }

    return oligosaccharides;
}

void Assembly::PopulateOntology(ofstream& main_stream, OligosaccharideVector oligos)
{
    stringstream pdb_stream;

    string pdb_resource = CreateURIResource(gmml::OntPDB, 0, "", "");
    //    CreateTitle(pdb_resource, pdb_stream);
    stringstream ss;
    ss << pdb_resource << "_";
    string id_prefix = ss.str();
    string pdb_uri = CreateURI(pdb_resource);

    //    pdb_stream << Ontology::ENTITY_COMMENT << pdb_resource << endl;
    AddTriple(pdb_uri, Ontology::TYPE, Ontology::PDB, pdb_stream);
    AddLiteral(pdb_uri, Ontology::id, pdb_resource, pdb_stream);
    //    AddLiteral(pdb_uri, Ontology::LABEL, pdb_resource, pdb_stream);
    //    AddLiteral(pdb_uri, Ontology::input_file_path, source_file_, pdb_stream);

    int link_id = 1;
    stringstream oligo_stream;
    stringstream mono_stream;
    stringstream linkage_stream;
    vector<string> side_or_ring_atoms = vector<string>();
    vector<int> visited_oligos = vector<int>();

    NoteVector notes = this->GetNotes();
    stringstream note_stream;
    if(notes.size() != 0)
    {
        int note_id = 1;
        PopulateNotes(pdb_stream, note_stream, pdb_uri, notes, id_prefix, note_id);
    }
    PopulateOligosaccharide(pdb_stream, oligo_stream, mono_stream, linkage_stream, pdb_uri, id_prefix, link_id, oligos, side_or_ring_atoms, visited_oligos);

    ResidueVector residues = this->GetResidues();
    stringstream residue_stream;
    PopulateResidue(pdb_stream, residue_stream, pdb_uri, id_prefix, residues, side_or_ring_atoms);

    main_stream << pdb_stream.str() << note_stream.str() << oligo_stream.str() << mono_stream.str() << linkage_stream.str() << residue_stream.str() << endl;
}
void Assembly::PopulateNotes(stringstream& pdb_stream, stringstream& note_stream, string pdb_uri, NoteVector notes, string id_prefix, int note_id)
{
    string note_resource = "";
    string note_uri = "";
    for(NoteVector::iterator it = notes.begin(); it != notes.end(); it++)
    {
        Note* note = (*it);
        note_resource = CreateURIResource(gmml::OntNote, note_id, id_prefix, "");
        note_uri = CreateURI(note_resource);
        AddTriple(pdb_uri, Ontology::hasNote, note_uri, pdb_stream);

        //        note_stream << Ontology::ENTITY_COMMENT << note_resource << endl;
        AddTriple(note_uri, Ontology::TYPE, Ontology::Note, note_stream);
        //        AddLiteral(note_uri, Ontology::LABEL, note_resource, note_stream);
        AddLiteral(note_uri, Ontology::note_type, note->ConvertGlycanNoteType2String(note->type_), note_stream);
        AddLiteral(note_uri, Ontology::note_category, note->ConvertGlycanNoteCat2String(note->category_), note_stream);
        AddLiteral(note_uri, Ontology::note_description, note->description_, note_stream);
        note_id++;
    }
}
void Assembly::PopulateOligosaccharide(stringstream& pdb_stream, stringstream& oligo_stream, stringstream& mono_stream, stringstream& linkage_stream, string pdb_uri, string id_prefix,
                                       int& link_id, OligosaccharideVector oligos, vector<string>& side_or_ring_atoms, vector<int>& visited_oligos)
{
    string oligo_resource = "";
    string oligo_uri = "";
    if(oligos.size() != NULL)
    {
        for(OligosaccharideVector::iterator it = oligos.begin(); it != oligos.end(); it++)
        {
            Oligosaccharide* oligo = (*it);

            oligo_resource = CreateURIResource(gmml::OntOligosaccharide, oligo->root_->mono_id, id_prefix, "");
            oligo_uri = CreateURI(oligo_resource);

            AddTriple(pdb_uri, Ontology::hasOligo, oligo_uri, pdb_stream);

            //            oligo_stream << Ontology::ENTITY_COMMENT << oligo_resource << endl;
            AddTriple(oligo_uri, Ontology::TYPE, Ontology::Oligosaccharide, oligo_stream);
            //            AddLiteral(oligo_uri, Ontology::LABEL, oligo_resource, oligo_stream);
            string o_name = oligo->oligosaccharide_name_;
            if(o_name.compare("") != 0)
                AddLiteral(oligo_uri, Ontology::oligo_name, o_name, oligo_stream);

            string o_residue_links = oligo->oligosaccharide_residue_linkages_;
            if(o_residue_links.compare("") != 0)
                AddLiteral(oligo_uri, Ontology::oligo_residue_linkages, o_residue_links, oligo_stream);

            if(oligo->child_oligos_.size() != 0 && (find(visited_oligos.begin(), visited_oligos.end(), oligo->root_->mono_id) == visited_oligos.end()))
            {
                PopulateLinkage(linkage_stream, oligo, oligo_uri, id_prefix, link_id, visited_oligos);
            }

            Monosaccharide* mono = oligo->root_;
            PopulateMonosaccharide(mono_stream, oligo_stream, oligo_uri, id_prefix, mono, side_or_ring_atoms);

            vector<Oligosaccharide*> child_oligos = oligo->child_oligos_;
            PopulateOligosaccharide(pdb_stream, oligo_stream, mono_stream, linkage_stream, pdb_uri, id_prefix, link_id, child_oligos, side_or_ring_atoms, visited_oligos);
        }
    }
}
void Assembly::PopulateLinkage(stringstream& linkage_stream, Oligosaccharide* oligo, string oligo_uri, string id_prefix, int& link_id, vector<int>& visited_oligos)
{
    string linkage_resource = "";
    string linkage_uri = "";
    string child_oligo_resource = "";
    string child_oligo_uri = "";
    string child_atom_resource = "";
    string child_atom_uri = "";
    string glycosidic_atom_resource = "";
    string glycosidic_atom_uri = "";
    string parent_atom_resource = "";
    string parent_atom_uri = "";
    stringstream linkage_str;
    stringstream glycosidic_linkage_str;


    visited_oligos.push_back(oligo->root_->mono_id);
    for(OligosaccharideVector::iterator it = oligo->child_oligos_.begin(); it != oligo->child_oligos_.end(); it++)
    {
        int index = distance(oligo->child_oligos_.begin(), it);

        Oligosaccharide* child_oligo = (*it);
        //        visited_oligos.push_back(child_oligo->root_->mono_id);

        linkage_resource = CreateURIResource(gmml::OntLinkage, link_id, id_prefix, "");
        linkage_uri = CreateURI(linkage_resource);
        //        linkage_stream << Ontology::ENTITY_COMMENT << linkage_resource << endl;
        AddTriple(linkage_uri, Ontology::TYPE, Ontology::Linkage, linkage_stream);
        //        AddLiteral(linkage_uri, Ontology::LABEL, linkage_resource, linkage_stream);
        link_id++;

        AddTriple(linkage_uri, Ontology::hasParent, oligo_uri, linkage_stream);
        child_oligo_resource = CreateURIResource(gmml::OntOligosaccharide, child_oligo->root_->mono_id, id_prefix, "");
        child_oligo_uri = CreateURI(child_oligo_resource);
        AddTriple(linkage_uri, Ontology::hasChild, child_oligo_uri, linkage_stream);

        vector<string> linkage_tokens = gmml::Split(oligo->child_oligos_linkages_.at(index), "-");
        string parent_atom_id = linkage_tokens.at(0);
        string glycosidic_atom_id = linkage_tokens.at(1);
        string child_atom_id = linkage_tokens.at(2);

        int parent_c_index = ExtractLinkageCarbonIndex(oligo, parent_atom_id);
        int child_c_index = ExtractLinkageCarbonIndex(child_oligo, child_atom_id);

        if(child_c_index != 0 && parent_c_index != 0)
        {
            stringstream link_indeces_str;
            link_indeces_str << child_c_index << "-" << parent_c_index;
            AddLiteral(linkage_uri, Ontology::linkageIndeces, link_indeces_str.str(), linkage_stream);
        }

        child_atom_resource = CreateURIResource(gmml::OntAtom, 0, id_prefix, child_atom_id);
        child_atom_uri = CreateURI(child_atom_resource);
        AddTriple(linkage_uri, Ontology::hasChildAtomLinkage, child_atom_uri, linkage_stream);

        glycosidic_atom_resource = CreateURIResource(gmml::OntAtom, 0, id_prefix, glycosidic_atom_id);
        glycosidic_atom_uri = CreateURI(glycosidic_atom_resource);
        AddTriple(linkage_uri, Ontology::hasGlycosidicLinkage, glycosidic_atom_uri, linkage_stream);

        parent_atom_resource = CreateURIResource(gmml::OntAtom, 0, id_prefix, parent_atom_id);
        parent_atom_uri = CreateURI(parent_atom_resource);
        AddTriple(linkage_uri, Ontology::hasParentAtomLinkage, parent_atom_uri, linkage_stream);

        vector<string> child_atom_id_tokens = gmml::Split(child_atom_id, "_");
        if(child_atom_id_tokens.at(3).at(0) == gmml::BLANK_SPACE)
            linkage_str << child_atom_id_tokens.at(2) << "(" << child_atom_id_tokens.at(4) << ")" << child_atom_id_tokens.at(0);
        else
            linkage_str << child_atom_id_tokens.at(2) << "(" << child_atom_id_tokens.at(4) << "_" << child_atom_id_tokens.at(3) << ")" << child_atom_id_tokens.at(0);

        std::vector<std::string> parent_atom_id_tokens = gmml::Split(parent_atom_id, "_");
        if(parent_atom_id_tokens.at(3).at(0) == gmml::BLANK_SPACE)
            linkage_str << "-" << parent_atom_id_tokens.at(2) << "(" << parent_atom_id_tokens.at(4) << ")"  << parent_atom_id_tokens.at(0);
        else
            linkage_str << "-" << parent_atom_id_tokens.at(2) << "(" << parent_atom_id_tokens.at(4) <<  "_" << parent_atom_id_tokens.at(3) << ")"  << parent_atom_id_tokens.at(0);

        //        AddLiteral(linkage_uri, Ontology::linkage_str, linkage_str.str(), linkage_stream);

        std::vector<std::string> glycosidic_atom_id_tokens = gmml::Split(glycosidic_atom_id, "_");
        if(glycosidic_atom_id_tokens.at(3).at(0) == gmml::BLANK_SPACE)
            glycosidic_linkage_str << glycosidic_atom_id_tokens.at(2) << "(" << glycosidic_atom_id_tokens.at(4) << ")" << glycosidic_atom_id_tokens.at(0);
        else
            glycosidic_linkage_str << glycosidic_atom_id_tokens.at(2) << "(" << glycosidic_atom_id_tokens.at(4) << "_" << glycosidic_atom_id_tokens.at(3)
                                   << ")"  << glycosidic_atom_id_tokens.at(0);

        AddLiteral(linkage_uri, Ontology::glycosidic_linkage, glycosidic_linkage_str.str(), linkage_stream);
    }
}
int Assembly::ExtractLinkageCarbonIndex(Oligosaccharide* oligo, string linkage_carbon_id)
{
    int c_index = 0;
    vector<string> cycle_atom_tokens = Split(oligo->root_->cycle_atoms_str_, "-");

    if(oligo->root_->side_atoms_.at(0).at(0) != NULL)
    {
        c_index++;
        Atom* anomeric_side_carbon = oligo->root_->side_atoms_.at(0).at(0);
        if(anomeric_side_carbon->GetId().compare(linkage_carbon_id) == 0)
            return c_index;
    }
    for(int i = 0; i < cycle_atom_tokens.size() - 1; i++) /// cycle_atom_tokens.size() - 1 > because the ring oxygen is not considered
    {
        c_index++;
        if(cycle_atom_tokens.at(i).compare(linkage_carbon_id) == 0)
            return c_index;
    }

    AtomVector side_atoms_of_last_ring_carbon = oligo->root_->side_atoms_.at(oligo->root_->side_atoms_.size() - 1);
    for(AtomVector::iterator it1 = side_atoms_of_last_ring_carbon.begin(); it1 != side_atoms_of_last_ring_carbon.end(); it1++)
    {
        Atom* side_atom = (*it1);
        c_index++;

        if(side_atom->GetId().compare(linkage_carbon_id) == 0)
            return c_index;
    }
    return c_index;
}
void Assembly::PopulateMonosaccharide(stringstream& mono_stream, stringstream& oligo_stream, string oligo_uri, string id_prefix, Monosaccharide* mono,
                                      vector<std::string>& side_or_ring_atoms)
{
    stringstream object;
    string mono_resource = "";
    string mono_uri = "";
    string ring_resource = "";
    string ring_uri = "";

    mono_resource = CreateURIResource(gmml::OntMonosaccharide, mono->mono_id, id_prefix, "");
    mono_uri = CreateURI(mono_resource);

    AddTriple(oligo_uri, Ontology::hasCore, mono_uri, oligo_stream);

    //    mono_stream << Ontology::ENTITY_COMMENT << mono_resource << endl;
    AddTriple(mono_uri, Ontology::TYPE, Ontology::Monosaccharide, mono_stream);
    AddLiteral(mono_uri, Ontology::id, mono_resource, mono_stream);
    //    AddLiteral(mono_uri, Ontology::LABEL, mono_resource, mono_stream);

    AtomVector ring_atoms = mono->cycle_atoms_;
    object.str(string());
    int ring_index = 1;
    stringstream ring_atom_stream;
    for(AtomVector::iterator it = ring_atoms.begin(); it != ring_atoms.end(); it++)
    {
        Atom* ring_atom = (*it);

        ring_resource = CreateURIResource(gmml::OntAtom, 0, id_prefix, ring_atom->GetId());
        ring_uri = CreateURI(ring_resource);
        AddTriple(mono_uri, Ontology::hasRingAtom, ring_uri, mono_stream);

        PopulateRingAtom(ring_atom_stream, id_prefix, ring_uri, ring_resource, ring_index, ring_atom, mono, side_or_ring_atoms);
        ring_index++;

        if(it == ring_atoms.end() - 1)
            object << ring_resource;
        else
            object << ring_resource << "-";
    }
    AddLiteral(mono_uri, Ontology::ring_atoms, object.str(), mono_stream);

    object.str(string());
    object << mono->anomeric_status_ << " " << CreateURIResource(gmml::OntAtom, 0, id_prefix, mono->cycle_atoms_.at(0)->GetId());
    AddLiteral(mono_uri, Ontology::anomeric_status, object.str(), mono_stream);

    AddLiteral(mono_uri, Ontology::stereochemistry_chemical_code, mono->sugar_name_.chemical_code_string_, mono_stream);
    if(mono->bfmp_ring_conformation_.compare("") != 0)
        AddLiteral(mono_uri, Ontology::bfmp_ring_conformation, mono->bfmp_ring_conformation_, mono_stream);

    SugarName sugar_name = mono->sugar_name_;
    PopulateSugarName(mono_stream, id_prefix, mono_uri, mono->mono_id, sugar_name);
    mono_stream << ring_atom_stream.str();

}
void Assembly::PopulateRingAtom(stringstream& ring_atom_stream, string id_prefix, string ring_uri, string ring_resource, int ring_index, Atom* ring_atom, Monosaccharide* mono,
                                vector<string>& side_or_ring_atoms)
{
    stringstream object;
    //    ring_atom_stream << Ontology::ENTITY_COMMENT << ring_resource << endl;
    AddTriple(ring_uri, Ontology::TYPE, Ontology::RingAtom, ring_atom_stream);
    object << ring_index;
    AddLiteral(ring_uri, Ontology::ring_index, object.str(), ring_atom_stream);
    AddLiteral(ring_uri, Ontology::id, ring_resource, ring_atom_stream);
    //    AddLiteral(ring_uri, Ontology::LABEL, ring_resource, ring_atom_stream);
    Coordinate* coords = ring_atom->GetCoordinates().at(model_index_);
    /*
    AddLiteral(ring_uri, Ontology::x, ConvertT<double>(coords->GetX()), ring_atom_stream);
    AddLiteral(ring_uri, Ontology::y, ConvertT<double>(coords->GetY()), ring_atom_stream);
    AddLiteral(ring_uri, Ontology::z, ConvertT<double>(coords->GetZ()), ring_atom_stream);
    */
    stringstream coord_stream;
    coord_stream << ConvertT<double>(coords->GetX()) << ", " << ConvertT<double>(coords->GetY()) << ", " << ConvertT<double>(coords->GetZ());
    AddLiteral(ring_uri, Ontology::coordinate, coord_stream.str(), ring_atom_stream);

    side_or_ring_atoms.push_back(ring_atom->GetId());

    string neighbor_resource = "";
    string neighbor_uri = "";
    AtomVector neighbors = ring_atom->GetNode()->GetNodeNeighbors();
    for(AtomVector::iterator it = neighbors.begin(); it != neighbors.end(); it++)
    {
        Atom* neighbor = (*it);
        neighbor_resource = CreateURIResource(gmml::OntAtom, 0, id_prefix, neighbor->GetId());
        neighbor_uri = CreateURI(neighbor_resource);
        AddTriple(ring_uri, Ontology::hasNeighbor, neighbor_uri, ring_atom_stream);
    }

    stringstream side_atom_stream;
    if((ring_atom->GetName().substr(0,1).compare("O") != 0 )) ///side atoms for the oxygen of the ring are not saved
    {
        vector<AtomVector> all_sides = mono->side_atoms_;
        AtomVector sides = all_sides.at(ring_index - 1);
        string side_resource = "";
        string side_uri = "";
        int side_index = 1;
        for(AtomVector::iterator it = sides.begin(); it != sides.end(); it++)
        {
            Atom* side_atom = (*it);
            if(side_atom != NULL)
            {
                side_resource = CreateURIResource(gmml::OntAtom, 0, id_prefix, side_atom->GetId());
                side_uri = CreateURI(side_resource);
                AddTriple(ring_uri, Ontology::hasSideAtom, side_uri, ring_atom_stream);

                PopulateSideAtom(side_atom_stream, id_prefix, side_uri, side_resource, ring_index, side_index, side_atom, mono, side_or_ring_atoms);
                side_index++;
            }
        }
    }
    ring_atom_stream << side_atom_stream.str();
}
void Assembly::PopulateSideAtom(stringstream& side_atom_stream, string id_prefix, string side_uri, string side_resource, int ring_index, int side_index, Atom* side_atom, Monosaccharide* mono,
                                vector<string>& side_or_ring_atoms)
{
    stringstream object;

    string chemical_code_str = mono->sugar_name_.chemical_code_string_;
    if((mono->sugar_name_.ring_type_.compare("P") == 0 && ring_index == 5) || (mono->sugar_name_.ring_type_.compare("F") == 0 && ring_index == 4))
        object << "+" << side_index;
    else if(ring_index == 1 && (side_atom->GetName().substr(0,1).compare("C") != 0))
        object << "1";
    else if(ring_index == 1 && (side_atom->GetName().substr(0,1).compare("C") == 0 ))
        object << "-1";
    else
        object << ring_index;

    if(find(side_or_ring_atoms.begin(), side_or_ring_atoms.end(), side_atom->GetId()) == side_or_ring_atoms.end())///if this side atom has not been added to the ontology as side atom of another mono
    {
        //    side_atom_stream << Ontology::ENTITY_COMMENT << side_resource << endl;
        AddTriple(side_uri, Ontology::TYPE, Ontology::SideAtom, side_atom_stream);
        AddLiteral(side_uri, Ontology::id, side_resource, side_atom_stream);
        //        AddLiteral(side_uri, Ontology::LABEL, side_resource, side_atom_stream);
        Coordinate* coords = side_atom->GetCoordinates().at(model_index_);
        /*AddLiteral(side_uri, Ontology::x, ConvertT<double>(coords->GetX()), side_atom_stream);
        AddLiteral(side_uri, Ontology::y, ConvertT<double>(coords->GetY()), side_atom_stream);
        AddLiteral(side_uri, Ontology::z, ConvertT<double>(coords->GetZ()), side_atom_stream);*/

        stringstream coord_stream;
        coord_stream << ConvertT<double>(coords->GetX()) << ", " << ConvertT<double>(coords->GetY()) << ", " << ConvertT<double>(coords->GetZ());
        AddLiteral(side_uri, Ontology::coordinate, coord_stream.str(), side_atom_stream);

        if(object.str().compare("1") == 0)
        {
            if(mono->derivatives_map_.find("a") != mono->derivatives_map_.end())
                AddLiteral(side_uri, Ontology::derivative, mono->derivatives_map_[object.str()], side_atom_stream);
        }
        else
        {
            if(mono->derivatives_map_.find(object.str()) != mono->derivatives_map_.end())
                AddLiteral(side_uri, Ontology::derivative, mono->derivatives_map_[object.str()], side_atom_stream);
        }

        string neighbor_resource = "";
        string neighbor_uri = "";
        AtomVector neighbors = side_atom->GetNode()->GetNodeNeighbors();
        for(AtomVector::iterator it = neighbors.begin(); it != neighbors.end(); it++)
        {
            Atom* neighbor = (*it);
            neighbor_resource = CreateURIResource(gmml::OntAtom, 0, id_prefix, neighbor->GetId());
            neighbor_uri = CreateURI(neighbor_resource);
            AddTriple(side_uri, Ontology::hasNeighbor, neighbor_uri, side_atom_stream);
        }

        side_or_ring_atoms.push_back(side_atom->GetId());
    }

    size_t index = string::npos;
    if(object.str().compare("1") == 0)
        index = chemical_code_str.find("a");
    else
        index = chemical_code_str.find(object.str().c_str());
    if (index != string::npos)
    {
        index--;
        if(chemical_code_str.at(index) == '^')
        {
            AddLiteral(side_uri, Ontology::orientation, "Up", side_atom_stream);
        }
        else if(chemical_code_str.at(index) == '_')
        {
            AddLiteral(side_uri, Ontology::orientation, "Down", side_atom_stream);
        }
        else
        {
            AddLiteral(side_uri, Ontology::orientation, "Deoxy", side_atom_stream);
        }
        AddLiteral(side_uri, Ontology::side_index, object.str(), side_atom_stream);
    }

}
void Assembly::PopulateSugarName(stringstream& mono_stream, string id_prefix, string mono_uri, int mono_id, SugarName sugar_name)
{
    stringstream sugar_name_stream;

    string sugar_name_resource = "";
    string sugar_name_uri = "";

    sugar_name_resource = CreateURIResource(gmml::OntSugarName, mono_id, id_prefix, "");
    sugar_name_uri = CreateURI(sugar_name_resource);

    AddTriple(mono_uri, Ontology::hasSugarName, sugar_name_uri, mono_stream);

    //    sugar_name_stream << Ontology::ENTITY_COMMENT << sugar_name_resource << endl;
    AddTriple(sugar_name_uri, Ontology::TYPE, Ontology::SugarName, sugar_name_stream);
    AddLiteral(sugar_name_uri, Ontology::id, sugar_name_resource, sugar_name_stream);
    //    AddLiteral(sugar_name_uri, Ontology::LABEL, sugar_name_resource, sugar_name_stream);
    AddLiteral(sugar_name_uri, Ontology::mono_stereo_name, sugar_name.monosaccharide_stereochemistry_name_, sugar_name_stream);
    AddLiteral(sugar_name_uri, Ontology::mono_stereo_short_name, sugar_name.monosaccharide_stereochemistry_short_name_, sugar_name_stream);
    AddLiteral(sugar_name_uri, Ontology::mono_name, sugar_name.monosaccharide_name_, sugar_name_stream);
    AddLiteral(sugar_name_uri, Ontology::mono_short_name, sugar_name.monosaccharide_short_name_, sugar_name_stream);
    AddLiteral(sugar_name_uri, Ontology::isomer, sugar_name.isomer_, sugar_name_stream);
    if(sugar_name.configuration_.compare("a") == 0)
        AddLiteral(sugar_name_uri, Ontology::configuration, "alpha", sugar_name_stream);
    else if(sugar_name.configuration_.compare("b") == 0)
        AddLiteral(sugar_name_uri, Ontology::configuration, "beta", sugar_name_stream);
    if(sugar_name.ring_type_.compare("") != 0)
        AddLiteral(sugar_name_uri, Ontology::ring_type, sugar_name.ring_type_, sugar_name_stream);

    mono_stream << sugar_name_stream.str();
}
void Assembly::PopulateResidue(stringstream& pdb_stream, stringstream& residue_stream, string pdb_uri, string id_prefix, ResidueVector residues, vector<string> side_or_ring_atoms)
{
    string res_resource = "";
    string res_uri = "";
    string atom_resource = "";
    string atom_uri = "";

    stringstream atom_stream;
    for(ResidueVector::iterator it = residues.begin(); it != residues.end(); it++)
    {
        Residue* residue = (*it);
        res_resource = CreateURIResource(gmml::OntResidue, 0, id_prefix, residue->GetId());
        res_uri = CreateURI(res_resource);

        //        residue_stream << Ontology::ENTITY_COMMENT << res_resource << endl;
        AddTriple(res_uri, Ontology::TYPE, Ontology::Residue, residue_stream);
        AddLiteral(res_uri, Ontology::id, res_resource, residue_stream);
        //        AddLiteral(res_uri, Ontology::LABEL, res_resource, residue_stream);

        AtomVector res_atoms = residue->GetAtoms();
        for(AtomVector::iterator it1 = res_atoms.begin(); it1 != res_atoms.end(); it1++)
        {
            Atom* atom = (*it1);
            atom_resource = CreateURIResource(gmml::OntAtom, 0, id_prefix, atom->GetId());
            atom_uri = CreateURI(atom_resource);

            AddTriple(res_uri, Ontology::hasAtom, atom_uri, residue_stream);

            if(find(side_or_ring_atoms.begin(), side_or_ring_atoms.end(), atom->GetId()) == side_or_ring_atoms.end())
            {
                PopulateAtom(atom_stream, atom_uri, atom_resource, id_prefix, atom);
            }
        }
        AddTriple(pdb_uri, Ontology::hasResidue, res_uri, pdb_stream);
    }
    residue_stream << atom_stream.str();
}
void Assembly::PopulateAtom(stringstream& atom_stream, string atom_uri, string atom_resource, string id_prefix, Atom* atom)
{
    //    atom_stream << Ontology::ENTITY_COMMENT << atom_resource << endl;
    AddTriple(atom_uri, Ontology::TYPE, Ontology::Atom, atom_stream);
    AddLiteral(atom_uri, Ontology::id, atom_resource, atom_stream);
    //    AddLiteral(atom_uri, Ontology::LABEL, atom_resource, atom_stream);
    Coordinate* coords = atom->GetCoordinates().at(model_index_);
    /*AddLiteral(atom_uri, Ontology::x, ConvertT<double>(coords->GetX()), atom_stream);
    AddLiteral(atom_uri, Ontology::y, ConvertT<double>(coords->GetY()), atom_stream);
    AddLiteral(atom_uri, Ontology::z, ConvertT<double>(coords->GetZ()), atom_stream);*/

    stringstream coord_stream;
    coord_stream << ConvertT<double>(coords->GetX()) << ", " << ConvertT<double>(coords->GetY()) << ", " << ConvertT<double>(coords->GetZ());
    AddLiteral(atom_uri, Ontology::coordinate, coord_stream.str(), atom_stream);

    string neighbor_resource = "";
    string neighbor_uri = "";
    AtomVector neighbors = atom->GetNode()->GetNodeNeighbors();
    for(AtomVector::iterator it = neighbors.begin(); it != neighbors.end(); it++)
    {
        Atom* neighbor = (*it);
        neighbor_resource = CreateURIResource(gmml::OntAtom, 0, id_prefix, neighbor->GetId());
        neighbor_uri = CreateURI(neighbor_resource);
        AddTriple(atom_uri, Ontology::hasNeighbor, neighbor_uri, atom_stream);
    }
}
void Assembly::CreateTitle(string pdb_resource, stringstream& pdb_stream)
{
    pdb_stream << endl << "####################################" << endl;
    pdb_stream << "#" << setw(9) << " " << pdb_resource << " Individuals" << setw(9) << " " << "#" << endl;
    pdb_stream << "####################################" << endl;
}
void Assembly::AddTriple(string s, string p, string o, stringstream& stream)
{
    stream << s << " " << p << " " << o << "." << endl;
}
void Assembly::AddLiteral(string s, string p, string o, stringstream& stream)
{
    stream << s << " " << p << " \"" << o << "\"." << endl;
}
string Assembly::CreateURI(string uri_resource)
{
    stringstream uri;
    uri << Ontology::ONT_PREFIX << uri_resource;
    return uri.str();
}
string Assembly::CreateURIResource(gmml::URIType resource , int number, string id_prefix, string id )
{
    stringstream uri_resource;
    switch(resource)
    {
        case gmml::OntPDB:
            uri_resource << (Split(this->GetSourceFile().substr(this->GetSourceFile().find_last_of('/') + 1), ".").at(0));
            break;
        case gmml::OntOligosaccharide:
            uri_resource << id_prefix << "oligo" << number;
            break;
        case gmml::OntNote:
            uri_resource << id_prefix << "note" << number;
            break;
        case gmml::OntLinkage:
            uri_resource << id_prefix << "link" << number;
            break;
        case gmml::OntMonosaccharide:
            uri_resource << id_prefix << "mono" << number;
            break;
        case gmml::OntSugarName:
            uri_resource << id_prefix << "mono" << number << "_sn";
            break;
        case gmml::OntAtom:
            replace( id.begin(), id.end(), '?', 'n'); // replace all '?' with 'n'
            FindReplaceString(id, "\'", "q");
            FindReplaceString(id, ",", "c");
            //FindReplaceString(id, "*", "s");
            replace( id.begin(), id.end(), '*', 's'); // replace all '*' with 's'
            uri_resource << id_prefix << id;
            break;
        case gmml::OntResidue:
            replace( id.begin(), id.end(), '?', 'n'); // replace all '?' with 'n'
            replace( id.begin(), id.end(), '*', 's'); // replace all '*' with 's'
            uri_resource << id_prefix << id;
            break;
    }
    return uri_resource.str();
}
void Assembly::FormulateCURL(string output_file_type, string query)
{
    cout << "GENERATED QUERY:" << endl;
    cout << query << endl;
    stringstream curl;
    curl << Ontology::CURL_PREFIX;

    if(output_file_type.compare("csv") == 0)
        curl << Ontology::CSV_OUTPUT_FORMAT;
    else if(output_file_type.compare("xml") == 0)
        curl << Ontology::XML_OUTPUT_FORMAT;
    else if(output_file_type.compare("json") == 0)
        curl << Ontology::JSON_OUTPUT_FORMAT;

    curl << Ontology::DATA_STORE_ADDRESS << Ontology::QUERY_PREFIX << query << Ontology::QUERY_POSTFIX;
    string tmp = curl.str();
    cout << endl << "RESULTS: " << endl;
    const char* cstr = tmp.c_str();
    system(cstr);
    cout << endl;
}
void Assembly::ExtractOntologyInfoByNameOfGlycan(string stereo_name, string stereo_condensed_name, string name, string condensed_name, string output_file_type)
{
    if(stereo_name.compare("") == 0 && stereo_condensed_name.compare("") == 0 && name.compare("") == 0 && condensed_name.compare("") == 0)
    {
        cout << "Please specify at least one of the arguments and set the others as \"\" " << endl;
        return;
    }
    stringstream query;
    query << Ontology::PREFIX << Ontology::SELECT_CLAUSE << " ?pdb ?residue_id ?stereo_name ?stereo_condensed_name ?name ?condensed_name ?shape " << Ontology::WHERE_CLAUSE;
    query << "?pdb_file     :hasOligo   ?oligo.\n";
    query << "?oligo        :hasCore    ?mono.\n";
    query << "?mono         :hasSugarName   ?sugarName.\n";
    if(stereo_name.compare("") != 0)
        query << "?sugarName    :monosaccharideStereochemName   \"" << stereo_name << "\".\n";
    if(stereo_condensed_name.compare("") != 0)
        query << "?sugarName    :monosaccharideStereochemShortName   \"" << stereo_condensed_name << "\".\n";
    if(name.compare("") != 0)
        query << "?sugarName    :monosaccharideName   \"" << name << "\".\n";
    if(condensed_name.compare("") != 0)
        query << "?sugarName    :monosaccharideShortName   \"" << condensed_name << "\".\n";
    query << "?pdb_file     :identifier   ?pdb.\n";

    query << "?mono         :hasRingAtom   ?ring_atoms.\n";
    query << "?residue      :hasAtom       ?ring_atoms.\n";
    query << "?residue      :identifier    ?residue_id.\n";

    query << "?sugarName    :monosaccharideStereochemName   ?stereo_name.\n";
    query << "?sugarName    :monosaccharideStereochemShortName   ?stereo_condensed_name.\n";
    query << "?sugarName    :monosaccharideName   ?name.\n";
    query << "?sugarName    :monosaccharideShortName   ?condensed_name.\n";
    query << "OPTIONAL { ?mono         :BFMPRingConformation   ?shape.}\n";
    query << Ontology::END_WHERE_CLAUSE;

    FormulateCURL(output_file_type, query.str());
}
void Assembly::ExtractOntologyInfoByNamePartsOfGlycan(string isomer, string ring_type, string configuration, string output_file_type)
{
    if(isomer.compare("") == 0 && ring_type.compare("") == 0 && configuration.compare("") == 0)
    {
        cout << "Please specify at least one of the arguments and set the others as \"\" " << endl;
        return;
    }
    stringstream query;
    query << Ontology::PREFIX << Ontology::SELECT_CLAUSE << " ?pdb ?stereo_name ?stereo_condensed_name ?name ?condensed_name ?shape " << Ontology::WHERE_CLAUSE;
    query << "?pdb_file     :hasOligo   ?oligo.\n";
    query << "?oligo        :hasCore    ?mono.\n";
    query << "?mono         :hasSugarName   ?sugarName.\n";
    if(isomer.compare("") != 0)
        query << "?sugarName    :isomer   \"" << isomer << "\".\n";
    if(ring_type.compare("") != 0)
        query << "?sugarName    :ringType   \"" << ring_type << "\".\n";
    if(configuration.compare("") != 0)
        query << "?sugarName    :configuration   \"" << configuration << "\".\n";
    query << "?pdb_file     :identifier   ?pdb.\n";

//    query << "?mono         :hasRingAtom   ?ring_atoms.\n";
//    query << "?residue      :hasAtom       ?ring_atoms.\n";
//    query << "?residue      :identifier    ?residue_id.\n";

    query << "?sugarName    :monosaccharideStereochemName   ?stereo_name.\n";
    query << "?sugarName    :monosaccharideStereochemShortName   ?stereo_condensed_name.\n";
    query << "?sugarName    :monosaccharideName   ?name.\n";
    query << "?sugarName    :monosaccharideShortName   ?condensed_name.\n";
    query << "OPTIONAL { ?mono         :BFMPRingConformation   ?shape.}\n";
    query << Ontology::END_WHERE_CLAUSE;

    FormulateCURL(output_file_type, query.str());
}
void Assembly::ExtractOntologyInfoByPDBID(string pdb_id, string output_file_type)
{
    if(pdb_id.compare("") == 0)
    {
        cout << "Please specify the input argument." << endl;
        return;
    }
    stringstream query;
    query << Ontology::PREFIX << Ontology::SELECT_CLAUSE << " ?oligo_sequence ?residue_links ?glycosidic_linkage " << Ontology::WHERE_CLAUSE;
    query <<  ":" << pdb_id << "    :hasOligo   ?oligo.\n";
    query << "?oligo    :oligoName 	?oligo_sequence.\n";
    query << "OPTIONAL { ?oligo	:oligoResidueLinks	?residue_links.\n";
    query << "?linkage 	:hasParent 	?oligo.\n";
    query << "?linkage	:glycosidicLinkage    ?glycosidic_linkage.}\n";

    //query << "?oligo	:hasCore	?mono.\n";
    //query << "?mono     :anomericStatus    ?anomeric_status.\n";

    query << Ontology::END_WHERE_CLAUSE;

    FormulateCURL(output_file_type, query.str());
}
void Assembly::ExtractOntologyInfoByStringChemicalCode(string chemical_code, string output_file_type)
{
    if(chemical_code.compare("") == 0)
    {
        cout << "Please specify the input argument." << endl;
        return;
    }
    stringstream query;
    query << Ontology::PREFIX << Ontology::SELECT_CLAUSE << " ?pdb ?name ?short_name ?stereo_name ?stereo_short_name " << Ontology::WHERE_CLAUSE;
    query << "?mono         :stereochemistryChemicalCode	   \"" << chemical_code << "\".\n";
    query << "?pdb_file     :hasOligo	?oligo.\n";
    query << "?oligo        :hasCore	?mono.\n";
    query << "?mono         :hasSugarName	?sn.\n";
    query << "?sn           :monosaccharideName 	?name.\n";
    query << "?sn           :monosaccharideShortName 	?short_name.\n";
    query << "?sn           :monosaccharideStereochemName 	?stereo_name.\n";
    query << "?sn           :monosaccharideStereochemShortName 	?stereo_short_name.\n";
    query << "?pdb_file     :identifier   ?pdb.\n";
    query << Ontology::END_WHERE_CLAUSE;

    FormulateCURL(output_file_type, query.str());
}
void Assembly::ExtractOntologyInfoByOligosaccharideNameSequence(string oligo_name, string output_file_type)
{
    stringstream query;
    query << Ontology::PREFIX << Ontology::SELECT_CLAUSE << " ?pdb ?residue_links ?glycosidic_linkage " << Ontology::WHERE_CLAUSE;

    query << "?pdb_file     :hasOligo	?oligo.\n";
    query << "?oligo        :oligoName	\"" << oligo_name << "\"\n";
    query << "OPTIONAL { ?oligo        :oligoResidueLinks	?residue_links.\n";
    query << "?linkage      :hasParent 	?oligo.\n";
    query << "?linkage      :glycosidicLinkage    ?glycosidic_linkage.}\n";
    query << "?pdb_file     :identifier   ?pdb.\n";

    ///To DO: string manipulation: split the oligo_name by _ and for each of them write the following to represent the names of the monos:
    //    query << "?oligo	:hasCore	?mono.\n";
    //    query << "?mono     :hasSugarName	?sn.\n";
    //    query << "?sn       :monosaccharideName 	?name.\n";
    //    query << "?sn       :monosaccharideShortName 	?short_name.\n";
    //    query << "?sn       :monosaccharideStereochemName 	?stereo_name.\n";
    //    query << "?sn       :monosaccharideStereochemShortName 	?stereo_short_name.\n";
    query << Ontology::END_WHERE_CLAUSE;

    FormulateCURL(output_file_type, query.str());
}
void Assembly::ExtractOntologyInfoByOligosaccharideNameSequenceByRegex(string oligo_name_pattern, string output_file_type)
{
    FindReplaceString(oligo_name_pattern, "[", "\\\\[");
    FindReplaceString(oligo_name_pattern, "]", "\\\\]");
    if(oligo_name_pattern.compare("") == 0)
    {
        cout << "Please specify the input argument. (you can use up to two * in the name pattern)" << endl;
        return;
    }
    if(count(oligo_name_pattern.begin(), oligo_name_pattern.end(), '*') > 3)
    {
        cout << "Wrong name pattern format. Please use only up tp three * in the input argument." << endl;
        return;
    }

    stringstream query;
    query << Ontology::PREFIX << Ontology::SELECT_CLAUSE << " ?pdb ?oligo_sequence ?residue_links ?glycosidic_linkage " << Ontology::WHERE_CLAUSE;
    query << "?oligo        :oligoName	?oligo_sequence.\n";

    size_t first = oligo_name_pattern.find_first_of("*");
    size_t last = oligo_name_pattern.find_last_of("*");

    string filter1 = oligo_name_pattern.substr(0, first);
    string filter2 = oligo_name_pattern.substr(first + 1, last - 1);
    string filter3 = oligo_name_pattern.substr(last + 1, oligo_name_pattern.size() - 1);
    if(count(oligo_name_pattern.begin(), oligo_name_pattern.end(), '*') == 0) ///No *
        query << "?oligo	:oligoName	\"" << oligo_name_pattern << "\".\n";
    else if(count(oligo_name_pattern.begin(), oligo_name_pattern.end(), '*') == 1) ///Only one *
    {
        if(first == 0) ///* at the beginning
            query << "FILTER regex(?oligo_sequence, \"" << filter2 << "$\", \"i\")\n";
        else if(first == oligo_name_pattern.size()-1) ///* at the end
            query << "FILTER regex(?oligo_sequence, \"^" << filter1 << "\", \"i\")\n";
        else if(first < oligo_name_pattern.size()-1 && first > 0) ///* in the middle
            query << "FILTER regex(?oligo_sequence, \"^" << filter1 << ".+" << filter3 << "$\", \"i\")\n";
    }
    else if (count(oligo_name_pattern.begin(), oligo_name_pattern.end(), '*') == 2)
    {
        if(first == 0 && last == oligo_name_pattern.size() - 1) ///* at the beginning and end
            query << "FILTER regex(?oligo_sequence, \"" << filter2 << "\", \"i\")\n";
        else if(first == 0 && last < oligo_name_pattern.size() - 1)///one * at the beginning another in the middle
            query << "FILTER regex(?oligo_sequence, \"" << filter2 << ".+" << filter3 << "$\", \"i\")\n";
        else if(first > 0 && last == oligo_name_pattern.size() - 1)///one * in the middle another at the end
            query << "FILTER regex(?oligo_sequence, \"^" << filter1 << ".+" << filter2 << "\", \"i\")\n";
    }
    else
    {
        vector<string> pattern_tokens = Split(filter2, "*");
        query << "FILTER regex(?oligo_sequence, \"" << pattern_tokens.at(0) << ".+" << pattern_tokens.at(1) << "\", \"i\")\n";
    }

    query << "?pdb_file     :hasOligo	?oligo.\n";
    query << "OPTIONAL { ?oligo        :oligoResidueLinks	?residue_links.\n";
    query << "?linkage      :hasParent 	?oligo.\n";
    query << "?linkage      :glycosidicLinkage    ?glycosidic_linkage.}\n";
    query << "?pdb_file     :identifier	?pdb.\n";

    query << Ontology::END_WHERE_CLAUSE;

    FormulateCURL(output_file_type, query.str());
}
void Assembly::ExtractOntologyInfoByGlycanStructure(string ring_type, string anomeric_orientation, string minus_one_orientation, string index_two_orientation, string index_three_orientation,
                                                    string index_four_orientation, string plus_one_orientation, string output_file_type)
{
    if(ring_type.compare("") == 0)
    {
        cout << "Please specify the ring type which is the first argument of the function as either \"P\" or \"F\" " << endl;
        return;
    }
    stringstream query;
    query << Ontology::PREFIX << Ontology::SELECT_CLAUSE << " ?stereo_name ?stereo_condensed_name ?condensed_name ?name ?oligo_sequence ?residue_links ?glycosidic_linkage " << Ontology::WHERE_CLAUSE;
    query << "?pdb_file     :hasOligo       ?oligo.\n";
    query << "?oligo	    :hasCore        ?mono.\n";

    if(anomeric_orientation.compare("") != 0)
    {
        query << "?mono         :hasRingAtom	?anomeric.\n";
        query << "?anomeric	    :ringIndex  	\"1\".\n";
        query << "?anomeric	    :hasSideAtom    ?a_side.\n";
        query << "?a_side	    :sideIndex      \"1\".\n";
        query << "?a_side	    :orientation	\"" << anomeric_orientation << "\".\n";
    }
    if(anomeric_orientation.compare("") != 0 && minus_one_orientation.compare("") != 0)
    {
        query << "?anomeric     :hasSideAtom    ?a_minus_side.\n";
        query << "?a_minus_side	:sideIndex      \"-1\".\n";
        query << "?a_side	    :orientation	\"" << anomeric_orientation << "\".\n";
    }
    else if(minus_one_orientation.compare("") != 0)
    {
        query << "?mono         :hasRingAtom	?anomeric.\n";
        query << "?anomeric     :ringIndex  	\"1\".\n";
        query << "?anomeric     :hasSideAtom    ?a_minus_side.\n";
        query << "?a_minus_side	:sideIndex      \"-1\".\n";
        query << "?a_minus_side	:orientation	\"" << minus_one_orientation << "\".\n";
    }
    if(index_two_orientation.compare("") != 0)
    {
        query << "?mono         :hasRingAtom	?two.\n";
        query << "?two  	    :ringIndex  	\"2\".\n";
        query << "?two          :hasSideAtom    ?two_side.\n";
        query << "?two_side	    :sideIndex      \"2\".\n";
        query << "?two_side	    :orientation	\"" << index_two_orientation << "\".\n";
    }
    if(index_three_orientation.compare("") != 0)
    {
        query << "?mono         :hasRingAtom	?three.\n";
        query << "?three        :ringIndex  	\"3\".\n";
        query << "?three        :hasSideAtom    ?three_side.\n";
        query << "?three_side	:sideIndex      \"3\".\n";
        query << "?three_side	:orientation	\"" << index_three_orientation << "\".\n";
    }
    if(index_four_orientation.compare("") != 0)
    {
        query << "?mono         :hasRingAtom	?four.\n";
        query << "?four         :ringIndex  	\"4\".\n";
        query << "?four         :hasSideAtom    ?four_side.\n";
        query << "?four_side    :sideIndex      \"4\".\n";
        query << "?four_side	:orientation	\"" << index_four_orientation << "\".\n";
    }
    if(plus_one_orientation.compare("") != 0)
    {
        query << "?mono         :hasRingAtom	?last_c.\n";
        if(ring_type.compare("P") == 0 )
            query << "?last_c       :ringIndex  	\"5\".\n";
        else
            query << "?last_c       :ringIndex  	\"4\".\n";
        query << "?last_c        :hasSideAtom    ?plus_one.\n";
        query << "?plus_one      :sideIndex      \"+1\".\n";
        query << "?plus_one    	 :orientation	\"" << plus_one_orientation << "\".\n";
    }

    query << "?oligo        :oligoName  ?oligo_sequence.\n";
    query << "?pdb_file     :identifier	?pdb.\n";
    query << "OPTIONAL { ?oligo        :oligoResidueLinks	?residue_links.\n";
    query << "?linkage      :hasParent 	?oligo.\n";
    query << "?linkage      :glycosidicLinkage    ?glycosidic_linkage.}\n";
    query << "?mono         :hasSugarName	?sn.\n";
    query << "?sn           :monosaccharideName 	?name.\n";
    query << "?sn           :monosaccharideShortName 	?condensed_name.\n";
    query << "?sn           :monosaccharideStereochemName 	?stereo_name.\n";
    query << "?sn           :monosaccharideStereochemShortName 	?stereo_condensed_name.\n";
    query << Ontology::END_WHERE_CLAUSE;

    FormulateCURL(output_file_type, query.str());

}
void Assembly::ExtractOntologyInfoByDerivativeModificationMap(string ring_type, DerivativeModificationMap derivative_modification_map, string output_file_type)
{
    if(ring_type.compare("") == 0)
    {
        cout << "Please specify the ring type as the first argument of the function as either \"P\" or \"F\" " << endl;
        return;
    }
    stringstream query;
    query << Ontology::PREFIX << Ontology::SELECT_CLAUSE << "?pdb ?stereo_name ?stereo_condensed_name ?condensed_name " << Ontology::WHERE_CLAUSE;
    for(DerivativeModificationMap::iterator it = derivative_modification_map.begin(); it != derivative_modification_map.end(); it++)
    {
        string index = (*it).first;
        string pattern = (*it).second;
        query << "?mono	       :hasRingAtom    ?ring_atom.\n";
        if(index.compare("-1") != 0 && index.compare("+1") != 0)
        {
            query << "?ring_atom    :ringIndex      \"" << index << "\".\n";
            query << "?ring_atom    :hasSideAtom    ?side.\n";
        }
        else if(index.compare("-1") == 0 )
        {
            query << "?ring_atom       :ringIndex      \"1\".\n";
            query << "?ring_atom       :hasSideAtom    ?side.\n";
        }
        else if(index.compare("+1") == 0 )
        {
            query << "?mono         :hasRingAtom	?ring_atom.\n";
            if(ring_type.compare("P") == 0 )
                query << "?ring_atom       :ringIndex  	\"5\".\n";
            else
                query << "?ring_atom       :ringIndex  	\"4\".\n";
            query << "?ring_atom        :hasSideAtom    ?side.\n";
        }
        query << "?side         :sideIndex      \"" << index << "\".\n";
        query << "?side         :derivative     \"" << pattern << "\".\n";
        query << "?mono         :hasSugarName    ?sn.\n";
        query << "?sn           :monosaccharideName 	?name.\n";
        query << "?sn           :monosaccharideShortName 	?condensed_name.\n";
        query << "?sn           :monosaccharideStereochemName 	?stereo_name.\n";
        query << "?sn           :monosaccharideStereochemShortName 	?stereo_condensed_name.\n";
        query << "?pdb_file     :hasOligo 	?oligo.\n";
        query << "?oligo        :hasCore 	?mono.\n";
        query << "?pdb_file     :identifier	?pdb.\n";
    }
    query << Ontology::END_WHERE_CLAUSE;

    FormulateCURL(output_file_type, query.str());
}
void Assembly::ExtractOntologyInfoByAttachedGlycanStructures(AttachedGlycanStructuresVector attached_structures, string output_file_type)
{
    stringstream query;
    query << Ontology::PREFIX << Ontology::SELECT_CLAUSE << " ?pdb ?linkage_indices ?stereo_short_name0 ?short_name0 ?stereo_short_name1 ?short_name1 ?oligo_sequence ?residue_links "
          << Ontology::WHERE_CLAUSE;
    int i = 0;
    vector<string> oligos = vector<string>();
    for(AttachedGlycanStructuresVector::iterator it = attached_structures.begin(); it != attached_structures.end(); it++)
    {
        vector<string> structure = (*it);
        if(structure.size() < 7)
        {
            cout << "Missing arguments! All should be set even as an empty value (for empty values set \"\")" << endl;
            return;
        }
        if(structure.at(0).compare("") == 0)
        {
            cout << "Please specify the ring type which is the first argument of the function as either \"P\" or \"F\" " << endl;
            return;
        }
        stringstream oligo;
        oligo << "?oligo" << i;
        stringstream mono;
        mono << "?mono" << i;
        query << oligo.str() << "		:hasCore	" << mono.str() << ".\n";

        for(int j = 1; j < 7; j++)
        {
            if(structure.at(j).compare("") != 0)
            {
                stringstream ring_atom;
                stringstream side_atom;
                switch (j)///anomeric, -1 and +1 are special cases
                {
                    case 1:///anomeric
                        ring_atom << mono.str() << "_anomeric";
                        side_atom << mono.str() << "_anomeric_side_atom";
                        query << mono.str() << "     :hasRingAtom	" << ring_atom.str() << ".\n";
                        query << ring_atom.str() << "  	:ringIndex  	\"" << j << "\".\n";
                        query << ring_atom.str() << "   :hasSideAtom    " << side_atom.str() << " .\n";
                        query << side_atom.str() << "	:sideIndex      \"" << j << "\".\n";
                        query << side_atom.str() << "	:orientation	\"" << structure.at(j) << "\".\n";
                        break;
                    case 2:///minus one
                        ring_atom << mono.str() << "_anomeric";
                        side_atom << mono.str() << "_anomeric_minus_one_side_atom";
                        if(structure.at(1).compare("") == 0)///anomeric has not been set
                        {
                            query << mono.str() << "     :hasRingAtom	" << ring_atom.str() << ".\n";
                            query << ring_atom.str() << "  	:ringIndex  	\""<< j - 1 << "\".\n";
                        }
                        query << ring_atom.str() << "   :hasSideAtom    " << side_atom.str() << " .\n";
                        query << ring_atom.str() << "	:sideIndex     \"-1\".\n";
                        query << side_atom.str() << "	:orientation	\"" << structure.at(j) << "\".\n";
                        break;
                    case 6:///plus one
                        ring_atom << mono.str() << "_last_c";
                        side_atom << mono.str() << "_last_c_side_atom";
                        query << mono.str() << "         :hasRingAtom	" << ring_atom.str() << ".\n";
                        if(structure.at(0).compare("P") == 0 )
                            query << ring_atom.str() << "       :ringIndex  	\"" << j - 1 << "\".\n";
                        else
                            query << ring_atom.str() << "       :ringIndex  	\"" << j - 2 << "\".\n";
                        query << ring_atom.str() << "         :hasSideAtom    " << side_atom.str() << ".\n";
                        query << side_atom.str() << "       :sideIndex      \"+1\".\n";
                        query << side_atom.str() << "    	  :orientation	\"" << structure.at(j) << "\".\n";
                        break;
                    default:
                        ring_atom << mono.str() << "_ring_atom" << j - 1;
                        side_atom << mono.str() << "_side_atom" << j - 1;
                        query << mono.str() << "     :hasRingAtom	" << ring_atom.str() << ".\n";
                        query << ring_atom.str() << "  	:ringIndex  	\"" << j - 1 << "\".\n";
                        query << ring_atom.str() << "   :hasSideAtom    " << side_atom.str() << " .\n";
                        query << side_atom.str() << "	:sideIndex      \"" << j - 1 << "\".\n";
                        query << side_atom.str() << "	:orientation	\"" << structure.at(j) << "\".\n";
                }
            }
        }
        i++;
        oligos.push_back(oligo.str());
    }
    for(i = 0; i < oligos.size(); i++)
    {
        if(i + 1 < oligos.size())
        {
            query << "{\n";
            query << "?linkage" << i << " :hasParent " << oligos.at(i) << ".\n";
            query << "?linkage" << i << " :hasChild " << oligos.at(i + 1) << ".\n";
            query << oligos.at(i) << " :oligoResidueLinks ?residue_links.\n";
            query << "OPTIONAL {" << oligos.at(i) << "  :oligoName ?oligo_sequence}\n";
            query << "} UNION {\n";
            query << "?linkage" << i << " :hasParent " << oligos.at(i + 1) << ".\n";
            query << "?linkage" << i << " :hasChild " << oligos.at(i) << ".\n";
            query << oligos.at(i + 1) << " :oligoResidueLinks ?residue_links.\n";
            query << "OPTIONAL {" << oligos.at(i + 1) << "  :oligoName ?oligo_sequence}\n";
            query << "}\n";

            query << "?linkage" << i << ":linkageIndeces ?linkage_indices.\n";

            query << "?mono0    :hasSugarName ?sn0.\n";
            query << "?sn0      :monosaccharideShortName ?short_name0.\n";
            query << "?sn0      :monosaccharideStereochemShortName ?stereo_short_name0.\n";

            query << "?mono1    :hasSugarName ?sn1.\n";
            query << "?sn1      :monosaccharideShortName ?short_name1.\n";
            query << "?sn1      :monosaccharideStereochemShortName ?stereo_short_name1.\n";
        }
    }
    for(i = 0; i < oligos.size(); i++)
    {
        query << "?pdb            :hasOligo       " << oligos.at(i) << ".\n";
        ///optional {?oligo0    :oligoName ?oligoName0.}???
        ///optional {?oligo1    :oligoName ?oligoName1.}??? only one of these should be matched
    }
    query << Ontology::END_WHERE_CLAUSE;

    FormulateCURL(output_file_type, query.str());
}
void Assembly::ExtractOntologyInfoByNote(string pdb_id, string note_type, string note_category, string output_file_type)
{
    stringstream query;
    query << Ontology::PREFIX << Ontology::SELECT_CLAUSE << " ?pdb ?note_type ?note_category ?description "<< Ontology::WHERE_CLAUSE;

    if(pdb_id.compare("") != 0)
    {
        query <<  ":" << pdb_id << "    :hasNote   ?note.\n";
        query <<  ":" << pdb_id << "    :identifier    ?pdb.\n";
    }
    else
    {
        query << "?pdb_file      :hasNote    ?note.\n";
        query << "?pdb_file      :identifier    ?pdb.\n";
    }
    if(note_type.compare("") != 0)
        query << "?note	         :NoteType      \"" << note_type << "\".\n";
    query << "?note	       :NoteType    ?note_type.\n";

    if(note_category.compare("") != 0)
        query << "?note	       :NoteCategory      \"" << note_category << "\".\n";
    query << "?note	       :NoteCategory    ?note_category.\n";
    query << "?note	       :description    ?description.\n";

    query << Ontology::END_WHERE_CLAUSE;
    FormulateCURL(output_file_type, query.str());
}

void Assembly::ExtractOntologyInfoByCustomQuery(string query_file, string output_file_type)
{
    string line;
    stringstream query;
    ifstream in(query_file.c_str());

    if (!in.is_open())
    {
        cout << "Error in reading the query file" << endl;
        return;
    }

    while (getline (in, line))
    {
        query << line << endl;
    }
    in.close();

    FormulateCURL(output_file_type, query.str());
}

void Assembly::ExtractAtomCoordinatesForTorsionAnglesFromOntologySlow(string disaccharide_pattern, string output_file_type)
{

    int link_index = disaccharide_pattern.find_first_of("-");
    string child_mono = disaccharide_pattern.substr(0, link_index - 1); /// e.g DNeupNAca in DNeupNAca2-3DGalpb
    string parent_mono = disaccharide_pattern.substr(link_index + 2, disaccharide_pattern.size()); /// e.g DGalpb in DNeupNAca2-3DGalpb
    string linkage_indeces = disaccharide_pattern.substr(link_index - 1, 3); /// e.g 2-3 in DNeupNAca2-3DGalpb
    bool omega = false;
    if(linkage_indeces.find("6") != string::npos || linkage_indeces.find("7") != string::npos
             || linkage_indeces.find("8") != string::npos || linkage_indeces.find("9") != string::npos) /// Disaccharides with any of 1-6, 1-7, 2-6, and 2-7 linkages might have omega torsion angles
        omega = true;

    stringstream query;
    query << Ontology::PREFIX << Ontology::SELECT_CLAUSE << " ?pdb ?O5_crd ?C1_crd ?Ox_crd ?Cx ?Cx_crd ?Cx_neighbor ?Cx_neighbor_crd ?O5_prime_crd "<< Ontology::WHERE_CLAUSE;

    if(child_mono.compare("*") == 0)
        query <<  "?oligo1        :hasCore    ?mono1. \n";
    else
    {
        if(child_mono.find("*") != string::npos)
        {
            query <<  "?sn1           :monosaccharideShortName    ?short_name1. \n";
            query << "FILTER regex(?short_name1, \"" << Split(child_mono, "*").at(0) << "\", \"i\")\n";
        }
        else
            query <<  "?sn1           :monosaccharideShortName    \"" << child_mono << "\". \n";
        query <<  "?mono1         :hasSugarName    ?sn1.\n";
        query <<  "?oligo1        :hasCore    ?mono1. \n";
    }

    if(parent_mono.compare("*") == 0)
        query <<  "?oligo2        :hasCore    ?mono2. \n";
    else
    {
        if(parent_mono.find("*") != string::npos)
        {
            query <<  "?sn2           :monosaccharideShortName    ?short_name2. \n";
            query << "FILTER regex(?short_name2, \"" << Split(parent_mono, "*").at(0) << "\", \"i\")\n";
        }
        else
            query <<  "?sn2           :monosaccharideShortName    \"" << parent_mono << "\". \n";
        query <<  "?mono2         :hasSugarName    ?sn2. \n";
        query <<  "?oligo2        :hasCore    ?mono2. \n";
    }

    query <<  "?link          :hasParent   ?oligo2. \n";
    query <<  "?link          :hasChild    ?oligo1. \n";

    if(linkage_indeces.find("?") == string::npos)
        query <<  "?link          :linkageIndeces   \"" << linkage_indeces << "\". \n";
    else if(linkage_indeces.compare("?-?") != 0 && linkage_indeces.find("?") != string::npos)
    {
        query <<  "?link          :linkageIndeces   ?linkage_indeces. \n";
        query << "FILTER regex(?linkage_indeces, \"" << Split(linkage_indeces, "?").at(0) << "\", \"i\")\n";
    }

    query <<  "?pdb           :hasOligo    ?oligo1. \n";
    query <<  "?pdb           :hasOligo    ?oligo2. \n";

    query <<  "?mono1         :hasRingAtom    ?O5. \n";
    query <<  "?O5            :ringIndex    \"6\". \n";
    query <<  "?O5            :coordinate    ?O5_crd. \n";

    query <<  "?link          :hasChildAtomLinkage    ?C1. \n";
    query <<  "?C1            :coordinate    ?C1_crd. \n";

    query <<  "?link          :hasGlycosidicLinkage    ?Ox. \n";
    query <<  "?Ox            :coordinate    ?Ox_crd. \n";

    query <<  "?link          :hasParentAtomLinkage    ?Cx. \n";
    query <<  "?Cx            :coordinate     ?Cx_crd. \n";

    query <<  "?Cx            :hasNeighbor    ?Cx_neighbor. \n";
    query <<  "?Cx_neighbor   :coordinate     ?Cx_neighbor_crd. \n";

    if(omega)
    {
        query <<  "?mono2         :hasRingAtom    ?O5_prime.\n";
        query <<  "?O5_prime      :ringIndex    \"6\".\n";
        query <<  "?O5_prime      :coordinate    ?O5_prime_crd.\n";
    }

    query << Ontology::END_WHERE_CLAUSE;
    stringstream curl;
    curl << Ontology::CURL_PREFIX;
    curl << Ontology::CSV_OUTPUT_FORMAT;

    curl << Ontology::DATA_STORE_ADDRESS << Ontology::QUERY_PREFIX << query.str() << Ontology::QUERY_POSTFIX << " \>\> result.txt";
    string tmp = curl.str();
    const char* cstr = tmp.c_str();
    system(cstr);
}

void Assembly::ExtractAtomCoordinatesForTorsionAnglesFromOntologyFast(string disaccharide_pattern, string output_file_type)
{
    int link_index = disaccharide_pattern.find_first_of("-");

    string child_mono = disaccharide_pattern.substr(0, link_index - 1); /// e.g DNeupNAca in DNeupNAca2-3DGalpb
    string parent_mono = disaccharide_pattern.substr(link_index + 2, disaccharide_pattern.size()); /// e.g DGalpb in DNeupNAca2-3DGalpb
    string linkage_indeces = disaccharide_pattern.substr(link_index - 1, 3);
    bool omega = false;
    if(linkage_indeces.find("-6") != string::npos || linkage_indeces.find("-7") != string::npos)
        omega = true;

    stringstream query;
    query << "sparql PREFIX : <http://gmmo.uga.edu/#> " <<
             "PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#> " <<
             "PREFIX owl: <http://www.w3.org/2002/07/owl#> " <<
             "PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#> " <<
             "PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>  " <<
             "SELECT ?pdb ?O5_crd ?C1_crd ?Ox_crd ?Cx ?Cx_crd ?Cx_neighbor ?Cx_neighbor_crd ?O5_prime_crd WHERE { " ;

    query <<  "?sn1           :monosaccharideShortName    \"" << child_mono << "\". ";
    query <<  "?mono1         :hasSugarName    ?sn1. ";
    query <<  "?oligo1        :hasCore    ?mono1. ";

    query <<  "?sn2           :monosaccharideShortName    \"" << parent_mono << "\". ";
    query <<  "?mono2         :hasSugarName    ?sn2. ";
    query <<  "?oligo2        :hasCore    ?mono2. ";

    query <<  "?link          :hasParent   ?oligo2. ";
    query <<  "?link          :hasChild    ?oligo1. ";
    query <<  "?link          :linkageIndeces   \"" << linkage_indeces << "\". ";

    query <<  "?pdb           :hasOligo    ?oligo1. ";
    query <<  "?pdb           :hasOligo    ?oligo2. ";

    query <<  "?mono1         :hasRingAtom    ?O5. ";
    query <<  "?O5            :ringIndex    \"6\". ";
    query <<  "?O5            :coordinate    ?O5_crd. ";

    query <<  "?link          :hasChildAtomLinkage    ?C1. ";
    query <<  "?C1            :coordinate    ?C1_crd. ";

    query <<  "?link          :hasGlycosidicLinkage    ?Ox. ";
    query <<  "?Ox            :coordinate    ?Ox_crd. ";

    query <<  "?link          :hasParentAtomLinkage    ?Cx. ";
    query <<  "?Cx            :coordinate     ?Cx_crd. ";

    query <<  "?Cx            :hasNeighbor    ?Cx_neighbor. ";
    query <<  "?Cx_neighbor   :coordinate     ?Cx_neighbor_crd. ";

    if(omega)
    {
        query <<  "?mono2         :hasRingAtom    ?O5_prime. ";
        query <<  "?O5_prime      :ringIndex    \"6\". ";
        query <<  "?O5_prime      :coordinate    ?O5_prime_crd. ";
    }

    query << "};";

    std::ofstream sparql;
    sparql.open("sparql.sparql", fstream::app);
    sparql << query.str() ;
    sparql.close();

//    stringstream ss;
//    ss << "/home/delaram/virtuoso-7.2.4/bin/isql 1111 dba dba \< sparql.sparql \>  result.txt";
//    cout << ss.str() << endl;
//    string tmp = ss.str();
//    const char* cstr = tmp.c_str();
    system("/home/delaram/virtuoso-7.2.4/bin/isql 1111 dba dba \< sparql.sparql \>  result.txt");
    remove("sparql.sparql");

//    stringstream curl;
//    curl << Ontology::CURL_PREFIX;
//    curl << Ontology::CSV_OUTPUT_FORMAT;

//    curl << Ontology::DATA_STORE_ADDRESS << Ontology::QUERY_PREFIX << query.str() << Ontology::QUERY_POSTFIX << " >> result.txt";
//    string tmp = curl.str();
//    const char* cstr = tmp.c_str();
//    system(cstr);
//    cout << endl;



 /*   int link_index = disaccharide_pattern.find_first_of("-");

    string child_mono = disaccharide_pattern.substr(0, link_index - 1); /// e.g DNeupNAca in DNeupNAca2-3DGalpb
    string parent_mono = disaccharide_pattern.substr(link_index + 2, disaccharide_pattern.size()); /// e.g DGalpb in DNeupNAca2-3DGalpb
    string linkage_indeces = disaccharide_pattern.substr(link_index - 1, 3);
    bool omega = false;
    if(linkage_indeces.find("6") != string::npos || linkage_indeces.find("7") != string::npos
             || linkage_indeces.find("8") != string::npos || linkage_indeces.find("9") != string::npos)
        omega = true;

    stringstream query;
    query << Ontology::PREFIX << Ontology::SELECT_CLAUSE << " ?pdb ?RA1 ?RA1_crd ?SA1 ?SA1_crd ?RA2 ?RA2_crd ?SA2 ?SA2_crd "<< Ontology::WHERE_CLAUSE;
    query << "?sn1           :monosaccharideShortName    \"" << child_mono << "\".\n";
    query << "?mono1         :hasSugarName    ?sn1.\n";
    query << "?oligo1        :hasCore    ?mono1.\n";
    query << "?sn2           :monosaccharideShortName    \"" << parent_mono << "\".\n";
    query << "?mono2         :hasSugarName    ?sn2.\n";
    query << "?oligo2        :hasCore    ?mono2.\n";
    query << "?link          :hasParent   ?oligo2.\n";
    query << "?link          :hasChild    ?oligo1.\n";
    query << "?link          :linkageIndeces   \"" << linkage_indeces << "\".\n";
    query << "?pdb           :hasOligo    ?oligo1.\n";
    query << "?pdb           :hasOligo    ?oligo2.\n";

    query << "?mono1      :hasRingAtom  ?RA1.\n";
    query << "?RA1        :coordinate    ?RA1_crd.\n";

    query << "optional {?RA1        :hasSideAtom    ?SA1.\n";
    query << "?SA1        :coordinate     ?SA1_crd.}\n";


    query << "?mono2      :hasRingAtom  ?RA2.\n";
    query << "?RA2        :coordinate    ?RA2_crd.\n";

    query << "optional {?RA2        :hasSideAtom    ?SA2.\n";
    query << "?SA2        :coordinate     ?SA2_crd.}\n";


    query << Ontology::END_WHERE_CLAUSE;
    FormulateCURL(output_file_type, query.str());
    */
}

void Assembly::ExtractTorsionAnglesFromSlowQueryResult()
{

    ///Uncomment the following section and substitute cout with out_file in order to write the results into a file
    /*
    ///Open file to append
    std::ofstream out_file;
    out_file.open("torsions.txt", fstream::app);

    ///If the file is empty, write the titles
    ifstream check_file("torsions.txt");
    size_t out_file_size = 0;
    check_file.seekg(0,ios_base::end);
    out_file_size = check_file.tellg();
    check_file.close();
    if(out_file_size == 0)
    {
        out_file << "ϕ (O5′-C1′-Ox-Cx)" << endl;
        out_file << "ψ (C1′-Ox-Cx-Cx−1)" << endl;
        out_file << "Ω (O1-C6′-C5′-O5′)" << endl;
        out_file << left << setw(15) << "PDB" << setw(15) << "Phi Angle" << setw(15) << "Psi Angle" << setw(15) << "Omega Angle" << endl;
    } */


    ///Outputting the result in std out. In order to write the results into a file comment the following section
    cout << "ϕ (O5′-C1′-Ox-Cx)" << endl;
    cout << "ψ (C1′-Ox-Cx-Cx−1)" << endl;
    cout << "Ω (O1-C6′-C5′-O5′)" << endl;
    cout << left << setw(15) << "PDB" << setw(15) << "Phi Angle" << setw(15) << "Psi Angle" << setw(15) << "Omega Angle" << endl;

    ///Read query result file
    string line;
    ifstream in("result.txt");

    while (getline (in, line))
    {
        if(line.find("\"pdb\",\"O5_crd\",\"C1_crd\",\"Ox_crd\",\"Cx\",\"Cx_crd\",\"Cx_neighbor\",\"Cx_neighbor_crd\",\"O5_prime_crd\"") != string::npos)
            break;
    }
    string last_Cx = "";

    Coordinate* O5_crd = new Coordinate();
    Coordinate* C1_crd = new Coordinate();
    Coordinate* Ox_crd = new Coordinate();
    Coordinate* Cx_crd = new Coordinate();
    Coordinate* Cx_neighbor_crd = new Coordinate();
    Coordinate* O5_prime_crd = NULL;
    string Cx_id = "";
    string Cx_neighbor_id = "";
    double phi_angle = 0.0;
    double psi_angle = 0.0;
    double omega_angle = 0.0;

    while (getline (in, line))
    {
        vector<string> line_tokens = Split(line, ",");

        if(last_Cx.compare("") == 0 || last_Cx.compare(line_tokens.at(10)) != 0)//If it is the first line or the line with info about a new oligosaccharide
        {
            ///e.g. "0.686, -15.194, 26.371" --after splits--> x=0.686 y=-15.194 z=26.371
            O5_crd->SetX(ConvertString<double>(Split(line_tokens.at(1), "\"").at(0)));
            O5_crd->SetY(ConvertString<double>(Split(line_tokens.at(2), " ").at(0)));
            O5_crd->SetZ(ConvertString<double>(Split(line_tokens.at(3), " \"").at(0)));

            C1_crd->SetX(ConvertString<double>(Split(line_tokens.at(4), "\"").at(0)));
            C1_crd->SetY(ConvertString<double>(Split(line_tokens.at(5), " ").at(0)));
            C1_crd->SetZ(ConvertString<double>(Split(line_tokens.at(6), " \"").at(0)));

            Ox_crd->SetX(ConvertString<double>(Split(line_tokens.at(7), "\"").at(0)));
            Ox_crd->SetY(ConvertString<double>(Split(line_tokens.at(8), " ").at(0)));
            Ox_crd->SetZ(ConvertString<double>(Split(line_tokens.at(9), " \"").at(0)));

            Cx_crd->SetX(ConvertString<double>(Split(line_tokens.at(11), "\"").at(0)));
            Cx_crd->SetY(ConvertString<double>(Split(line_tokens.at(12), " ").at(0)));
            Cx_crd->SetZ(ConvertString<double>(Split(line_tokens.at(13), " \"").at(0)));

            if(line_tokens.size() > 18)
            {
                O5_prime_crd = new Coordinate();
                O5_prime_crd->SetX(ConvertString<double>(Split(line_tokens.at(18), "\"").at(0)));
                O5_prime_crd->SetY(ConvertString<double>(Split(line_tokens.at(19), " ").at(0)));
                O5_prime_crd->SetZ(ConvertString<double>(Split(line_tokens.at(20), " \"").at(0)));
            }

            Cx_id = Split(line_tokens.at(10), "#").at(1); ///e.g. spliting 5BO9_C3_4773_GAL_A_410_n_n_1 from http://gmmo.uga.edu/#5BO9_C3_4773_GAL_A_410_n_n_1
        }

        last_Cx = line_tokens.at(10);
        Cx_neighbor_id = Split(line_tokens.at(14), "#").at(1);

        if(Split(Cx_neighbor_id, "_").at(1).find("C") != string::npos)///If the neighbor is a carbon
        {
            /// e.g. removing C * , and ' from the atom name to get the index
            int Cx_index = ConvertString<int>(Split(Split(Cx_id, "_").at(1), "C*,\'").at(0)); ///e.g. 5BO9_C3_4773_GAL_A_410_n_n_1 -split-> C3 -split-> 3
            int Cx_neighbor_index = ConvertString<int>(Split(Split(Cx_neighbor_id, "_").at(1), "C*,\'").at(0));

            if(Cx_index > Cx_neighbor_index) ///if the neighbor is Cx-1
            {
                Cx_neighbor_crd->SetX(ConvertString<double>(Split(line_tokens.at(15), "\"").at(0)));
                Cx_neighbor_crd->SetY(ConvertString<double>(Split(line_tokens.at(16), " ").at(0)));
                Cx_neighbor_crd->SetZ(ConvertString<double>(Split(line_tokens.at(17), "\"").at(0)));

                phi_angle = CalculateTorsionAngleByCoordinates(O5_crd, C1_crd, Ox_crd, Cx_crd); /// ϕ (O5′-C1′-Ox-Cx)
                psi_angle = CalculateTorsionAngleByCoordinates(C1_crd, Ox_crd, Cx_crd,Cx_neighbor_crd); /// ψ (C1′-Ox-Cx-Cx−1)

                cout << left << setw(15) << Split(Split(line_tokens.at(0), "#").at(1), "\"").at(0)
                         << setw(15) << ConvertRadian2Degree(phi_angle) << setw(15)
                         << ConvertRadian2Degree(psi_angle);

                if(O5_prime_crd != NULL)
                {
                    omega_angle = CalculateTorsionAngleByCoordinates(Ox_crd, Cx_crd,Cx_neighbor_crd, O5_prime_crd); /// Ω (O1-C6′-C5′-O5′)
                    cout << setw(15) << ConvertRadian2Degree(omega_angle);
                }
                cout << endl;
            }
        }
    }
    in.close();
}

void Assembly::ExtractTorsionAnglesFromFastQueryResult()
{

    /* Sample query result

     */

    ///Uncomment the following section and substitute cout with out_file in order to write the results into a file
    /*
    ///Open file to append
    std::ofstream out_file;
    out_file.open("torsions.txt", fstream::app);

    ///If the file is empty, write the titles
    ifstream check_file("torsions.txt");
    size_t out_file_size = 0;
    check_file.seekg(0,ios_base::end);
    out_file_size = check_file.tellg();
    check_file.close();
    if(out_file_size == 0)
    {
        out_file << "ϕ (O5′-C1′-Ox-Cx)" << endl;
        out_file << "ψ (C1′-Ox-Cx-Cx−1)" << endl;
        out_file << "Ω (O1-C6′-C5′-O5′)" << endl;
        out_file << left << setw(15) << "PDB" << setw(15) << "Phi Angle" << setw(15) << "Psi Angle" << setw(15) << "Omega Angle" << endl;
    } */


    ///Outputting the result in std out. In order to write the results into a file comment the following section
    cout << "ϕ (O5′-C1′-Ox-Cx)" << endl;
    cout << "ψ (C1′-Ox-Cx-Cx−1)" << endl;
    cout << "Ω (O1-C6′-C5′-O5′)" << endl;
    cout << left << setw(15) << "PDB" << setw(15) << "Phi Angle" << setw(15) << "Psi Angle" << setw(15) << "Omega Angle" << endl;

    ///Read query result file
    string line;
    ifstream in("result.txt");

    while (getline (in, line))
    {
        if(line.find("http://gmmo.uga.edu/") != string::npos)
            break;
    }
    string last_Cx = "";

    Coordinate* O5_crd = new Coordinate();
    Coordinate* C1_crd = new Coordinate();
    Coordinate* Ox_crd = new Coordinate();
    Coordinate* Cx_crd = new Coordinate();
    Coordinate* Cx_neighbor_crd = new Coordinate();
    Coordinate* O5_prime_crd = NULL;
    string Cx_id = "";
    string Cx_neighbor_id = "";
    double phi_angle = 0.0;
    double psi_angle = 0.0;
    double omega_angle = 0.0;

    do
    {
        vector<string> line_tokens = Split(line, " ");
        if(last_Cx.compare("") == 0 || last_Cx.compare(line_tokens.at(10)) != 0)//If it is the first line or the line with info about a new oligosaccharide
        {
            ///e.g. "0.686, -15.194, 26.371" --after splits--> x=0.686 y=-15.194 z=26.371

            O5_crd->SetX(ConvertString<double>(Split(line_tokens.at(1), ",").at(0)));
            O5_crd->SetY(ConvertString<double>(Split(line_tokens.at(2), ",").at(0)));
            O5_crd->SetZ(ConvertString<double>(line_tokens.at(3)));

            C1_crd->SetX(ConvertString<double>(Split(line_tokens.at(4), ",").at(0)));
            C1_crd->SetY(ConvertString<double>(Split(line_tokens.at(5), ",").at(0)));
            C1_crd->SetZ(ConvertString<double>(line_tokens.at(6)));

            Ox_crd->SetX(ConvertString<double>(Split(line_tokens.at(7), ",").at(0)));
            Ox_crd->SetY(ConvertString<double>(Split(line_tokens.at(8), ",").at(0)));
            Ox_crd->SetZ(ConvertString<double>(line_tokens.at(9)));

            Cx_crd->SetX(ConvertString<double>(Split(line_tokens.at(11), ",").at(0)));
            Cx_crd->SetY(ConvertString<double>(Split(line_tokens.at(12), ",").at(0)));
            Cx_crd->SetZ(ConvertString<double>(line_tokens.at(13)));

            if(line_tokens.at(18).compare("") != 0)
            {
                O5_prime_crd->SetX(ConvertString<double>(Split(line_tokens.at(18), ",").at(0)));
                O5_prime_crd->SetY(ConvertString<double>(Split(line_tokens.at(19), ",").at(0)));
                O5_prime_crd->SetZ(ConvertString<double>(line_tokens.at(20)));
            }

            Cx_id = Split(line_tokens.at(10), "#").at(1); ///e.g. spliting 5BO9_C3_4773_GAL_A_410_n_n_1 from http://gmmo.uga.edu/#5BO9_C3_4773_GAL_A_410_n_n_1
        }

        last_Cx = line_tokens.at(10);
        Cx_neighbor_id = Split(line_tokens.at(14), "#").at(1);

        if(Split(Cx_neighbor_id, "_").at(1).find("C") != string::npos)///If the neighbor is a carbon
        {
            /// e.g. removing C * , and ' from the atom name to get the index
            int Cx_index = ConvertString<int>(Split(Split(Cx_id, "_").at(1), "C*,\'").at(0)); ///e.g. 5BO9_C3_4773_GAL_A_410_n_n_1 -split-> C3 -split-> 3
            int Cx_neighbor_index = ConvertString<int>(Split(Split(Cx_neighbor_id, "_").at(1), "C*,\'").at(0));

            if(Cx_index > Cx_neighbor_index) ///if the neighbor is Cx-1
            {
                Cx_neighbor_crd->SetX(ConvertString<double>(Split(line_tokens.at(15), ",").at(0)));
                Cx_neighbor_crd->SetY(ConvertString<double>(Split(line_tokens.at(16), ",").at(0)));
                Cx_neighbor_crd->SetZ(ConvertString<double>(line_tokens.at(17)));

                phi_angle = CalculateTorsionAngleByCoordinates(O5_crd, C1_crd, Ox_crd, Cx_crd); /// ϕ (O5′-C1′-Ox-Cx)
                psi_angle = CalculateTorsionAngleByCoordinates(C1_crd, Ox_crd, Cx_crd,Cx_neighbor_crd); /// ψ (C1′-Ox-Cx-Cx−1)

                cout << left << setw(15) << Split(Split(line_tokens.at(0), "#").at(1), "\"").at(0)
                         << setw(15) << ConvertRadian2Degree(phi_angle) << setw(15)
                         << ConvertRadian2Degree(psi_angle);

                if(O5_prime_crd != NULL)
                {
                    omega_angle = CalculateTorsionAngleByCoordinates(Ox_crd, Cx_crd,Cx_neighbor_crd, O5_prime_crd); /// Ω (O1-C6′-C5′-O5′)
                    cout << setw(15) << ConvertRadian2Degree(omega_angle);
                }
                cout << endl;
            }
        }
    }
    while (getline (in, line) && !line.empty());
    in.close();
}

vector<double> Assembly::CalculateBondlengthsStatisticsBasedOnOntologyInfo(string atom_name1, string atom_name2, bool is_atom2_ring, string mono_name)
{
    stringstream query;
    query << "sparql PREFIX : <http://gmmo.uga.edu/#> " <<
             "PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#> " <<
             "PREFIX owl: <http://www.w3.org/2002/07/owl#> " <<
             "PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#> " <<
             "PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>" <<
             "SELECT ?atom1_crd ?atom2_crd WHERE {" ;

    query <<  "?mono          :hasSugarName   ?sn.";
    query <<  "?sn            :monosaccharideShortName   \"" << mono_name << "\".\n";
    query <<  "?mono          :hasRingAtom   ?atom1.";
    if(is_atom2_ring)
    {
        query <<  "?mono          :hasRingAtom   ?atom2.";
        query <<  "?atom1         :hasNeighbor    ?atom2.";
    }
    else
        query <<  "?atom1          :hasSideAtom   ?atom2.";
    query <<  "?atom1         :identifier    ?atom1_id.";
    query <<  "?atom2         :identifier    ?atom2_id.";
    query << "FILTER regex(?atom1_id, \"" << atom_name1 << "_" << "\", \"i\")";
    query << "FILTER regex(?atom2_id, \"" << atom_name2 << "_" << "\", \"i\")";

    query <<  "?atom1         :coordinate    ?atom1_crd.";
    query <<  "?atom2         :coordinate    ?atom2_crd.";
    query << "};";

    std::ofstream sparql;
    sparql.open("bonds.sparql", fstream::app);
    sparql << query.str();
    sparql.close();

    system("/home/delaram/virtuoso-7.2.4/bin/isql 1111 dba dba \< bonds.sparql \>  bond_results.txt");
    remove("bonds.sparql");

    ///Read query result file
    string line;
    ifstream in("bond_results.txt");
    while (getline (in, line))///skip the first lines until the coordinates
    {
        if(line.find("____") != string::npos)
        {
            getline (in, line);
            break;
        }
    }
    Coordinate* atom1_crd = new Coordinate();
    Coordinate* atom2_crd = new Coordinate();
    double distance = 0.0;
    double sum_of_bond_lengths = 0.0;
    double sum_of_bond_lengths_squared = 0.0;
    int number_of_bond_lengths = 0;

    ///Reading the coordinates from the result file, calculating the distance
    while (getline (in, line) && !line.empty())
    {
        vector<string> line_tokens = Split(line," ");
        atom1_crd->SetX(ConvertString<double>(Split(line_tokens.at(0), ",").at(0)));
        atom1_crd->SetY(ConvertString<double>(Split(line_tokens.at(1), ",").at(0)));
        atom1_crd->SetZ(ConvertString<double>(line_tokens.at(2)));

        atom2_crd->SetX(ConvertString<double>(Split(line_tokens.at(3), ",").at(0)));
        atom2_crd->SetY(ConvertString<double>(Split(line_tokens.at(4), ",").at(0)));
        atom2_crd->SetZ(ConvertString<double>(line_tokens.at(5)));

        distance = atom1_crd->Distance(*(atom2_crd));
        sum_of_bond_lengths += distance;
        sum_of_bond_lengths_squared += (distance*distance);
        number_of_bond_lengths++;
    }
    in.close();
    remove("bond_results.txt");

    double mean  = sum_of_bond_lengths/number_of_bond_lengths;
    double standard_deviation = sqrt(((sum_of_bond_lengths_squared - ((sum_of_bond_lengths*sum_of_bond_lengths)/number_of_bond_lengths))/number_of_bond_lengths));

    vector<double> statistics = vector<double>();
    statistics.push_back(mean);
    statistics.push_back(standard_deviation);

    return statistics;
}
vector<double> Assembly::CalculateBondAnglesStatisticsBasedOnOntologyInfo(string atom_name1, string atom_name2, string atom_name3, bool is_atom3_ring, string mono_name)
{
    stringstream query;
    query << "sparql PREFIX : <http://gmmo.uga.edu/#> " <<
             "PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#> " <<
             "PREFIX owl: <http://www.w3.org/2002/07/owl#> " <<
             "PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#> " <<
             "PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>" <<
             "SELECT ?atom1_crd ?atom2_crd ?atom3_crd WHERE {" ;

    query <<  "?mono          :hasSugarName   ?sn.";
    query <<  "?sn            :monosaccharideShortName   \"" << mono_name << "\".\n";
    query <<  "?mono          :hasRingAtom   ?atom1.";
    query <<  "?mono          :hasRingAtom   ?atom2.";
    query <<  "?atom1         :hasNeighbor    ?atom2.";
    if(is_atom3_ring)
    {
        query <<  "?mono          :hasRingAtom   ?atom3.";
        query <<  "?atom2         :hasNeighbor    ?atom3.";
        query <<  "FILTER (?atom1 != ?atom3)";
    }
    else
        query <<  "?atom2          :hasSideAtom   ?atom3.";
    query <<  "?atom1         :identifier    ?atom1_id.";
    query <<  "?atom2         :identifier    ?atom2_id.";
    query <<  "?atom3         :identifier    ?atom3_id.";
    query << "FILTER regex(?atom1_id, \"" << atom_name1 << "_" << "\", \"i\")";
    query << "FILTER regex(?atom2_id, \"" << atom_name2 << "_" << "\", \"i\")";
    query << "FILTER regex(?atom3_id, \"" << atom_name3 << "_" << "\", \"i\")";

    query <<  "?atom1         :coordinate    ?atom1_crd.";
    query <<  "?atom2         :coordinate    ?atom2_crd.";
    query <<  "?atom3         :coordinate    ?atom3_crd.";
    query << "};";

    std::ofstream sparql;
    sparql.open("bond_angles.sparql", fstream::app);
    sparql << query.str();
    sparql.close();

    system("/home/delaram/virtuoso-7.2.4/bin/isql 1111 dba dba \< bond_angles.sparql \>  bond_angle_results.txt");
    remove("bond_angles.sparql");

    ///Read query result file
    string line;
    ifstream in("bond_angle_results.txt");
    while (getline (in, line))///skip the first lines until the coordinates
    {
        if(line.find("____") != string::npos)
        {
            getline (in, line);
            break;
        }
    }
    Coordinate* atom1_crd = new Coordinate();
    Coordinate* atom2_crd = new Coordinate();
    Coordinate* atom3_crd = new Coordinate();
    double bond_angle = 0.0;
    double sum_of_bond_angles = 0.0;
    double sum_of_bond_angles_squared = 0.0;
    int number_of_bond_angles = 0;

    ///Reading the coordinates from the result file, calculating the distance
    while (getline (in, line) && !line.empty())
    {
        vector<string> line_tokens = Split(line," ");
        atom1_crd->SetX(ConvertString<double>(Split(line_tokens.at(0), ",").at(0)));
        atom1_crd->SetY(ConvertString<double>(Split(line_tokens.at(1), ",").at(0)));
        atom1_crd->SetZ(ConvertString<double>(line_tokens.at(2)));

        atom2_crd->SetX(ConvertString<double>(Split(line_tokens.at(3), ",").at(0)));
        atom2_crd->SetY(ConvertString<double>(Split(line_tokens.at(4), ",").at(0)));
        atom2_crd->SetZ(ConvertString<double>(line_tokens.at(5)));

        atom3_crd->SetX(ConvertString<double>(Split(line_tokens.at(6), ",").at(0)));
        atom3_crd->SetY(ConvertString<double>(Split(line_tokens.at(7), ",").at(0)));
        atom3_crd->SetZ(ConvertString<double>(line_tokens.at(8)));

        bond_angle = CalculateBondAngleByCoordinates(atom1_crd, atom2_crd, atom3_crd);
        sum_of_bond_angles += bond_angle;
        sum_of_bond_angles_squared += (bond_angle*bond_angle);
        number_of_bond_angles++;
    }

    in.close();
    remove("bond_angle_results.txt");

    double mean  = sum_of_bond_angles/number_of_bond_angles;
    double standard_deviation = sqrt(((sum_of_bond_angles_squared - ((sum_of_bond_angles*sum_of_bond_angles)/number_of_bond_angles))/number_of_bond_angles));

    vector<double> statistics = vector<double>();
    statistics.push_back(mean);
    statistics.push_back(standard_deviation);

    return statistics;
}


void Assembly::ExtractTorsionAnglesFromPDB(vector<string> amino_lib_files, string disaccharide)
{
    string pdb_file_path = this->GetSourceFile();
    string pdb = pdb_file_path.substr(pdb_file_path.find_last_of("/") + 1, pdb_file_path.size());
    OligosaccharideVector oligos = this->ExtractSugars(amino_lib_files);
    double phi_angle = 0.0;
    double psi_angle = 0.0;
    int link_index = disaccharide.find_first_of("-");
    char mono1_carbon_index = disaccharide.at(link_index - 1);
    char mono2_carbon_index = disaccharide.at(link_index + 1);
    string first_mono = disaccharide.substr(0, link_index - 1);
    string second_mono = disaccharide.substr(link_index + 2, disaccharide.size());

    std::ofstream out_file;
    out_file.open("torsions.txt", fstream::app);

    ifstream in("torsions.txt");///Checking if the file is empty
    size_t out_file_size = 0;
    in.seekg(0,ios_base::end);
    out_file_size = in.tellg();
    in.close();

    if(out_file_size == 0)
        out_file << left << setw(15) << "PDB" << setw(15) << "Phi Angle" << setw(15) << "Psi Angle" << setw(25) << "Pattern" << "Oligosaccharide" << endl;

    for(OligosaccharideVector::iterator it = oligos.begin(); it != oligos.end(); it++)
    {
        Oligosaccharide* oligo = (*it);
        string oligo_name = oligo->oligosaccharide_name_;
        ///The root of the sequence, last monosaccharide in the oligo name. sequence created backward from root to first mono
        if(oligo_name.compare("") != 0 && oligo_name.find(disaccharide) != string::npos)///found root and the oligo contains the disaccharide.
        {
            queue<Oligosaccharide*> oligo_queue;
            oligo_queue.push(oligo);
            if(MatchDisaccharide(oligo_queue, phi_angle, psi_angle, first_mono, mono1_carbon_index, second_mono, mono2_carbon_index))
                out_file << left << setw(15) << pdb << setw(15) << ConvertRadian2Degree(phi_angle) << setw(15) << ConvertRadian2Degree(psi_angle) << setw(25) << disaccharide << oligo_name << endl;
        }
    }
    out_file.close();
}

bool Assembly::MatchDisaccharide(queue<Oligosaccharide*> oligo_queue, double &phi_angle, double &psi_angle, string first_mono, char mono1_carbon_index, string second_mono, char mono2_carbon_index)
{
    Oligosaccharide* oligo = oligo_queue.front();
    oligo_queue.pop();
    bool found_disaccharide = false;

    Oligosaccharide* corresponding_second_oligo = NULL;
    ///second mono of the input disaccharide comes first in the tree structure of oligo (parent of first mono)
    if(oligo->root_->sugar_name_.monosaccharide_short_name_.compare(second_mono) == 0 )///found a mono with the same name as the right side mono of the disaccharide
        corresponding_second_oligo = oligo;

    OligosaccharideVector child_oligos = oligo->child_oligos_;
    for(OligosaccharideVector::iterator it1 = child_oligos.begin(); it1 != child_oligos.end(); it1++)
    {
        oligo_queue.push(*it1);

        if(corresponding_second_oligo != NULL) ///if current mono matches the disaccharide, look for a linked mono that macthes the other mono in the disaccharide
        {
            Oligosaccharide* corresponding_first_oligo = (*it1);
            oligo_queue.push(corresponding_first_oligo);

            if(corresponding_first_oligo->root_->sugar_name_.monosaccharide_short_name_.compare(first_mono) == 0)///found a mono with the same name as the left side mono of the disaccharide
            {
                vector<string> child_links = corresponding_second_oligo->child_oligos_linkages_;///links from right mono of the disaccharide to child monos

                vector<string> mono2_cycle_atom_tokens = Split(corresponding_second_oligo->root_->cycle_atoms_str_, "-");
                for(vector<string>::iterator it2 = child_links.begin(); it2 != child_links.end(); it2++)
                {
                    //                int index = distance(child_links.begin(), it2);
                    string link = (*it2);
                    vector<string> link_tokens = Split(link, "-");
                    int parent_c_index = 0;
                    int child_c_index = 0;

                    if(corresponding_second_oligo->root_->side_atoms_.at(0).at(0) != NULL)
                        parent_c_index++;
                    for(int i = 0; i < mono2_cycle_atom_tokens.size(); i++)
                    {
                        parent_c_index++;
                        if(mono2_cycle_atom_tokens.at(i).compare(link_tokens.at(0)) == 0)
                            break;
                    }
                    vector<string> mono1_cycle_atom_tokens = Split(corresponding_first_oligo->root_->cycle_atoms_str_, "-");
                    if(corresponding_first_oligo->root_->side_atoms_.at(0).at(0) != NULL)
                        child_c_index++;
                    for(int i = 0; i < mono1_cycle_atom_tokens.size(); i++)
                    {
                        child_c_index++;
                        if(mono1_cycle_atom_tokens.at(i).compare(link_tokens.at(2)) == 0)
                            break;
                    }

                    if(mono2_carbon_index == parent_c_index + '0' && mono1_carbon_index == child_c_index + '0')///indeces matched the indeces of the given disaccharide
                    {
                        found_disaccharide = true;
                        AtomVector mono1_ring_atoms = corresponding_first_oligo->root_->cycle_atoms_;

                        ///Preparing atoms for phi and psi angle
                        Atom* phi_atom1 = mono1_ring_atoms.at(mono1_ring_atoms.size() - 1); ///O5 or N5
                        Atom* phi_atom2 = new Atom(); ///C1
                        phi_atom2 = corresponding_first_oligo->root_->cycle_atoms_.at(0); ///anomeric carbon
                        Atom* phi_atom3 = NULL;
                        AtomVector phi_atom2_neighbors = phi_atom2->GetNode()->GetNodeNeighbors();
                        for(AtomVector::iterator it3 = phi_atom2_neighbors.begin(); it3 != phi_atom2_neighbors.end(); it3++)
                        {
                            Atom* atom2_neighbor= (*it3);
                            ///If the neighbor id is the same as the linkage intermediate atom id
                            if(atom2_neighbor->GetId().compare(link_tokens.at(1)) == 0)
                            {
                                phi_atom3 = atom2_neighbor; ///Ox
                                break;
                            }
                        }
                        if(phi_atom3 != NULL)
                        {
                            Atom* phi_atom4 = NULL;
                            AtomVector phi_atom3_neighbors = phi_atom3->GetNode()->GetNodeNeighbors();
                            for(AtomVector::iterator it4 = phi_atom3_neighbors.begin(); it4 != phi_atom3_neighbors.end(); it4++)
                            {
                                Atom* atom3_neighbor= (*it4);
                                ///If the neighbor id is the same as the linkage carbon at second mono side (first linkage atom from second mono side)
                                if(atom3_neighbor->GetId().compare(link_tokens.at(0)) == 0)
                                {
                                    phi_atom4 = atom3_neighbor; ///Cx
                                    phi_angle = CalculateTorsionAngleByAtoms(phi_atom1, phi_atom2, phi_atom3, phi_atom4); /// ϕ (O5′-C1′-Ox-Cx)
                                    AtomVector phi_atom4_neighbors = phi_atom4->GetNode()->GetNodeNeighbors();

                                    int atom4_index = ConvertString<int>(Split(phi_atom4->GetName(), "C*,\'").at(0));
                                    for(AtomVector::iterator it5 = phi_atom4_neighbors.begin(); it5 != phi_atom4_neighbors.end(); it5++)
                                    {
                                        Atom* atom4_neighbor = (*it5);
                                        string neighbor_name = atom4_neighbor->GetName();
                                        if(neighbor_name.find("C") != string::npos)
                                        {
                                            int neighbor_index = ConvertString<int>(Split(neighbor_name, "C*,\'").at(0));
                                            if(atom4_index > neighbor_index) ///Cx-1
                                            {
                                                psi_angle = CalculateTorsionAngleByAtoms(phi_atom2, phi_atom3, phi_atom4, atom4_neighbor); /// ψ (C1′-Ox-Cx-Cx−1)
                                                break;
                                            }
                                        }
                                    }
                                    break;
                                }
                            }
                        }
                    }
                    if(found_disaccharide)
                        break;
                }
                if(found_disaccharide)
                    break;
            }
        }
    }
    if(!found_disaccharide && !oligo_queue.empty())
        return MatchDisaccharide(oligo_queue, phi_angle, psi_angle, first_mono, mono1_carbon_index, second_mono, mono2_carbon_index);
    else
        return found_disaccharide;
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

void Assembly::CalculateTorsionStatistics(string torsion_file, int low_range, int high_range)
{
    /* Sample output of the function
   |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  1 |  0 |  0 |  0 |  0 |  0 |  0 |  0 | 180
   |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |
   |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  1 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |
   |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  3 |  0 |  1 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |
   |  0 |  0 |  0 |  0 |  0 |  0 |  1 |  0 |  0 |  1 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |
   |  0 |  0 |  0 |  0 |  0 |  1 |  1 |  0 |  1 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |
   |  0 |  0 |  0 |  0 |  0 |  1 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |
   |  0 |  0 |  0 |  0 |  1 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |
Psi|  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |
   |  0 |  0 |  0 |  1 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |
   |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |
   |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |
   |  2 |  0 |  0 |  0 |  0 |  1 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |
   |  0 |  0 |  0 |  0 |  0 |  1 |  1 |  0 |  0 |  5 |  1 |  0 |  0 |  1 |  0 |  1 |  0 |  0 |
   |  0 |  0 |  1 |  0 | 15 | 24 | 15 |  0 |  0 |  0 |  0 |  2 |  3 |  5 |  0 |  0 |  1 |  0 |
   |  0 |  0 |  0 |  1 | 18 | 41 | 43 |  0 |  0 |  0 |  1 | 14 |  6 |  0 |  0 |  0 |  1 |  0 |
   |  0 |  0 |  1 |  0 |  1 | 23 |  9 |  0 |  1 |  2 |  4 | 31 |  1 |  1 |  0 |  0 |  0 |  0 |
   |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  2 |  0 |  5 |  1 |  0 |  0 |  0 |  0 |  0 |  0 | -180
-180                                         Phi                                          180
     */
    vector<vector<int> > psi_phi_matrix(18, vector<int> (18, 0)); ///matrix 18x18
    vector<vector<int> > psi_omega_matrix(18, vector<int> (18, 0)); ///matrix 18x18
    vector<vector<int> > phi_omega_matrix(18, vector<int> (18, 0)); ///matrix 18x18
    int scale = (high_range - low_range) / 18; ///range of each cell in the matrix
    int x = 0;
    int y = 0;
    bool has_omega = false;

    string line;
    string input_file = "";
    if(torsion_file.compare("") == 0)
        input_file = "torsions.txt";
    else
        input_file = torsion_file;

    ifstream in(input_file.c_str());
    if (!in.is_open())
        cout << "Error in reading the torsion results file" << endl;

    ///Read the head of the file untill the data about torsions
    while (getline (in, line) )
    {
        if(line.find("PDB") != string::npos)
            break;
    }

    while (getline (in, line) && !line.empty())
    {
        x = 0;
        y = 0;
        vector<string> line_tokens = Split(line, " ");
        x = (abs(low_range) + ConvertString<double>(line_tokens.at(1))) / scale; ///x axis in psi_phi_matrix
        y = (abs(low_range) + ConvertString<double>(line_tokens.at(2))) / scale; ///y axis in psi_phi_matrix
        psi_phi_matrix.at(y).at(x) += 1;

        if(line_tokens.size() > 3 && line_tokens.at(3).compare("") != 0) ///x axis in psi_omega_matrix
        {
            has_omega = true;

            x = (abs(low_range) + ConvertString<double>(line_tokens.at(3))) / scale; ///x axis in psi_omega_matrix and phi_omega_matrix
            psi_omega_matrix.at(y).at(x) += 1;

            y = (abs(low_range) + ConvertString<double>(line_tokens.at(1))) / scale; ///y axis in phi_omega_matrix
            phi_omega_matrix.at(y).at(x) += 1;
        }
    }
    in.close();

    for(int j = 17; j >= 0; j-- )
    {
        if(j == 9)
            cout << "Psi| ";
        else
            cout << setw(5) << right << "| ";
        for(int k = 0; k <= 17; k++ )
        {
            cout << setw(2) << psi_phi_matrix.at(j).at(k) << " | ";
        }
        if(j == 17)
            cout << high_range;
        if(j == 0)
            cout << low_range;
        cout << endl;
    }
    cout << setw(5) << left << low_range << setw(43) << right << "Phi" << setw(45) << right << high_range << endl << endl;

    if(has_omega)
    {
        for(int j = 17; j >= 0; j-- )
        {
            if(j == 9)
                cout << "Psi| ";
            else
                cout << setw(5) << right << "| ";
            for(int k = 0; k <= 17; k++ )
            {
                cout << setw(2) << psi_omega_matrix.at(j).at(k) << " | ";
            }
            if(j == 17)
                cout << high_range;
            if(j == 0)
                cout << low_range;
            cout << endl;
        }
        cout << setw(5) << left << low_range << setw(43) << right << "Omega" << setw(45) << right << high_range << endl << endl;

        for(int j = 17; j >= 0; j-- )
        {
            if(j == 9)
                cout << "Phi| ";
            else
                cout << setw(5) << right << "| ";
            for(int k = 0; k <= 17; k++ )
            {
                cout << setw(2) << phi_omega_matrix.at(j).at(k) << " | ";
            }
            if(j == 17)
                cout << high_range;
            if(j == 0)
                cout << low_range;
            cout << endl;
        }
        cout << setw(5) << left << low_range << setw(43) << right << "Omega" << setw(45) << right << high_range << endl;
    }
}

void Assembly::ExtractRingAtomsInformation()
{
    ///CYCLE DETECTION
    CycleMap cycles = DetectCyclesByExhaustiveRingPerception();
    ///FILTERING OUT FUSED CYCLES
    RemoveFusedCycles(cycles);
    ///FILTERING OUT OXYGENLESS CYCLES
    FilterAllCarbonCycles(cycles);
    CycleMap sorted_cycles = CycleMap();
    vector<string> anomeric_carbons_status = vector<string>();
    ///ANOMERIC DETECTION and SORTING
    for(CycleMap::iterator it = cycles.begin(); it != cycles.end(); it++)
    {
        string cycle_atoms_str = (*it).first;
        AtomVector cycle_atoms = (*it).second;
        Note* anomeric_note = new Note();
        Atom* anomeric = FindAnomericCarbon(anomeric_note, anomeric_carbons_status, cycle_atoms, cycle_atoms_str);

        if(anomeric != NULL)
        {
            AtomVector sorted_cycle_atoms = AtomVector();
            stringstream sorted_cycle_stream;
            sorted_cycle_atoms = SortCycle(cycle_atoms, anomeric, sorted_cycle_stream);
            sorted_cycles[sorted_cycle_stream.str()] = sorted_cycle_atoms;
        }
    }
    cycles = sorted_cycles;
    cout << endl << "Detailed information of sorted cycles after discarding fused or oxygenless rings: " << endl;
    gmml::log(__LINE__, __FILE__,  gmml::INF, "Detailed information of sorted cycles after discarding fused or oxygenless rings: ");

    for(CycleMap::iterator it = cycles.begin(); it != cycles.end(); it++)
    {
        Monosaccharide* mono = new Monosaccharide();
        int status_index = distance(cycles.begin(), it);
        mono->anomeric_status_ = anomeric_carbons_status.at(status_index);
        string cycle_atoms_str = (*it).first;
        AtomVector cycle = (*it).second;

        stringstream ring_atoms;
        ring_atoms << "Ring atoms: " << cycle_atoms_str;
        cout << ring_atoms.str() << endl;
        gmml::log(__LINE__, __FILE__,  gmml::INF, ring_atoms.str());

        stringstream ring_atom_names;
        ring_atom_names << "Ring atom names: " ;
        for(AtomVector::iterator it1 = cycle.begin(); it1 != cycle.end(); it1++)
        {
            Atom* atom = (*it1);
            vector<string> atom_tokens = Split(atom->GetId(), "_");
            if(it1 == cycle.end() - 1)
                ring_atom_names << atom_tokens.at(0);
            else
                ring_atom_names << atom_tokens.at(0) << ",";
        }
        cout << ring_atom_names.str() << endl;
        gmml::log(__LINE__, __FILE__,  gmml::INF, ring_atom_names.str());
    }
}

Assembly::CycleMap Assembly::DetectCyclesByExhaustiveRingPerception()
{
    CycleMap cycles = CycleMap();
    AtomVector atoms = GetAllAtomsOfAssemblyExceptProteinWaterResiduesAtoms();
    vector<string> path_graph_edges = vector<string> (); ///The list of edges in the molecular graph
    vector<string> path_graph_labels = vector<string> (); ///The list of labels of edges in the molecular graph
    vector<string> cycless = vector<string>();

    ///Initializing the map
    map<string, Atom*> IdAtom = map<string, Atom*>(); ///A map from atom ID to Assembly atom object
    for(AtomVector::iterator it = atoms.begin(); it != atoms.end(); it++)
    {
        Atom* atom = (*it);
        IdAtom[atom->GetId()] = atom;
    }
    ///Pruning the graph (filter out atoms with less than 2 neighbors)
    PruneGraph(atoms);

    ///Converting the molecular graph into a path graph
    ConvertIntoPathGraph(path_graph_edges, path_graph_labels, atoms);
    vector<string> reduced_path_graph_edges = path_graph_edges;
    vector<string> reduced_path_graph_labels = path_graph_labels;

    int neighbor_counter = 2;
    ///Reducing the path graph
    ///Whenever a walk a-b-c is found it should be reduced to a-c and the label should be changed from [a-b], [b-c] to [a-b-c]
    /// the node with lowest number of edges to other nodes should be examined first
    while(atoms.size() > 1 && path_graph_edges.size() != 0)
    {
        AtomVector::iterator common_atom_it;
        bool neighbor_counter_update = true;
        for(AtomVector::iterator it = atoms.begin(); it != atoms.end(); it++)
        {
            common_atom_it = it;
            int counter = 0;
            for(vector<string>::iterator it1 = path_graph_edges.begin(); it1 != path_graph_edges.end(); it1++)
            {
                string edge = (*it1);
                if(edge.find((*common_atom_it)->GetId()) != string::npos)
                    counter++;
            }
            if(counter <= neighbor_counter)
            {
                neighbor_counter_update = false;
                break;
            }
        }
        if(neighbor_counter_update)
        {
            neighbor_counter++;
            continue;
        }

        ReducePathGraph(path_graph_edges, path_graph_labels, reduced_path_graph_edges, reduced_path_graph_labels, (*common_atom_it)->GetId(), cycless);

        atoms.erase(common_atom_it);

        path_graph_edges = reduced_path_graph_edges;
        path_graph_labels = reduced_path_graph_labels;
    }

    for(vector<string>::iterator it = cycless.begin(); it != cycless.end(); it++)
    {
        string cycle = (*it);
        vector<string> splitted_cycle = Split(cycle, "-");
        if(splitted_cycle.size() <= 7)
        {
            AtomVector atomvector = AtomVector();
            stringstream ss;
            for(vector<string>::iterator it1 = splitted_cycle.begin(); it1 != splitted_cycle.end() - 1; it1++)
            {
                string atom_str = (*it1);
                map<string, Atom*>::iterator mit = IdAtom.find(atom_str);
                atomvector.push_back((*mit).second);
                if(it1 == splitted_cycle.end() - 2)
                    ss << atom_str;
                else
                    ss << atom_str << "-";
            }
            cycles[ss.str()] = atomvector;
        }
    }
    return cycles;

    /*
    CycleMap cycles = CycleMap();
    AtomVector atoms = GetAllAtomsOfAssemblyExceptProteinWaterResiduesAtoms();
    vector<string> path_graph_edges = vector<string> (); ///The list of edges in the molecular graph
    vector<string> path_graph_labels = vector<string> (); ///The list of labels of edges in the molecular graph
    vector<string> cycless = vector<string>();

    ///Initializing the map
    map<string, Atom*> IdAtom = map<string, Atom*>(); ///A map from atom ID to Assembly aTom object
    for(AtomVector::iterator it = atoms.begin(); it != atoms.end(); it++)
    {
        Atom* atom = (*it);
        IdAtom[atom->GetId()] = atom;
    }

    ///Pruning the graph (filter out atoms with less than 2 neighbors)
    PruneGraph(atoms);

    ///Converting the molecular graph into a path graph
    ConvertIntoPathGraph(path_graph_edges, path_graph_labels, atoms);
    vector<string> reduced_path_graph_edges = path_graph_edges;
    vector<string> reduced_path_graph_labels = path_graph_labels;

    ///Reducing the path graph
    int neighbor_counter = 2;
    ///Whenever a walk a-b-c is found it should be reduced to a-c and the lable should be changed from [a-b], [b-c] to [a-b-c]
    /// the node with lowest number of connected edges should be examined first
    while(atoms.size() > 1 && path_graph_edges.size() != 0)
    {
        AtomVector::iterator common_atom_it;
        bool neighbor_counter_update = true;
        for(AtomVector::iterator it = atoms.begin(); it != atoms.end(); it++)
        {
            common_atom_it = it;
            int counter = 0;
            for(vector<string>::iterator it1 = path_graph_edges.begin(); it1 != path_graph_edges.end(); it1++)
            {
                string edge = (*it1);
                if(edge.find((*common_atom_it)->GetId()) != string::npos) ///finding edges connected to the node
                    counter++;
            }
            if(counter <= neighbor_counter)///A node with lower number of edges has been found
            {
                neighbor_counter_update = false;
                break;
            }
        }
        if(neighbor_counter_update)
        {
            neighbor_counter++;
            continue;
        }

        ReducePathGraph(path_graph_edges, path_graph_labels, reduced_path_graph_edges, reduced_path_graph_labels, (*common_atom_it)->GetId(), cycless);

        atoms.erase(common_atom_it);

        path_graph_edges = reduced_path_graph_edges;
    }

    for(vector<string>::iterator it = cycless.begin(); it != cycless.end(); it++)
    {
        string cycle = (*it);
        vector<string> splitted_cycle = Split(cycle, "-");
        if(splitted_cycle.size() <= 7)
        {
            AtomVector atomvector = AtomVector();
            stringstream ss;
            for(vector<string>::iterator it1 = splitted_cycle.begin(); it1 != splitted_cycle.end() - 1; it1++)
            {
                string atom_str = (*it1);
                map<string, Atom*>::iterator mit = IdAtom.find(atom_str);
                atomvector.push_back((*mit).second);
                if(it1 == splitted_cycle.end() - 2)
                    ss << atom_str;
                else
                    ss << atom_str << "-";
            }
            cycles[ss.str()] = atomvector;
        }
    }
    return cycles; */
}

void Assembly::ReducePathGraph(vector<string> path_graph_edges, vector<string> path_graph_labels, vector<string>& reduced_path_graph_edges,
                               vector<string>& reduced_path_graph_labels, string common_atom, vector<string>& cycles)
{
    vector<int> to_be_deleted_edges = vector<int>();
    for(vector<string>::iterator it = path_graph_edges.begin(); it != path_graph_edges.end() - 1; it++)
    {
        int source_index = distance(path_graph_edges.begin(), it);
        string source_edge = (*it);
        string source_label = path_graph_labels.at(source_index);
        vector<string> source_edge_atoms = Split(source_edge, ",");

        for(vector<string>::iterator it1 = it + 1; it1 != path_graph_edges.end(); it1++)
        {
            bool walk_found = false;
            int target_index = distance(path_graph_edges.begin(), it1);
            string target_edge = (*it1);
            string target_label = path_graph_labels.at(target_index);
            vector<string> target_edge_atoms = Split(target_edge, ",");
            stringstream new_edge;
            stringstream new_label;
            if(source_edge_atoms.at(0).compare(source_edge_atoms.at(1)) != 0 && target_edge_atoms.at(0).compare(target_edge_atoms.at(1)) != 0)
            {///if the edges a != b and b != c, so there might be a walk a_b_c

                if(source_edge_atoms.at(1).compare(target_edge_atoms.at(0)) == 0 && source_edge_atoms.at(1).compare(common_atom) == 0)///if there is a walk a_b_c in the graph (edges: a,b and b,c)
                {
                    new_edge << source_edge_atoms.at(0) << "," << target_edge_atoms.at(1);
                    new_label << source_label;
                    vector<string> target_path_values = Split(target_label,"-");
                    for(int i = 1; i < target_path_values.size(); i++)
                        new_label << "-" << target_path_values.at(i);
                    walk_found = true;
                }
                else if(source_edge_atoms.at(1).compare(target_edge_atoms.at(1)) == 0 && source_edge_atoms.at(1).compare(common_atom) == 0)///if there is a walk a_b_c in the graph (edges: a,b and c,b)
                {
                    new_edge << source_edge_atoms.at(0) << "," << target_edge_atoms.at(0);
                    new_label << source_label;
                    vector<string> target_path_values = Split(target_label,"-");
                    for(int i = target_path_values.size() - 2 ; i >= 0; i--)
                        new_label << "-" << target_path_values.at(i);
                    walk_found = true;
                }
                else if(source_edge_atoms.at(0).compare(target_edge_atoms.at(0)) == 0 && source_edge_atoms.at(0).compare(common_atom) == 0)///if there is a walk a_b_c in the graph (edges: b,a and b,c)
                {
                    new_edge << source_edge_atoms.at(1) << "," << target_edge_atoms.at(1);
                    vector<string> source_path_values = Split(source_label,"-");
                    for(int i = source_path_values.size() - 1 ; i >= 1; i--)
                        new_label << source_path_values.at(i) << "-";
                    new_label << target_label;
                    walk_found = true;
                }
                else if(source_edge_atoms.at(0).compare(target_edge_atoms.at(1)) == 0 && source_edge_atoms.at(0).compare(common_atom) == 0)///if there is a walk a_b_c in the graph (edges: b,a and c,b)
                {
                    new_edge << source_edge_atoms.at(1) << "," << target_edge_atoms.at(0);
                    vector<string> source_path_values = Split(source_label,"-");
                    for(int i = source_path_values.size() - 1 ; i >= 0; i--)
                        new_label << source_path_values.at(i) << "-";
                    vector<string> target_path_values = Split(target_label,"-");
                    for(int i = target_path_values.size() - 2 ; i >= 0; i--)
                    {
                        if(i == 0)
                            new_label << target_path_values.at(i);
                        else
                            new_label << target_path_values.at(i) << "-";
                    }
                    walk_found = true;
                }
            }
            if(walk_found)
            {
                ///checking the new edge for cycle
                vector<string> new_edge_atoms = Split(new_edge.str(), ",");
                if(new_edge_atoms.at(0).compare(new_edge_atoms.at(1)) == 0) ///edge is a,a
                {
                    cycles.push_back(new_label.str());///label shows the atom involved in a cycle
                }
                ///adding the newly-formed edge (a,c) and label(a-b-c)
                else if(find(reduced_path_graph_labels.begin(), reduced_path_graph_labels.end(), new_label.str()) == reduced_path_graph_labels.end())
                {
                    reduced_path_graph_edges.push_back(new_edge.str());
                    reduced_path_graph_labels.push_back(new_label.str());
                }

                ///adding to be deleted edges with the common atom b
                if(find(to_be_deleted_edges.begin(), to_be_deleted_edges.end(), source_index) == to_be_deleted_edges.end())
                    to_be_deleted_edges.push_back(source_index);
                if(find(to_be_deleted_edges.begin(), to_be_deleted_edges.end(), target_index) == to_be_deleted_edges.end())
                    to_be_deleted_edges.push_back(target_index);
            }
        }
    }

    vector<string> temp_reduced_path_graph_edges = vector<string>();
    vector<string> temp_reduced_path_graph_labels = vector<string>();
    for(int i = 0; i < reduced_path_graph_edges.size(); i++)
    {
        if(find(to_be_deleted_edges.begin(), to_be_deleted_edges.end(), i) == to_be_deleted_edges.end())
        {
            temp_reduced_path_graph_edges.push_back(reduced_path_graph_edges.at(i));
            temp_reduced_path_graph_labels.push_back(reduced_path_graph_labels.at(i));
        }
    }
    reduced_path_graph_edges = temp_reduced_path_graph_edges;
    reduced_path_graph_labels = temp_reduced_path_graph_labels;
}

void Assembly::PruneGraph(AtomVector& all_atoms)
{
    AtomVector atoms_with_more_than_two_neighbors = AtomVector();
    vector<string> het_atom_ids = vector<string>();
    for(AtomVector::iterator it = all_atoms.begin(); it != all_atoms.end(); it++)
    {
        Atom* atom = *it;
        het_atom_ids.push_back(atom->GetId());
    }
    for(AtomVector::iterator it = all_atoms.begin(); it != all_atoms.end(); it++)
    {
        Atom* atom = (*it);
        AtomNode* node = atom->GetNode();
        int count = 0;
        for(int i = 0; i < node->GetNodeNeighbors().size(); i++)
        {
            Atom* neighbor = node->GetNodeNeighbors().at(i);
            if(find(het_atom_ids.begin(), het_atom_ids.end(), neighbor->GetId()) != het_atom_ids.end())
                count++;
        }
        if(count > 1)
            atoms_with_more_than_two_neighbors.push_back(atom);
    }
    if(atoms_with_more_than_two_neighbors.size() != all_atoms.size())
    {
        all_atoms = atoms_with_more_than_two_neighbors;
        PruneGraph(all_atoms);
    }
    else
        return;
}

void Assembly::ConvertIntoPathGraph(vector<string>& path_graph_edges, vector<string>& path_graph_labels, AtomVector atoms)
{
    vector<string> atoms_id = vector<string>();
    for(AtomVector::iterator it = atoms.begin(); it != atoms.end(); it++)
    {
        Atom* atom = *it;
        atoms_id.push_back(atom->GetId());
    }
    for(AtomVector::iterator it = atoms.begin(); it != atoms.end(); it++)
    {
        Atom* atom = (*it);
        AtomNode* node = atom->GetNode();
        AtomVector neighbors = node->GetNodeNeighbors();
        for(AtomVector::iterator it1 = neighbors.begin(); it1 != neighbors.end(); it1++)
        {
            Atom* neighbor = (*it1);
            if(find(atoms_id.begin(), atoms_id.end(), neighbor->GetId()) != atoms_id.end())
            {
                stringstream ss;
                ss << atom->GetId() << "," << neighbor->GetId();
                stringstream reverse_ss;
                reverse_ss << neighbor->GetId() << "," << atom->GetId();
                if(find(path_graph_edges.begin(), path_graph_edges.end(), ss.str()) == path_graph_edges.end() &&
                        find(path_graph_edges.begin(), path_graph_edges.end(), reverse_ss.str()) == path_graph_edges.end()) ///path not existed before
                {
                    stringstream path;
                    path << atom->GetId() << "-" << neighbor->GetId();
                    path_graph_edges.push_back(ss.str());
                    path_graph_labels.push_back(path.str());
                }
            }
        }
    }
}

Assembly::CycleMap Assembly::DetectCyclesByDFS()
{
    int counter = 0;

    AtomStatusMap atom_status_map = AtomStatusMap();
    AtomIdAtomMap atom_parent_map = AtomIdAtomMap();
    AtomIdAtomMap src_dest_map = AtomIdAtomMap();
    AtomVector cycle = AtomVector();
    CycleMap cycles = CycleMap();

    AtomVector atoms = GetAllAtomsOfAssemblyExceptProteinWaterResiduesAtoms();
    for(AtomVector::iterator it = atoms.begin(); it != atoms.end(); it++)
    {
        Atom* atom = (*it);
        atom_status_map[atom->GetId()] = gmml::UNVISITED;
        Atom* parent = new Atom();
        parent->SetId("null");
        atom_parent_map[atom->GetId()] = parent;
    }
    for(AtomVector::iterator it = atoms.begin(); it != atoms.end(); it++)
    {
        Atom* atom = (*it);
        if(atom_status_map[atom->GetId()] == gmml::UNVISITED)
        {
            DFSVisit(atoms, atom_status_map, atom_parent_map, atom, counter, src_dest_map);
        }
    }

    stringstream n_of_cycle;
    n_of_cycle << "Number of cycles found: " << counter;
    cout << n_of_cycle.str() << endl;
    gmml::log(__LINE__, __FILE__,  gmml::INF, n_of_cycle.str());
    for(AtomIdAtomMap::iterator it = src_dest_map.begin(); it != src_dest_map.end(); it++)
    {
        string src_dest = (*it).first;
        Atom* destination = (*it).second;
        cycle.clear();
        stringstream cycle_stream;
        vector<string> key = Split(src_dest, "-");
        ReturnCycleAtoms(key.at(0), destination, atom_parent_map, cycle, cycle_stream);
        cycles[cycle_stream.str()] = cycle;
    }
    return cycles;
}

void Assembly::DFSVisit(AtomVector atoms, AtomStatusMap& atom_status_map, AtomIdAtomMap& atom_parent_map, Atom *atom, int& counter, AtomIdAtomMap& src_dest_map)
{
    atom_status_map[atom->GetId()] = gmml::VISITED;
    AtomNode* node = atom->GetNode();
    AtomVector neighbors = node->GetNodeNeighbors();

    for(AtomVector::iterator it = neighbors.begin(); it != neighbors.end(); it++)
    {
        Atom* neighbor = (*it);
        if(neighbor->GetDescription().find("Het;") != string::npos)
        {
            if(atom_status_map[neighbor->GetId()] == gmml::UNVISITED)
            {
                atom_parent_map[neighbor->GetId()] = atom;
                DFSVisit(atoms, atom_status_map, atom_parent_map, neighbor, counter, src_dest_map);
            }
            if(atom_status_map[neighbor->GetId()] == gmml::VISITED)
            {
                Atom* parent = atom_parent_map[atom->GetId()];
                if(neighbor->GetId().compare(parent->GetId()) != 0)///making sure we are not tracking back to the previous atom which is the parent of neigbor (current atom)
                {
                    counter++;
                    stringstream key;
                    key << neighbor->GetId() << "-" << atom->GetId();
                    src_dest_map[key.str()] = atom;
                }
            }
        }
    }
    atom_status_map[atom->GetId()] = gmml::DONE;
}

void Assembly::ReturnCycleAtoms(string src_id, Atom *current_atom, AtomIdAtomMap &atom_parent_map, AtomVector &cycle, stringstream &cycle_stream)
{
    cycle.push_back(current_atom);
    cycle_stream << current_atom->GetId() << "-";
    Atom* parent = atom_parent_map[current_atom->GetId()];
    if(src_id.compare(parent->GetId()) == 0)
    {
        cycle.push_back(parent);
        cycle_stream << parent->GetId();
        return;
    }
    ReturnCycleAtoms(src_id, parent, atom_parent_map, cycle, cycle_stream);
}

void Assembly::FilterAllCarbonCycles(CycleMap &cycles)
{
    map<string, bool> to_be_deleted_cycles = map<string, bool>();
    bool all_carbons = true;
    for(CycleMap::iterator it = cycles.begin(); it != cycles.end(); it++)
    {
        string cycle_str = (*it).first;
        AtomVector cycle_atoms = (*it).second;
        all_carbons = true;
        for(AtomVector::iterator it1 = cycle_atoms.begin(); it1 != cycle_atoms.end(); it1++)
        {
            Atom* atom = (*it1);
            atom->GetName().find("C") == string::npos;
            {
                all_carbons = false;
                break;
            }
        }
        if(all_carbons)
            to_be_deleted_cycles[cycle_str] = true;
    }
    CycleMap all_carbons_filtered_cycles = CycleMap();
    for(CycleMap::iterator it = cycles.begin(); it != cycles.end(); it++)
    {
        string cycle_str = (*it).first;
        AtomVector cycle_atoms = (*it).second;
        if(to_be_deleted_cycles.find(cycle_str) == to_be_deleted_cycles.end())
            all_carbons_filtered_cycles[cycle_str] = cycle_atoms;
    }
    cycles.clear();
    cycles = all_carbons_filtered_cycles;
}

void Assembly::RemoveFusedCycles(CycleMap &cycles)
{
    map<string, bool> to_be_deleted_cycles = map<string, bool>();
    for(CycleMap::iterator it = cycles.begin(); it != cycles.end(); it++)
    {
        string cycle_i_str = (*it).first; ///cycle i will be compared with all the other cycles
        AtomVector cycle_i_atoms = (*it).second;
        for(CycleMap::iterator it1 = cycles.begin(); it1 != cycles.end(); it1++)
        {
            if(it != it1)
            {
                string cycle_j_str = (*it1).first; ///cycle j to be compared with cycle i
                for(int i = 0; i < cycle_i_atoms.size(); i++)
                {
                    stringstream mutual_edge;
                    stringstream mutual_edge_reverse;
                    Atom* a1 = new Atom();
                    Atom* a2 = new Atom();
                    if(i == cycle_i_atoms.size() - 1)
                    {
                        a1 = cycle_i_atoms.at(i);
                        a2 = cycle_i_atoms.at(0);
                    }
                    else
                    {
                        a1 = cycle_i_atoms.at(i);
                        a2 = cycle_i_atoms.at(i + 1);
                    }
                    mutual_edge << a1->GetId() << "-" << a2->GetId();
                    mutual_edge_reverse << a2->GetId() << "-" << a1->GetId();
                    if(cycle_j_str.find(mutual_edge.str()) != string::npos || cycle_j_str.find(mutual_edge_reverse.str()) != string::npos)///mutual edge found
                    {
                        to_be_deleted_cycles[cycle_i_str] = true;
                        to_be_deleted_cycles[cycle_j_str] = true;
                        break;
                    }
                }
            }
        }
    }
    CycleMap fused_filtered_cycles = CycleMap();
    for(CycleMap::iterator it = cycles.begin(); it != cycles.end(); it++)
    {
        string cycle_str = (*it).first;
        AtomVector cycle_atoms = (*it).second;
        if(to_be_deleted_cycles.find(cycle_str) == to_be_deleted_cycles.end())
            fused_filtered_cycles[cycle_str] = cycle_atoms;
    }
    cycles.clear();
    cycles = fused_filtered_cycles;
}

Atom* Assembly::FindAnomericCarbon(Note* anomeric_note, vector<string>& anomeric_carbons_status, AtomVector cycle, string cycle_atoms_str)
{
    Atom* anomeric_carbon = new Atom();
    for(AtomVector::iterator it = cycle.begin(); it != cycle.end(); it++)
    {
        Atom* cycle_atom = (*it);
        if((cycle_atom->GetName().substr(0,1).compare("O") == 0 ))///find oxygen in ring
            //                && isdigit(ConvertString<char>(cycle_atom->GetName().substr(1,1)))))
        {
            AtomNode* node = cycle_atom->GetNode();
            AtomVector neighbors = node->GetNodeNeighbors();

            ///Check the first neighbor of oxygen
            Atom* o_neighbor1 = neighbors.at(0);
            AtomNode* o_neighbor1_node = o_neighbor1->GetNode();
            AtomVector o_neighbor1_neighbors = o_neighbor1_node->GetNodeNeighbors();
            for(AtomVector::iterator it1 = o_neighbor1_neighbors.begin(); it1 != o_neighbor1_neighbors.end(); it1++)///check if neighbor1 of oxygen has another oxygen or nitrogen neighbor
            {
                Atom* neighbor1_neighbor = (*it1);
                if(cycle_atoms_str.find(neighbor1_neighbor->GetId()) == string::npos ///if the neighbor is not one of the cycle atoms
                        && (neighbor1_neighbor->GetName().substr(0,1).compare("O") == 0 || neighbor1_neighbor->GetName().substr(0,1).compare("N") == 0)) ///if first element is "O" or "N"
                    //                        && isdigit(ConvertString<char>(neighbor1_neighbor->GetName().substr(1,1))))///if second element is a digit
                {
                    anomeric_carbon = o_neighbor1;
                    anomeric_carbons_status.push_back("Anomeric carbon: ");
                    anomeric_note->description_ = "";

                    return anomeric_carbon;
                }
            }

            ///Check the second neighbor of oxygen
            Atom* o_neighbor2 = neighbors.at(1);
            AtomNode* o_neighbor2_node = o_neighbor2->GetNode();
            AtomVector o_neighbor2_neighbors = o_neighbor2_node->GetNodeNeighbors();
            for(AtomVector::iterator it2 = o_neighbor2_neighbors.begin(); it2 != o_neighbor2_neighbors.end(); it2++)///check if neighbor2 of oxygen has another oxygen or nitrogen neighbor
            {
                Atom* neighbor2_neighbor = (*it2);
                if( cycle_atoms_str.find(neighbor2_neighbor->GetId()) == string::npos
                        && (neighbor2_neighbor->GetName().substr(0,1).compare("O") == 0 || neighbor2_neighbor->GetName().substr(0,1).compare("N") == 0))
                    //                        && isdigit(ConvertString<char>(neighbor2_neighbor->GetName().substr(1,1))))
                {
                    anomeric_carbon = o_neighbor2;
                    //                    cout << "Anomeric carbon is: " << anomeric_carbon->GetId() << endl;
                    anomeric_carbons_status.push_back("Anomeric carbon: ");
                    anomeric_note->description_ = "";
                    return anomeric_carbon;
                }
            }

            ///Check the order of the carbons based on their names to locate the anomeric
            stringstream ss1;
            for(int i = 0; i < o_neighbor1->GetName().size(); i++)
            {
                if(isdigit(o_neighbor1->GetName().at(i)) != 0)
                    ss1 << o_neighbor1->GetName().at(i);
            }
            stringstream ss2;
            for(int i = 0; i < o_neighbor2->GetName().size(); i++)
            {
                if(isdigit(o_neighbor2->GetName().at(i)) != 0)
                    ss2 << o_neighbor2->GetName().at(i);
            }
            if(ConvertString<int>(ss1.str()) < ConvertString<int>(ss2.str()))
            {
                anomeric_note->type_ = Glycan::WARNING;
                anomeric_note->category_ = Glycan::ANOMERIC;
                anomeric_note->description_ = "Anomeric oxygen is missing";
                anomeric_carbons_status.push_back("Anomeric carbon probably is: ");
                return o_neighbor1;
            }
            if(ConvertString<int>(ss2.str()) < ConvertString<int>(ss1.str()))
            {
                anomeric_note->type_ = Glycan::WARNING;
                anomeric_note->category_ = Glycan::ANOMERIC;
                anomeric_note->description_ = "Anomeric oxygen is missing";
                anomeric_carbons_status.push_back("Anomeric carbon probably is: ");
                return o_neighbor2;
            }

            ///Check non-ring neighbors of oxygen neighbors (the one without non-ring carbon is anomeric)
            bool neighbor2_is_anomeric = false;
            for(AtomVector::iterator it1 = o_neighbor1_neighbors.begin(); it1 != o_neighbor1_neighbors.end(); it1++)///check if neighbor1 of oxygen has non-ring carbon neighbor
            {
                Atom* neighbor1_neighbor = (*it1);
                if(cycle_atoms_str.find(neighbor1_neighbor->GetId()) == string::npos ///if the neighbor is not one of the cycle atoms
                        && (neighbor1_neighbor->GetName().substr(0,1).compare("C") == 0 )) ///if first element is "C"
                {
                    neighbor2_is_anomeric = true;
                    break;
                }
            }
            bool neighbor1_is_anomeric = false;
            for(AtomVector::iterator it1 = o_neighbor2_neighbors.begin(); it1 != o_neighbor2_neighbors.end(); it1++)///check if neighbor1 of oxygen has non-ring carbon neighbor
            {
                Atom* neighbor2_neighbor = (*it1);
                if(cycle_atoms_str.find(neighbor2_neighbor->GetId()) == string::npos ///if the neighbor is not one of the cycle atoms
                        && (neighbor2_neighbor->GetName().substr(0,1).compare("C") == 0 )) ///if first element is "C"
                {
                    neighbor1_is_anomeric = true;
                }
            }
            if(!neighbor1_is_anomeric)
            {
                anomeric_note->type_ = Glycan::WARNING;
                anomeric_note->category_ = Glycan::ANOMERIC;
                anomeric_note->description_ = "Anomeric oxygen is missing";
                anomeric_carbons_status.push_back("Anomeric carbon probably is: ");
                return o_neighbor2;
            }
            else if(!neighbor2_is_anomeric)
            {
                anomeric_note->type_ = Glycan::WARNING;
                anomeric_note->category_ = Glycan::ANOMERIC;
                anomeric_note->description_ = "Anomeric oxygen is missing";
                anomeric_carbons_status.push_back("Anomeric carbon probably is: ");
                return o_neighbor1;
            }
            anomeric_note->type_ = Glycan::WARNING;
            anomeric_note->category_ = Glycan::ANOMERIC;
            anomeric_note->description_ = "Anomeric oxygen is missing";
            anomeric_carbons_status.push_back("Not enough information to detect the anomeric carbon, it has been chosen randomely: ");
            return o_neighbor1;
        }
    }
    return NULL;
}

Assembly::AtomVector Assembly::SortCycle(AtomVector cycle, Atom *anomeric_atom, stringstream &sorted_cycle_stream)
{
    AtomVector sorted_cycle = AtomVector();
    for(AtomVector::iterator it = cycle.begin(); it != cycle.end(); it++)
    {
        Atom* atom = (*it);
        int index = distance(cycle.begin(), it);
        if(atom->GetId().compare(anomeric_atom->GetId()) == 0)
        {
            if(index == cycle.size() - 1)///anomeric atom is at the end of the cycle
            {
                sorted_cycle.push_back(anomeric_atom);///anomeric atom as the first atom of the cycle
                sorted_cycle_stream << anomeric_atom->GetId() << "-";
                anomeric_atom->SetIsRing(true);
                Atom* a0 = cycle.at(0);
                if(a0->GetName().substr(0,1).compare("O") == 0)///a0 is oxygen so the vector is in reverse order
                {
                    for(AtomVector::iterator it1 = it - 1; it1 != cycle.begin(); it1--)///atoms before the anomeric atom in reverse order
                    {
                        Atom* a = (*it1);
                        sorted_cycle.push_back(a);
                        sorted_cycle_stream << a->GetId() << "-";
                        a->SetIsRing(true);
                    }
                    sorted_cycle.push_back((*cycle.begin()));
                    sorted_cycle_stream << (*cycle.begin())->GetId();
                    (*cycle.begin())->SetIsRing(true);
                }
                else
                {
                    for(AtomVector::iterator it1 = cycle.begin(); it1 != it; it1++)///atoms before the anomeric atom from beginning of vector
                    {
                        Atom* a = (*it1);
                        sorted_cycle.push_back(a);
                        a->SetIsRing(true);
                        if(it1 == it - 1)
                            sorted_cycle_stream << a->GetId();
                        else
                            sorted_cycle_stream << a->GetId() << "-";
                    }
                }
            }
            else///anomeric is not at the end of the cycle
            {
                Atom* next_atom = cycle.at(index + 1);
                if(next_atom->GetName().substr(0,1).compare("O") == 0)///next atom is oxygen so the vector is in reverse order
                {
                    for(AtomVector::iterator it1 = it; it1 != cycle.begin(); it1--) ///atoms befor anomeric atom down to beginning of the vector
                    {
                        Atom* a_before = (*it1);
                        sorted_cycle.push_back(a_before);
                        sorted_cycle_stream << a_before->GetId() << "-";
                        a_before->SetIsRing(true);
                    }
                    sorted_cycle.push_back((*cycle.begin()));
                    sorted_cycle_stream << (*cycle.begin())->GetId() << "-";
                    (*cycle.begin())->SetIsRing(true);
                    for(AtomVector::iterator it2 = cycle.end() - 1; it2 != it; it2--)///atoms from end of the vector down to anomeric atom
                    {
                        Atom* atom_after = (*it2);
                        sorted_cycle.push_back(atom_after);
                        atom_after->SetIsRing(true);
                        if(it2 == it + 1)
                            sorted_cycle_stream << atom_after->GetId();
                        else
                            sorted_cycle_stream << atom_after->GetId() << "-";
                    }
                }
                else///oxygen is before the anomeric atom so the vector is in normal order
                {
                    for(AtomVector::iterator it1 = it; it1 != cycle.end(); it1++) ///atoms after anomeric atom to the end of the vector
                    {
                        Atom* atom_after = (*it1);
                        sorted_cycle.push_back(atom_after);
                        atom_after->SetIsRing(true);
                        //                        if(it1 == cycle.end() - 1)
                        //                            sorted_cycle_stream << atom_after->GetId();
                        //                        else
                        sorted_cycle_stream << atom_after->GetId() << "-";
                    }
                    for(AtomVector::iterator it2 = cycle.begin(); it2 != it; it2++)///atoms befor the anomeric atom from beginning of vector
                    {
                        Atom* atom_before = (*it2);
                        sorted_cycle.push_back(atom_before);
                        atom_before->SetIsRing(true);
                        if(it2 == it - 1)
                            sorted_cycle_stream << atom_before->GetId();
                        else
                            sorted_cycle_stream << atom_before->GetId() << "-";
                    }
                }
            }
        }
    }
    return sorted_cycle;
}

void Assembly::CalculateGlyprobityGeometryOutliers(Monosaccharide* mono)
{
    ///EXTRACTING BOND LENGTHS OF MONOSACCHARIDE ATOM PAIRS AND BOND ANGLES
    vector<string> visited_bonds = vector<string>();
    vector<string> visited_angles = vector<string>();
    vector<double> bond_statistics = vector<double>();
    vector<double> bond_angle_statistics = vector<double>();
    stringstream bond_lengths_stream;
    stringstream bond_angles_stream;
    double bond_angle = 0.0;
    bool is_atom2_ring = false;
    bool is_atom3_ring = false;

    bond_lengths_stream << "GLYPROBITY REPORT" << endl <<
                           "<--Geometry Outliers-->" << endl <<
                           "Bond lengths (current bond length, ontology mean, ontology standard deviation)" << endl;

    bond_angles_stream << "Bond angles (current bond angle, ontology mean, ontology standard deviation)" << endl;

    for(AtomVector::iterator ring_atom_it = mono->cycle_atoms_.begin(); ring_atom_it != mono->cycle_atoms_.end(); ring_atom_it++)
    {
        Atom* atom1 = (*ring_atom_it);
        AtomVector atom1_neighbors = atom1->GetNode()->GetNodeNeighbors();
        for(AtomVector::iterator atom1_neighbors_it = atom1_neighbors.begin(); atom1_neighbors_it != atom1_neighbors.end(); atom1_neighbors_it++)
        {
            Atom* atom2 = (*atom1_neighbors_it);
            stringstream check_bond;
            stringstream check_bond_reverse;

            check_bond << atom1->GetId() << "-" << atom2->GetId();
            check_bond_reverse << atom2->GetId() << "-" << atom1->GetId();
            if(find(visited_bonds.begin(), visited_bonds.end(), check_bond.str()) == visited_bonds.end() &&
                    find(visited_bonds.begin(), visited_bonds.end(), check_bond_reverse.str()) == visited_bonds.end())///if the bond has not been visited before
            {
                visited_bonds.push_back(check_bond.str());
                visited_bonds.push_back(check_bond_reverse.str());
                bond_lengths_stream << atom1->GetName() << "-" << atom2->GetName() << ": " <<
                                       atom1->GetCoordinates().at(model_index_)->Distance(*(atom2->GetCoordinates().at(model_index_)));

                ///if atom 2 is not a ring atom set the flag to true
                if(mono->cycle_atoms_str_.find(atom2->GetId()) == string::npos)
                    is_atom2_ring = false;
                else
                    is_atom2_ring = true;

                ///Find same bonds in the same monosaccharide in the ontology, calculate mean and standard deviation from ontology
                bond_statistics = CalculateBondlengthsStatisticsBasedOnOntologyInfo(atom1->GetName(), atom2->GetName(), is_atom2_ring, mono->sugar_name_.monosaccharide_short_name_);
                bond_lengths_stream << ", " << bond_statistics.at(0) << ", " << bond_statistics.at(1) << endl;
            }

            if(is_atom2_ring)///Calculating bond angles only for ring atoms and combination of ring atoms anf their immediate neighbors
            {
                ///EXTRACTING BOND ANGLES
                AtomVector atom2_neighbors = atom2->GetNode()->GetNodeNeighbors();
                for(AtomVector::iterator atom2_neighbors_it = atom2_neighbors.begin(); atom2_neighbors_it != atom2_neighbors.end(); atom2_neighbors_it++)
                {
                    if((*ring_atom_it) != (*atom2_neighbors_it))///if the neighbor of the second atom is not the first atom we chose
                    {
                        Atom* atom3 = (*atom2_neighbors_it);
                        stringstream check_angle;
                        stringstream check_angle_reverse;
                        check_angle << atom1->GetId() << "-" << atom2->GetId() << "-" << atom3->GetId();
                        check_angle_reverse  << atom3->GetId() << "-" << atom2->GetId() << "-" << atom1->GetId();
                        if(find(visited_angles.begin(), visited_angles.end(), check_angle.str()) == visited_angles.end() &&
                                find(visited_angles.begin(), visited_angles.end(), check_angle_reverse.str()) == visited_angles.end())///if the bond angle has not been visited before
                        {
                            visited_angles.push_back(check_angle.str());
                            visited_angles.push_back(check_angle_reverse.str());

                            bond_angle = CalculateBondAngleByAtoms(atom1, atom2, atom3);
                            bond_angles_stream << atom1->GetName() << "-" << atom2->GetName() << "-" << atom3->GetName() << ": " << bond_angle;

                            ///if atom 3 is not a ring atom set the flag to true
                            if(mono->cycle_atoms_str_.find(atom3->GetId()) == string::npos)
                                is_atom3_ring = false;
                            else
                                is_atom3_ring = true;

                            ///Find same bond angles in the same monosaccharide in the ontology, calculate mean and standard deviation from ontology
                            bond_angle_statistics = CalculateBondAnglesStatisticsBasedOnOntologyInfo(atom1->GetName(), atom2->GetName(), atom3->GetName(), is_atom3_ring,
                                                                                                     mono->sugar_name_.monosaccharide_short_name_);
                            bond_angles_stream << ", " << bond_angle_statistics.at(0) << ", " << bond_angle_statistics.at(1) << endl;
                        }
                    }
                }
            }
        }
    }
    cout << endl << bond_lengths_stream.str() << endl << bond_angles_stream.str() << "<--------------------->" << endl << endl;
}

vector<string> Assembly::GetSideGroupOrientations(Monosaccharide* mono, string cycle_atoms_str)
{
    vector<string> orientations = vector<string>();
    vector<AtomVector> side_atoms = vector<AtomVector>();
    AtomVector default_atom_vector = AtomVector(3, NULL);

    for(AtomVector::iterator it = mono->cycle_atoms_.begin(); it != mono->cycle_atoms_.end() - 1; it++) ///iterate on cycle atoms except the oxygen in the ring
    {
        orientations.push_back("N");
        side_atoms.push_back(default_atom_vector);
        int index = distance(mono->cycle_atoms_.begin(), it);
        Atom* prev_atom = new Atom();
        Atom* current_atom = (*it);
        Atom* next_atom = new Atom();
        if(index == 0)///if the current atom is the anomeric atom
            prev_atom = mono->cycle_atoms_.at(mono->cycle_atoms_.size() - 1); ///previous atom is the oxygen(last atom of the sorted cycle)
        else
            prev_atom = mono->cycle_atoms_.at(index - 1);
        next_atom = mono->cycle_atoms_.at(index + 1);

        ///Calculating the plane based on the two ring neighbors of the current atom
        Coordinate prev_atom_coord = Coordinate(*prev_atom->GetCoordinates().at(model_index_));
        Coordinate current_atom_coord = Coordinate(*current_atom->GetCoordinates().at(model_index_));
        Coordinate next_atom_coord = Coordinate(*next_atom->GetCoordinates().at(model_index_));
        prev_atom_coord.operator -(current_atom_coord) ;
        next_atom_coord.operator -(current_atom_coord) ;
        Plane plane = Plane();
        plane.SetV1(prev_atom_coord);
        plane.SetV2(next_atom_coord);
        Coordinate normal_v = plane.GetUnitNormalVector();

        ///Calculating the orientation of the side atoms
        AtomNode* node = current_atom->GetNode();
        AtomVector neighbors = node->GetNodeNeighbors();
        int not_h_neighbors = 0;
        for(AtomVector::iterator it1 = neighbors.begin(); it1 != neighbors.end(); it1++)
        {
            Atom* neighbor = (*it1);
            string neighbor_id = neighbor->GetId();
            if(cycle_atoms_str.find(neighbor_id) == string::npos) ///if not one of the cycle atoms
            {
                if(neighbor->GetName().at(0) != 'H') ///deoxy check
                    not_h_neighbors++;
            }
        }
        if(not_h_neighbors != 0)
        {
            for(AtomVector::iterator it1 = neighbors.begin(); it1 != neighbors.end(); it1++)
            {
                Atom* neighbor = (*it1);
                string neighbor_id = neighbor->GetId();
                if(cycle_atoms_str.find(neighbor_id) == string::npos) ///if not one of the cycle atoms
                {
                    if(neighbor->GetName().at(0) != 'H') ///deoxy check
                        not_h_neighbors++;
                    Coordinate side_atom_coord = Coordinate(*neighbor->GetCoordinates().at(model_index_));
                    side_atom_coord.operator -(current_atom_coord);
                    side_atom_coord.Normalize();
                    double theta = acos(normal_v.DotProduct(side_atom_coord));

                    if(index == 0 && neighbor_id.at(0) == 'C')///if anomeric atom has a non-ring carbon neighbor
                    {
                        if(orientations.at(index).compare("N") == 0) ///if the position of non-ring oxygen or nitrogen hasn't been set yet
                        {
                            if(theta > (gmml::PI_RADIAN/2))
                            {
                                orientations.at(index) = "-1D";
                                side_atoms.at(index).at(0) = neighbor;
                            }
                            else
                            {
                                orientations.at(index) = "-1U";
                                side_atoms.at(index).at(0) = neighbor;
                            }
                            continue;
                        }
                        else
                        { ///position of non-ring oxygen or nitrogen + the non-ring carbon
                            stringstream ss;
                            if(theta > (gmml::PI_RADIAN/2))
                            {
                                ss << orientations.at(index) << "-1D";
                                orientations.at(index) = ss.str();
                                side_atoms.at(index).at(0) = neighbor;
                            }
                            else
                            {
                                ss << orientations.at(index) << "-1U";
                                orientations.at(index) = ss.str();
                                side_atoms.at(index).at(0) = neighbor;
                            }
                            break;
                        }
                    }
                    else if(neighbor_id.at(0) == 'O' || neighbor_id.at(0) == 'N')///if neighbor is a non-ring oxygen or nitrogen
                    {
                        if(index == 0)///current atom is anomeric
                        {
                            if(orientations.at(index).compare("N") == 0) ///if the position of non-ring carbon neighbor (if exist) hasn't been set yet
                            {
                                if(theta > (gmml::PI_RADIAN/2))
                                {
                                    orientations.at(index) = "D";
                                    side_atoms.at(index).at(1) = neighbor;
                                }
                                else
                                {
                                    orientations.at(index) = "U";
                                    side_atoms.at(index).at(1) = neighbor;
                                }
                                continue;
                            }
                            else
                            { ///position of non-ring oxygen or nitrogen + the non-ring carbon
                                stringstream ss;
                                if(theta > (gmml::PI_RADIAN/2))
                                {
                                    ss << "D" << orientations.at(index);
                                    orientations.at(index) = ss.str();
                                    side_atoms.at(index).at(1) = neighbor;
                                }
                                else
                                {
                                    ss << "U" << orientations.at(index);
                                    orientations.at(index) = ss.str();
                                    side_atoms.at(index).at(1) = neighbor;
                                }
                                break;
                            }
                        }
                        else
                        {
                            if(theta > (gmml::PI_RADIAN/2))
                            {
                                orientations.at(index) = "D";
                                side_atoms.at(index).at(1) = neighbor;
                            }
                            else
                            {
                                orientations.at(index) = "U";
                                side_atoms.at(index).at(1) = neighbor;
                            }
                            break;
                        }
                    }
                    else if(index == mono->cycle_atoms_.size() - 2 && neighbor_id.at(0) == 'C')///if the last ring carbon has a non-ring carbon neighbor
                    {
                        ///Check if neighbor of neighbor is oxygen or nitrogen
                        AtomNode* neighbor_node = neighbor->GetNode();
                        AtomVector neighbors_of_neighbor = neighbor_node->GetNodeNeighbors();
                        int o_neighbors = 0;
                        int not_h_neighbors = 0;
                        for(AtomVector::iterator it2 = neighbors_of_neighbor.begin(); it2 != neighbors_of_neighbor.end(); it2++)
                        {
                            Atom* neighbor_of_neighbor = (*it2);
                            string neighbor_of_neighbor_id = neighbor_of_neighbor->GetId();
                            if(cycle_atoms_str.find(neighbor_of_neighbor_id) == string::npos &&
                                    (neighbor_of_neighbor_id.at(0) == 'O' || neighbor_of_neighbor_id.at(0) == 'N'))///if neighbor of neighbor is a non-ring oxygen or nitrogen
                                o_neighbors++;
                            if(cycle_atoms_str.find(neighbor_of_neighbor_id) == string::npos &&
                                    (neighbor_of_neighbor_id.at(0) != 'H'))///if neighbor of neighbor is any non-ring atom other than hydrogen
                                not_h_neighbors++;
                        }
                        if (o_neighbors >= 1)
                        {
                            if(theta > (gmml::PI_RADIAN/2))
                            {
                                orientations.at(index) = "D";
                                side_atoms.at(index).at(0) = neighbor;
                            }
                            else
                            {
                                orientations.at(index) = "U";
                                side_atoms.at(index).at(0) = neighbor;
                            }
                        }
                        else if(not_h_neighbors == 0)///Type Deoxy
                        {
                            if(theta > (gmml::PI_RADIAN/2))
                            {
                                orientations.at(index) = "Dd";
                                side_atoms.at(index).at(0) = neighbor;
                            }
                            else
                            {
                                orientations.at(index) = "Ud";
                                side_atoms.at(index).at(0) = neighbor;
                            }
                            break;
                        }
                        if(orientations.at(index).compare("N") != 0)
                            break;
                    }
                }
            }
        }
    }
    mono->side_atoms_ = side_atoms;

    return orientations;
}

ChemicalCode* Assembly::BuildChemicalCode(vector<string> orientations)
{
    ChemicalCode* code = new ChemicalCode();
    if(orientations.size() == 5 )
        code->base_ = "P";
    else if(orientations.size() == 4 )
        code->base_ = "F";
    else
        code->base_ = "?";

    ///Side atom(s) of anomeric
    ///Has only non-ring oxygen neighbor
    if(orientations.at(0).compare("U") == 0)
        code->right_up_.push_back("a");
    else if(orientations.at(0).compare("D") == 0)
        code->right_down_.push_back("a");

    ///Has non-ring oxygen and carbon neighbors
    else if(orientations.at(0).compare("U-1U") == 0)
    {
        code->right_up_.push_back("a");
        code->right_up_.push_back("-1");
    }
    else if(orientations.at(0).compare("D-1U") == 0)
    {
        code->right_down_.push_back("a");
        code->right_up_.push_back("-1");
    }
    else if(orientations.at(0).compare("U-1D") == 0)
    {
        code->right_up_.push_back("a");
        code->right_down_.push_back("-1");
    }
    else if(orientations.at(0).compare("D-1D") == 0)
    {
        code->right_down_.push_back("a");
        code->right_down_.push_back("-1");
    }

    ///Has only non-ring carbon neighbor
    else if(orientations.at(0).compare("-1U") == 0)
        code->right_up_.push_back("-1");
    else if(orientations.at(0).compare("-1D") == 0)
        code->right_down_.push_back("-1");

    ///Side atom of other carbons of the ring
    for(vector<string>::iterator it = orientations.begin() + 1; it != orientations.end() - 1; it++)
    {
        string orientation = (*it);
        int index = distance(orientations.begin(), it);
        if(orientation.compare("U") == 0)
            code->left_up_.push_back(gmml::ConvertT(index + 1));
        else if(orientation.compare("D") == 0)
            code->left_down_.push_back(gmml::ConvertT(index + 1));
        else if(orientation.compare("N") == 0)
        {
            stringstream ss;
            ss << gmml::ConvertT(index + 1) << "d";
            code->left_middle_.push_back(ss.str() );
        }
    }

    ///Side atom(s) of last carbon
    if(orientations.at(orientations.size() - 1).compare("U") == 0)
        code->right_up_.push_back("+1");
    else if(orientations.at(orientations.size() - 1).compare("D") == 0)
        code->right_down_.push_back("+1");
    ///Type Deoxy
    else if(orientations.at(orientations.size() - 1).compare("Ud") == 0)
        code->right_up_.push_back("+1d");
    else if(orientations.at(orientations.size() - 1).compare("Dd") == 0)
        code->right_down_.push_back("+1d");

    return code;
}

Assembly::AtomVector Assembly::ExtractAdditionalSideAtoms(Monosaccharide *mono)
{
    AtomVector plus_sides = AtomVector();
    if(mono->side_atoms_.at(mono->side_atoms_.size() - 1).at(0) != NULL)///if there exist a +1 carbon atom. in side_atoms_ structure (vector<AtomVector>) the first index of the last element is dedicated to +1 atom
    {
        plus_sides.push_back(mono->side_atoms_.at(mono->side_atoms_.size() - 1).at(0));
        AtomVector plus_one_atom_neighbors = mono->side_atoms_.at(mono->side_atoms_.size() - 1).at(0)->GetNode()->GetNodeNeighbors();
        for(AtomVector::iterator it1 = plus_one_atom_neighbors.begin(); it1 != plus_one_atom_neighbors.end(); it1++)
        {
            if((*it1)->GetName().at(0) == 'C' && mono->cycle_atoms_str_.find((*it1)->GetId()) == string::npos)///+2 carbon atom found
            {
                Atom* plus_two = (*it1);
                plus_sides.push_back(plus_two);
                mono->side_atoms_.at(mono->side_atoms_.size() - 1).at(1) = plus_two;///in side_atoms_ structure (vector<AtomVector>) the second index of the last element is dedicated to +2 atom

                AtomVector plus_two_atom_neighbors = plus_two->GetNode()->GetNodeNeighbors();
                for(AtomVector::iterator it2 = plus_two_atom_neighbors.begin(); it2 != plus_two_atom_neighbors.end(); it2++)
                {
                    Atom* plus_three = (*it2);
                    if(plus_three->GetName().at(0) == 'C' && plus_three->GetId().compare(mono->side_atoms_.at(mono->side_atoms_.size() - 1).at(0)->GetId()) != 0)///+3 carbon atom found
                    {
                        plus_sides.push_back(plus_three);
                        mono->side_atoms_.at(mono->side_atoms_.size() - 1).at(2) = plus_three;///in side_atoms_ structure (vector<AtomVector>) the third index of the last element is dedicated to +3 atom
                        break;
                    }
                }
                break;
            }
        }
    }
    return plus_sides;
}

void Assembly::ExtractDerivatives(Monosaccharide * mono)
{
    for(AtomVector::iterator it = mono->cycle_atoms_.begin(); it != mono->cycle_atoms_.end() - 1; it++) ///iterate on cycle atoms except the oxygen of the ring
    {
        int index = distance(mono->cycle_atoms_.begin(), it);
        Atom* target = (*it);
        string key = "";
        string value = "";
        AtomVector t_neighbors = target->GetNode()->GetNodeNeighbors();
        for(AtomVector::iterator it1 = t_neighbors.begin(); it1 != t_neighbors.end(); it1++)
        {
            AtomVector pattern_atoms = AtomVector();
            Atom* t_neighbor = (*it1);
            if(t_neighbor->GetName().at(0) == 'N' && mono->cycle_atoms_str_.find(t_neighbor->GetId()) == string::npos)///check formulas with nitrogen
            {
                if((value = CheckxC_N(target, mono->cycle_atoms_str_, pattern_atoms)).compare("") != 0)///xCH-N
                    break;
                if((value = CheckxC_NxO_CO_C(target, mono->cycle_atoms_str_, 'N', pattern_atoms)).compare("") != 0)///xC-N-C=OCH3
                    break;
                if((value = CheckxC_NxO_CO_CO(target, mono->cycle_atoms_str_, 'N', pattern_atoms)).compare("") != 0)///xC-N-C=OCH2OH
                    break;
                if((value = CheckxC_NxO_SO3(target, mono->cycle_atoms_str_, 'N', pattern_atoms)).compare("") != 0)///xC-N-SO3
                    break;
                if((value = CheckxC_NxO_PO3(target, mono->cycle_atoms_str_, 'N', pattern_atoms)).compare("") != 0)///xC-N-PO3
                    break;
                if((value = CheckxC_NxO_C(target, mono->cycle_atoms_str_, 'N', pattern_atoms)).compare("") != 0)///xC-N-CH3
                    break;
            }
            if(t_neighbor->GetName().at(0) == 'O' && mono->cycle_atoms_str_.find(t_neighbor->GetId()) == string::npos)///check formulas with oxygen
            {
                if((value = CheckxC_NxO_CO_C(target, mono->cycle_atoms_str_, 'O', pattern_atoms)).compare("") != 0)///xC-O-C=OCH3
                    break;
                if((value = CheckxC_NxO_CO_CO(target, mono->cycle_atoms_str_, 'O', pattern_atoms)).compare("") != 0)///xC-O-C=OCH2OH
                    break;
                if((value = CheckxC_NxO_SO3(target, mono->cycle_atoms_str_, 'O', pattern_atoms)).compare("") != 0)///xC-O-SO3
                    break;
                if((value = CheckxC_NxO_PO3(target, mono->cycle_atoms_str_, 'O', pattern_atoms)).compare("") != 0)///xC-O-PO3
                    break;
                if((value = CheckxC_NxO_C(target, mono->cycle_atoms_str_, 'O', pattern_atoms)).compare("") != 0)///xC-O-CH3
                    break;
                if((value = CheckxCOO(target, mono->cycle_atoms_str_, pattern_atoms)).compare("") != 0)///xC-(O,O) and xC-(O,OH)
                    break;
            }
        }
        if(value.compare("") != 0)///if any pattern matched add it to the index-derivative map
        {
            if(index == 0)
                key = "a";
            else
                key = gmml::ConvertT(index + 1);
            mono->derivatives_map_[key] = value;
        }
    }
    for(vector<AtomVector>::iterator it = mono->side_atoms_.begin(); it != mono->side_atoms_.end(); it++) ///iterate on side atoms
    {
        int index = distance(mono->side_atoms_.begin(), it);
        int side_branch_last_carbon_index = 0;
        AtomVector sides = (*it);
        Atom* target = NULL;
        string key = "";
        string value = "";
        if(it == mono->side_atoms_.begin())///side atoms of anomeric carbon
        {
            if(sides.at(0) != NULL)
                target = sides.at(0);///first index of each side is for carbon atoms in the vector<AtomVector> structure
        }
        if(it == mono->side_atoms_.end() - 1)//side atoms of last carbon of the ring
        {
            if(sides.at(0) != NULL)
            {
                for(side_branch_last_carbon_index = sides.size() - 1; sides.at(side_branch_last_carbon_index) == NULL; side_branch_last_carbon_index-- ){}
                target = sides.at(side_branch_last_carbon_index);
            }
        }
        if(target != NULL)
        {
            AtomVector t_neighbors = target->GetNode()->GetNodeNeighbors();
            for(AtomVector::iterator it1 = t_neighbors.begin(); it1 != t_neighbors.end(); it1++)
            {
                AtomVector pattern_atoms = AtomVector();
                Atom* t_neighbor = (*it1);
                if(t_neighbor->GetName().at(0) == 'N' && mono->cycle_atoms_str_.find(t_neighbor->GetId()) == string::npos)///check formulas with nitrogen
                {
                    if((value = CheckxC_N(target, mono->cycle_atoms_str_, pattern_atoms)).compare("") != 0)///xCH-N
                        break;
                    if((value = CheckxC_NxO_CO_C(target, mono->cycle_atoms_str_, 'N', pattern_atoms)).compare("") != 0)///xC-N-C=OCH3
                        break;
                    if((value = CheckxC_NxO_CO_CO(target, mono->cycle_atoms_str_, 'N', pattern_atoms)).compare("") != 0)///xC-N-C=OCH2OH
                        break;
                    if((value = CheckxC_NxO_SO3(target, mono->cycle_atoms_str_, 'N', pattern_atoms)).compare("") != 0)///xC-N-SO3
                        break;
                    if((value = CheckxC_NxO_PO3(target, mono->cycle_atoms_str_, 'N', pattern_atoms)).compare("") != 0)///xC-N-PO3
                        break;
                    if((value = CheckxC_NxO_C(target, mono->cycle_atoms_str_, 'N', pattern_atoms)).compare("") != 0)///xC-N-CH3
                        break;
                }
                if(t_neighbor->GetName().at(0) == 'O' && mono->cycle_atoms_str_.find(t_neighbor->GetId()) == string::npos)///check formulas with oxygen
                {
                    if((value = CheckxC_NxO_CO_C(target, mono->cycle_atoms_str_, 'O', pattern_atoms)).compare("") != 0)///xC-O-C=OCH3
                        break;
                    if((value = CheckxC_NxO_CO_CO(target, mono->cycle_atoms_str_, 'O', pattern_atoms)).compare("") != 0)///xC-O-C=OCH2OH
                        break;
                    if((value = CheckxC_NxO_SO3(target, mono->cycle_atoms_str_, 'O', pattern_atoms)).compare("") != 0)///xC-O-SO3
                        break;
                    if((value = CheckxC_NxO_PO3(target, mono->cycle_atoms_str_, 'O', pattern_atoms)).compare("") != 0)///xC-O-PO3
                        break;
                    if((value = CheckxC_NxO_C(target, mono->cycle_atoms_str_, 'O', pattern_atoms)).compare("") != 0)///xC-O-CH3
                        break;
                    if((value = CheckxCOO(target, mono->cycle_atoms_str_, pattern_atoms)).compare("") != 0)///xC-(O,O) and xC-(O,OH)
                        break;
                }
            }
            if(value.compare("") != 0)///if any pattern matched add it to the index-derivative map
            {
                if(index == 0)
                    key = "-1";
                else
                {
                    switch (side_branch_last_carbon_index)
                    {
                        case 0:
                            key = "+1";
                            break;
                        case 1:
                            key = "+2";
                            break;
                        case 2:
                            key = "+3";
                            break;
                    }

                }
                mono->derivatives_map_[key] = value;
            }
        }
    }
}

string Assembly::CheckxC_N(Atom* target, string cycle_atoms_str, AtomVector& pattern_atoms)
{
    stringstream pattern;
    pattern << "xC";
    Atom* N = NULL;
    AtomVector t_neighbors = target->GetNode()->GetNodeNeighbors();
    for(AtomVector::iterator it1 = t_neighbors.begin(); it1 != t_neighbors.end(); it1++)
    {
        Atom* t_neighbor = (*it1);
        if(cycle_atoms_str.find(t_neighbor->GetId()) == string::npos)
        {
            if(t_neighbor->GetName().at(0) != 'N')
                pattern << t_neighbor->GetName().at(0);
            else if(t_neighbor->GetName().at(0) == 'N' && N != NULL)
                pattern << t_neighbor->GetName().at(0);
            else if(t_neighbor->GetName().at(0) == 'N' && N == NULL)
                N = t_neighbor;
        }
    }
    if(N != NULL)
    {
        pattern << "-N";
        AtomVector n_neighbors = N->GetNode()->GetNodeNeighbors();
        for(AtomVector::iterator it = n_neighbors.begin(); it != n_neighbors.end(); it++)
        {
            Atom* n_neighbor = (*it);
            if(n_neighbor->GetId().compare(target->GetId()) != 0)
                pattern << n_neighbor->GetName().at(0);
        }
    }
    //        cout << "CheckxC_N:" << pattern.str() << endl;
    if(pattern.str().compare("xCH-NHH") == 0 || pattern.str().compare("xC-N") == 0 || pattern.str().compare("xCHH-NHH") == 0 || pattern.str().compare("xCH-N") == 0  ||
            pattern.str().compare("xCHH-N") == 0 || pattern.str().compare("xC-NHH") == 0)
        return "xCH-N";
    else
        return "";
}

string Assembly::CheckxC_NxO_CO_C(Atom *target, string cycle_atoms_str, char NxO, AtomVector& pattern_atoms)
{
    stringstream pattern;
    pattern << "xC";
    Atom* N_or_O = NULL;
    AtomVector t_neighbors = target->GetNode()->GetNodeNeighbors();
    for(AtomVector::iterator it = t_neighbors.begin(); it != t_neighbors.end(); it++)
    {
        Atom* t_neighbor = (*it);
        if(cycle_atoms_str.find(t_neighbor->GetId()) == string::npos)
        {
            if(t_neighbor->GetName().at(0) != NxO)
                pattern << t_neighbor->GetName().at(0);
            else if(t_neighbor->GetName().at(0) == NxO && N_or_O != NULL)
                pattern << t_neighbor->GetName().at(0);
            else if(t_neighbor->GetName().at(0) == NxO && N_or_O == NULL)
                N_or_O = t_neighbor;
        }
    }
    if(N_or_O != NULL)
    {
        Atom* C = NULL;
        pattern << "-" << NxO;
        AtomVector n_neighbors = N_or_O->GetNode()->GetNodeNeighbors();
        for(AtomVector::iterator it = n_neighbors.begin(); it != n_neighbors.end(); it++)
        {
            Atom* n_neighbor = (*it);
            if(n_neighbor->GetId().compare(target->GetId()) != 0 && n_neighbor->GetName().at(0) != 'C')
                pattern << n_neighbor->GetName().at(0);
            else if(n_neighbor->GetId().compare(target->GetId()) != 0 && n_neighbor->GetName().at(0) == 'C' && C != NULL)
                pattern << n_neighbor->GetName().at(0);
            else if(n_neighbor->GetId().compare(target->GetId()) != 0 && n_neighbor->GetName().at(0) == 'C' && C == NULL)
                C = n_neighbor;
        }
        if(C != NULL)
        {
            pattern_atoms.push_back(C);
            Atom* CC = NULL;
            Atom* CO = NULL;
            pattern << "-C";
            AtomVector c_neighbors = C->GetNode()->GetNodeNeighbors();
            for(AtomVector::iterator it = c_neighbors.begin(); it != c_neighbors.end(); it++)
            {
                Atom* c_neighbor = (*it);
                if(c_neighbor->GetId().compare(N_or_O->GetId()) != 0 && c_neighbor->GetName().at(0) != 'C' && c_neighbor->GetName().at(0) != 'O')
                    pattern << c_neighbor->GetName().at(0);
                else if(c_neighbor->GetId().compare(N_or_O->GetId()) != 0 &&  c_neighbor->GetName().at(0) == 'C' && CC != NULL)
                    pattern << c_neighbor->GetName().at(0);
                else if(c_neighbor->GetId().compare(N_or_O->GetId()) != 0 &&  c_neighbor->GetName().at(0) == 'O' && CO != NULL)
                    pattern << c_neighbor->GetName().at(0);
                else if(c_neighbor->GetId().compare(N_or_O->GetId()) != 0 &&  c_neighbor->GetName().at(0) == 'O' && CO == NULL)
                    CO = c_neighbor;
                else if(c_neighbor->GetId().compare(N_or_O->GetId()) != 0 &&  c_neighbor->GetName().at(0) == 'C' && CC == NULL)
                    CC = c_neighbor;
            }
            if(CO != NULL)
            {
                pattern_atoms.push_back(CO);
                pattern << "O";
                AtomVector co_neighbors = CO->GetNode()->GetNodeNeighbors();
                for(AtomVector::iterator it = co_neighbors.begin(); it != co_neighbors.end(); it++)
                {
                    Atom* co_neighbor = (*it);
                    if(co_neighbor->GetId().compare(C->GetId()) != 0)
                        pattern << co_neighbor->GetName().at(0);
                }
            }
            if(CC != NULL)
            {
                pattern_atoms.push_back(CC);
                pattern << "-C";
                AtomVector cc_neighbors = CC->GetNode()->GetNodeNeighbors();
                for(AtomVector::iterator it = cc_neighbors.begin(); it != cc_neighbors.end(); it++)
                {
                    Atom* cc_neighbor = (*it);
                    if(cc_neighbor->GetId().compare(C->GetId()) != 0)
                        pattern << cc_neighbor->GetName().at(0);
                }
            }
        }
    }
    //            cout << "CheckxC_NxO_CO_C:" << pattern.str() << endl;
    if(NxO == 'N')
    {
        if(pattern.str().compare("xCH-NH-CO-CHHH") == 0 || pattern.str().compare("xC-N-CO-C") == 0 || pattern.str().compare("xCHH-NH-CO-CHHH") == 0 || pattern.str().compare("xC-NH-CO-CHHH") == 0 ||
                pattern.str().compare("xC-N-CO-CHHH") == 0 || pattern.str().compare("xC-NH-CO-C") == 0 || pattern.str().compare("xCHH-N-CO-CHHH") == 0 || pattern.str().compare("xCHH-NH-CO-C") == 0 ||
                pattern.str().compare("xCHH-N-CO-C") == 0 || pattern.str().compare("xCH-N-CO-CHHH") == 0 || pattern.str().compare("xCH-NH-CO-C") == 0 || pattern.str().compare("xCH-N-CO-C") == 0 )
            return "xC-N-C=OCH3";
        else
            return "";
    }
    else if(NxO == 'O')
    {
        if(pattern.str().compare("xCH-OH-CO-CHHH") == 0 || pattern.str().compare("xC-O-CO-C") == 0 || pattern.str().compare("xCHH-OH-CO-CHHH") == 0 || pattern.str().compare("xC-OH-CO-CHHH") == 0 ||
                pattern.str().compare("xC-O-CO-CHHH") == 0 || pattern.str().compare("xC-OH-CO-C") == 0 || pattern.str().compare("xCHH-O-CO-CHHH") == 0 || pattern.str().compare("xCHH-OH-CO-C") == 0 ||
                pattern.str().compare("xCHH-O-CO-C") == 0 || pattern.str().compare("xCH-O-CO-CHHH") == 0 || pattern.str().compare("xCH-OH-CO-C") == 0 || pattern.str().compare("xCH-O-CO-C") == 0)
            return "xC-O-C=OCH3";
        else
            return "";
    }
    else
        return "";
}

string Assembly::CheckxC_NxO_CO_CO(Atom *target, string cycle_atoms_str, char NxO, AtomVector& pattern_atoms)
{
    stringstream pattern;
    pattern << "xC";
    Atom* N_or_O = NULL;
    AtomVector t_neighbors = target->GetNode()->GetNodeNeighbors();
    for(AtomVector::iterator it = t_neighbors.begin(); it != t_neighbors.end(); it++)
    {
        Atom* t_neighbor = (*it);
        if(cycle_atoms_str.find(t_neighbor->GetId()) == string::npos)
        {
            if(t_neighbor->GetName().at(0) != NxO)
                pattern << t_neighbor->GetName().at(0);
            else if(t_neighbor->GetName().at(0) == NxO && N_or_O != NULL)
                pattern << t_neighbor->GetName().at(0);
            else if(t_neighbor->GetName().at(0) == NxO && N_or_O == NULL)
                N_or_O = t_neighbor;
        }
    }
    if(N_or_O != NULL)
    {
        Atom* C = NULL;
        pattern << "-" << NxO;
        AtomVector n_neighbors = N_or_O->GetNode()->GetNodeNeighbors();
        for(AtomVector::iterator it = n_neighbors.begin(); it != n_neighbors.end(); it++)
        {
            Atom* n_neighbor = (*it);
            if(n_neighbor->GetId().compare(target->GetId()) != 0 && n_neighbor->GetName().at(0) != 'C')
                pattern << n_neighbor->GetName().at(0);
            else if(n_neighbor->GetId().compare(target->GetId()) != 0 && n_neighbor->GetName().at(0) == 'C' && C != NULL)
                pattern << n_neighbor->GetName().at(0);
            else if(n_neighbor->GetId().compare(target->GetId()) != 0 && n_neighbor->GetName().at(0) == 'C' && C == NULL)
                C = n_neighbor;
        }
        if(C != NULL)
        {
            pattern_atoms.push_back(C);
            Atom* CC = NULL;
            Atom* CO = NULL;
            pattern << "-C";
            AtomVector c_neighbors = C->GetNode()->GetNodeNeighbors();
            for(AtomVector::iterator it = c_neighbors.begin(); it != c_neighbors.end(); it++)
            {
                Atom* c_neighbor = (*it);
                if(c_neighbor->GetId().compare(N_or_O->GetId()) != 0 && c_neighbor->GetName().at(0) != 'C' && c_neighbor->GetName().at(0) != 'O')
                    pattern << c_neighbor->GetName().at(0);
                else if(c_neighbor->GetId().compare(N_or_O->GetId()) != 0 &&  c_neighbor->GetName().at(0) == 'C' && CC != NULL)
                    pattern << c_neighbor->GetName().at(0);
                else if(c_neighbor->GetId().compare(N_or_O->GetId()) != 0 &&  c_neighbor->GetName().at(0) == 'O' && CO != NULL)
                    pattern << c_neighbor->GetName().at(0);
                else if(c_neighbor->GetId().compare(N_or_O->GetId()) != 0 &&  c_neighbor->GetName().at(0) == 'O' && CO == NULL)
                    CO = c_neighbor;
                else if(c_neighbor->GetId().compare(N_or_O->GetId()) != 0 &&  c_neighbor->GetName().at(0) == 'C' && CC == NULL)
                    CC = c_neighbor;
            }
            if(CO != NULL)
            {
                pattern_atoms.push_back(CO);
                pattern << "O";
                AtomVector co_neighbors = CO->GetNode()->GetNodeNeighbors();
                for(AtomVector::iterator it = co_neighbors.begin(); it != co_neighbors.end(); it++)
                {
                    Atom* co_neighbor = (*it);
                    if(co_neighbor->GetId().compare(C->GetId()) != 0)
                        pattern << co_neighbor->GetName().at(0);
                }
            }
            if(CC != NULL)
            {
                pattern_atoms.push_back(CC);
                Atom* O = NULL;
                pattern << "-C";
                AtomVector cc_neighbors = CC->GetNode()->GetNodeNeighbors();
                for(AtomVector::iterator it = cc_neighbors.begin(); it != cc_neighbors.end(); it++)
                {
                    Atom* cc_neighbor = (*it);
                    if(cc_neighbor->GetId().compare(C->GetId()) != 0 && cc_neighbor->GetName().at(0) != 'O')
                        pattern << cc_neighbor->GetName().at(0);
                    else if(cc_neighbor->GetId().compare(C->GetId()) != 0 && cc_neighbor->GetName().at(0) == 'O' && O != NULL)
                        pattern << cc_neighbor->GetName().at(0);
                    else if(cc_neighbor->GetId().compare(C->GetId()) != 0 && cc_neighbor->GetName().at(0) == 'O' && O == NULL)
                        O = cc_neighbor;
                }
                if(O != NULL)
                {
                    pattern_atoms.push_back(O);
                    pattern << "-O";
                    AtomVector o_neighbors = O->GetNode()->GetNodeNeighbors();
                    for(AtomVector::iterator it = o_neighbors.begin(); it != o_neighbors.end(); it++)
                    {
                        Atom* o_neighbor = (*it);
                        if(o_neighbor->GetId().compare(CC->GetId()) != 0)
                            pattern << o_neighbor->GetName().at(0);
                    }
                }
            }
        }
    }
    //        cout << "CheckxC_NxO_CO_CO: " << pattern.str() << endl;
    if(NxO == 'N')
    {
        if(pattern.str().compare("xCH-NH-CO-CHH-OH") == 0 || pattern.str().compare("xCH-N-CO-CHH-OH") == 0 || pattern.str().compare("xCH-NH-CO-C-OH") == 0 ||
                pattern.str().compare("xCH-NH-CO-CHH-O") == 0 || pattern.str().compare("xCH-N-CO-C-OH") == 0 || pattern.str().compare("xCH-N-CO-CHH-O") == 0 ||
                pattern.str().compare("xCH-NH-CO-C-OH") == 0 || pattern.str().compare("xCH-N-CO-C-O") == 0 || pattern.str().compare("xC-N-CO-C-O") == 0 ||
                pattern.str().compare("xC-NH-CO-C-O") == 0 || pattern.str().compare("xC-NH-CO-C-OH") == 0 || pattern.str().compare("xC-NH-CO-CHH-O") == 0 ||
                pattern.str().compare("xC-N-CO-CHH-O") == 0 || pattern.str().compare("xC-N-CO-C-OH") == 0 || pattern.str().compare("xC-NH-CO-CHH-OH") == 0 ||
                pattern.str().compare("xC-N-CO-CHH-OH") == 0 || pattern.str().compare("xCHH-NH-CO-CHH-OH") == 0 || pattern.str().compare("xCHH-N-CO-C-O") == 0 ||
                pattern.str().compare("xCHH-N-CO-CHH-OH") == 0 || pattern.str().compare("xCHH-N-CO-CHH-O") == 0 || pattern.str().compare("xCHH-N-CO-C-OH") == 0 ||
                pattern.str().compare("xCHH-NH-CO-C-OH") == 0 || pattern.str().compare("xCHH-NH-CO-CHH-O") == 0 || pattern.str().compare("xCHH-NH-CO-C-O") == 0 )
            return "xC-N-C=OCH2OH";
        else
            return "";
    }
    else if(NxO == 'O')
    {
        if(pattern.str().compare("xCH-OH-CO-CHH-OH") == 0 || pattern.str().compare("xCH-O-CO-CHH-OH") == 0 || pattern.str().compare("xCH-OH-CO-C-OH") == 0 ||
                pattern.str().compare("xCH-OH-CO-CHH-O") == 0 || pattern.str().compare("xCH-O-CO-C-OH") == 0 || pattern.str().compare("xCH-O-CO-CHH-O") == 0 ||
                pattern.str().compare("xCH-OH-CO-C-OH") == 0 || pattern.str().compare("xCH-O-CO-C-O") == 0 || pattern.str().compare("xC-O-CO-C-O") == 0 ||
                pattern.str().compare("xC-OH-CO-C-O") == 0 || pattern.str().compare("xC-OH-CO-C-OH") == 0 || pattern.str().compare("xC-OH-CO-CHH-O") == 0 ||
                pattern.str().compare("xC-O-CO-CHH-O") == 0 || pattern.str().compare("xC-O-CO-C-OH") == 0 || pattern.str().compare("xC-OH-CO-CHH-OH") == 0 ||
                pattern.str().compare("xC-O-CO-CHH-OH") == 0 || pattern.str().compare("xCHH-OH-CO-CHH-OH") == 0 || pattern.str().compare("xCHH-O-CO-C-O") == 0 ||
                pattern.str().compare("xCHH-O-CO-CHH-OH") == 0 || pattern.str().compare("xCHH-O-CO-CHH-O") == 0 || pattern.str().compare("xCHH-O-CO-C-OH") == 0 ||
                pattern.str().compare("xCHH-OH-CO-C-OH") == 0 || pattern.str().compare("xCHH-OH-CO-CHH-O") == 0 || pattern.str().compare("xCHH-OH-CO-C-O") == 0 )

            return "xC-O-C=OCH2OH";
        else
            return "";
    }
    else
        return "";
}

string Assembly::CheckxC_NxO_SO3(Atom *target, string cycle_atoms_str, char NxO, AtomVector& pattern_atoms)
{
    stringstream pattern;
    pattern << "xC";
    Atom* N_or_O = NULL;
    AtomVector t_neighbors = target->GetNode()->GetNodeNeighbors();
    for(AtomVector::iterator it = t_neighbors.begin(); it != t_neighbors.end(); it++)
    {
        Atom* t_neighbor = (*it);
        if(cycle_atoms_str.find(t_neighbor->GetId()) == string::npos)
        {
            if(t_neighbor->GetName().at(0) != NxO)
                pattern << t_neighbor->GetName().at(0);
            else if(t_neighbor->GetName().at(0) == NxO && N_or_O != NULL)
                pattern << t_neighbor->GetName().at(0);
            else if(t_neighbor->GetName().at(0) == NxO && N_or_O == NULL)
                N_or_O = t_neighbor;
        }
    }
    if(N_or_O != NULL)
    {
        Atom* S = NULL;
        pattern << "-" << NxO;
        AtomVector n_neighbors = N_or_O->GetNode()->GetNodeNeighbors();
        for(AtomVector::iterator it = n_neighbors.begin(); it != n_neighbors.end(); it++)
        {
            Atom* n_neighbor = (*it);
            if(n_neighbor->GetId().compare(target->GetId()) != 0 && n_neighbor->GetName().at(0) != 'S')
                pattern << n_neighbor->GetName().at(0);
            else if(n_neighbor->GetId().compare(target->GetId()) != 0 && n_neighbor->GetName().at(0) == 'S' && S != NULL)
                pattern << n_neighbor->GetName().at(0);
            else if(n_neighbor->GetId().compare(target->GetId()) != 0 && n_neighbor->GetName().at(0) == 'S' && S == NULL)
                S = n_neighbor;
        }
        if(S != NULL)
        {
            pattern_atoms.push_back(S);
            Atom* O1 = NULL;
            Atom* O2 = NULL;
            Atom* O3 = NULL;
            pattern << "-S";
            AtomVector s_neighbors = S->GetNode()->GetNodeNeighbors();
            for(AtomVector::iterator it = s_neighbors.begin(); it != s_neighbors.end(); it++)
            {
                Atom* s_neighbor = (*it);
                if(s_neighbor->GetId().compare(N_or_O->GetId()) != 0 && s_neighbor->GetName().at(0) != 'O')
                    pattern << s_neighbor->GetName().at(0);
                else if(s_neighbor->GetId().compare(N_or_O->GetId()) != 0 &&  s_neighbor->GetName().at(0) == 'O' && O1 == NULL)
                    O1 = s_neighbor;
                else if(s_neighbor->GetId().compare(N_or_O->GetId()) != 0 &&  s_neighbor->GetName().at(0) == 'O' && O2 == NULL)
                    O2 = s_neighbor;
                else if(s_neighbor->GetId().compare(N_or_O->GetId()) != 0 &&  s_neighbor->GetName().at(0) == 'O' && O3 == NULL)
                    O3 = s_neighbor;
            }
            if(O1 != NULL && O2 != NULL && O3 != NULL)
            {
                pattern_atoms.push_back(O1);
                pattern_atoms.push_back(O2);
                pattern_atoms.push_back(O3);
                pattern << "OOO";
                AtomVector o1_neighbors = O1->GetNode()->GetNodeNeighbors();
                for(AtomVector::iterator it = o1_neighbors.begin(); it != o1_neighbors.end(); it++)
                {
                    Atom* o1_neighbor = (*it);
                    if(o1_neighbor->GetId().compare(S->GetId()) != 0)
                        pattern << o1_neighbor->GetName().at(0);
                }
                AtomVector o2_neighbors = O2->GetNode()->GetNodeNeighbors();
                for(AtomVector::iterator it = o2_neighbors.begin(); it != o2_neighbors.end(); it++)
                {
                    Atom* o2_neighbor = (*it);
                    if(o2_neighbor->GetId().compare(S->GetId()) != 0)
                        pattern << o2_neighbor->GetName().at(0);
                }
                AtomVector o3_neighbors = O3->GetNode()->GetNodeNeighbors();
                for(AtomVector::iterator it = o3_neighbors.begin(); it != o3_neighbors.end(); it++)
                {
                    Atom* o3_neighbor = (*it);
                    if(o3_neighbor->GetId().compare(S->GetId()) != 0)
                        pattern << o3_neighbor->GetName().at(0);
                }
            }
        }
    }
    //        cout << "CheckxC_NxO_SO3: " << pattern.str() << endl;
    if(NxO == 'N')
    {
        if(pattern.str().compare("xCH-NH-SOOOH") == 0 || pattern.str().compare("xCH-N-SOOOH") == 0 || pattern.str().compare("xCH-NH-SOOO") == 0 || pattern.str().compare("xCH-N-SOOO") == 0 ||
                pattern.str().compare("xCHH-NH-SOOOH") == 0 || pattern.str().compare("xCHH-NH-SOOO") == 0 || pattern.str().compare("xCHH-N-SOOOH") == 0 || pattern.str().compare("xCHH-N-SOOO") == 0 ||
                pattern.str().compare("xC-N-SOOO") == 0 || pattern.str().compare("xC-NH-SOOOH") == 0 || pattern.str().compare("xC-N-SOOOH") == 0 || pattern.str().compare("xC-NH-SOOO") == 0 )
            return "xC-N-SO3";
        else
            return "";
    }
    else if(NxO == 'O')
    {
        if(pattern.str().compare("xCH-OH-SOOOH") == 0 || pattern.str().compare("xCH-O-SOOOH") == 0 || pattern.str().compare("xCH-OH-SOOO") == 0 || pattern.str().compare("xCH-O-SOOO") == 0 ||
                pattern.str().compare("xCHH-OH-SOOOH") == 0 || pattern.str().compare("xCHH-OH-SOOO") == 0 || pattern.str().compare("xCHH-O-SOOOH") == 0 || pattern.str().compare("xCHH-O-SOOO") == 0 ||
                pattern.str().compare("xC-O-SOOO") == 0 || pattern.str().compare("xC-OH-SOOOH") == 0 || pattern.str().compare("xC-O-SOOOH") == 0 || pattern.str().compare("xC-OH-SOOO") == 0 )
            return "xC-O-SO3";
        else
            return "";
    }
    else
        return "";
}

string Assembly::CheckxC_NxO_PO3(Atom *target, string cycle_atoms_str, char NxO, AtomVector& pattern_atoms)
{
    stringstream pattern;
    pattern << "xC";
    Atom* N_or_O = NULL;
    AtomVector t_neighbors = target->GetNode()->GetNodeNeighbors();
    for(AtomVector::iterator it = t_neighbors.begin(); it != t_neighbors.end(); it++)
    {
        Atom* t_neighbor = (*it);
        if(cycle_atoms_str.find(t_neighbor->GetId()) == string::npos)
        {
            if(t_neighbor->GetName().at(0) != NxO)
                pattern << t_neighbor->GetName().at(0);
            else if(t_neighbor->GetName().at(0) == NxO && N_or_O != NULL)
                pattern << t_neighbor->GetName().at(0);
            else if(t_neighbor->GetName().at(0) == NxO && N_or_O == NULL)
                N_or_O = t_neighbor;
        }
    }
    if(N_or_O != NULL)
    {
        Atom* P = NULL;
        pattern << "-" << NxO;
        AtomVector n_neighbors = N_or_O->GetNode()->GetNodeNeighbors();
        for(AtomVector::iterator it = n_neighbors.begin(); it != n_neighbors.end(); it++)
        {
            Atom* n_neighbor = (*it);
            if(n_neighbor->GetId().compare(target->GetId()) != 0 && n_neighbor->GetName().at(0) != 'P')
                pattern << n_neighbor->GetName().at(0);
            else if(n_neighbor->GetId().compare(target->GetId()) != 0 && n_neighbor->GetName().at(0) == 'P' && P != NULL)
                pattern << n_neighbor->GetName().at(0);
            else if(n_neighbor->GetId().compare(target->GetId()) != 0 && n_neighbor->GetName().at(0) == 'P' && P == NULL)
                P = n_neighbor;
        }
        if(P != NULL)
        {
            pattern_atoms.push_back(P);
            Atom* O1 = NULL;
            Atom* O2 = NULL;
            Atom* O3 = NULL;
            pattern << "-P";
            AtomVector p_neighbors = P->GetNode()->GetNodeNeighbors();
            for(AtomVector::iterator it = p_neighbors.begin(); it != p_neighbors.end(); it++)
            {
                Atom* p_neighbor = (*it);
                if(p_neighbor->GetId().compare(N_or_O->GetId()) != 0 && p_neighbor->GetName().at(0) != 'O')
                    pattern << p_neighbor->GetName().at(0);
                else if(p_neighbor->GetId().compare(N_or_O->GetId()) != 0 &&  p_neighbor->GetName().at(0) == 'O' && O1 == NULL)
                    O1 = p_neighbor;
                else if(p_neighbor->GetId().compare(N_or_O->GetId()) != 0 &&  p_neighbor->GetName().at(0) == 'O' && O2 == NULL)
                    O2 = p_neighbor;
                else if(p_neighbor->GetId().compare(N_or_O->GetId()) != 0 &&  p_neighbor->GetName().at(0) == 'O' && O3 == NULL)
                    O3 = p_neighbor;
            }
            if(O1 != NULL && O2 != NULL && O3 != NULL)
            {
                pattern_atoms.push_back(O1);
                pattern_atoms.push_back(O2);
                pattern_atoms.push_back(O3);
                pattern << "OOO";
                AtomVector o1_neighbors = O1->GetNode()->GetNodeNeighbors();
                for(AtomVector::iterator it = o1_neighbors.begin(); it != o1_neighbors.end(); it++)
                {
                    Atom* o1_neighbor = (*it);
                    if(o1_neighbor->GetId().compare(P->GetId()) != 0)
                        pattern << o1_neighbor->GetName().at(0);
                }
                AtomVector o2_neighbors = O2->GetNode()->GetNodeNeighbors();
                for(AtomVector::iterator it = o2_neighbors.begin(); it != o2_neighbors.end(); it++)
                {
                    Atom* o2_neighbor = (*it);
                    if(o2_neighbor->GetId().compare(P->GetId()) != 0)
                        pattern << o2_neighbor->GetName().at(0);
                }
                AtomVector o3_neighbors = O3->GetNode()->GetNodeNeighbors();
                for(AtomVector::iterator it = o3_neighbors.begin(); it != o3_neighbors.end(); it++)
                {
                    Atom* o3_neighbor = (*it);
                    if(o3_neighbor->GetId().compare(P->GetId()) != 0)
                        pattern << o3_neighbor->GetName().at(0);
                }
            }
        }
    }
    //        cout << "CheckxC_NxO_PO3: " << pattern.str() << endl;
    if(NxO == 'N')
    {
        if(pattern.str().compare("xCH-NH-POOOH") == 0 || pattern.str().compare("xCH-N-POOOH") == 0 || pattern.str().compare("xCH-NH-POOO") == 0 || pattern.str().compare("xCH-N-POOO") == 0 ||
                pattern.str().compare("xCHH-NH-POOOH") == 0 || pattern.str().compare("xCHH-NH-POOO") == 0 || pattern.str().compare("xCHH-N-POOOH") == 0 || pattern.str().compare("xCHH-N-POOO") == 0 ||
                pattern.str().compare("xC-N-POOO") == 0 || pattern.str().compare("xC-NH-POOOH") == 0 || pattern.str().compare("xC-N-POOOH") == 0 || pattern.str().compare("xC-NH-POOO") == 0 )
            return "xC-N-PO3";
        else
            return "";
    }
    else if(NxO = 'O')
    {
        if(pattern.str().compare("xCH-OH-POOOH") == 0 || pattern.str().compare("xCH-O-POOOH") == 0 || pattern.str().compare("xCH-OH-POOO") == 0 || pattern.str().compare("xCH-O-POOO") == 0 ||
                pattern.str().compare("xCHH-OH-POOOH") == 0 || pattern.str().compare("xCHH-OH-POOO") == 0 || pattern.str().compare("xCHH-O-POOOH") == 0 || pattern.str().compare("xCHH-O-POOO") == 0 ||
                pattern.str().compare("xC-O-POOO") == 0 || pattern.str().compare("xC-OH-POOOH") == 0 || pattern.str().compare("xC-O-POOOH") == 0 || pattern.str().compare("xC-OH-POOO") == 0 )            return "xC-O-PO3";
        else
            return "";
    }
    else
        return "";
}

string Assembly::CheckxC_NxO_C(Atom *target, string cycle_atoms_str, char NxO, AtomVector& pattern_atoms)
{
    stringstream pattern;
    pattern << "xC";
    Atom* N_or_O = NULL;
    AtomVector t_neighbors = target->GetNode()->GetNodeNeighbors();
    for(AtomVector::iterator it = t_neighbors.begin(); it != t_neighbors.end(); it++)
    {
        Atom* t_neighbor = (*it);
        if(cycle_atoms_str.find(t_neighbor->GetId()) == string::npos)
        {
            if(t_neighbor->GetName().at(0) != NxO)
                pattern << t_neighbor->GetName().at(0);
            else if(t_neighbor->GetName().at(0) == NxO && N_or_O != NULL)
                pattern << t_neighbor->GetName().at(0);
            else if(t_neighbor->GetName().at(0) == NxO && N_or_O == NULL)
                N_or_O = t_neighbor;
        }
    }
    if(N_or_O != NULL)
    {
        Atom* C = NULL;
        pattern << "-" << NxO;
        AtomVector n_neighbors = N_or_O->GetNode()->GetNodeNeighbors();
        for(AtomVector::iterator it = n_neighbors.begin(); it != n_neighbors.end(); it++)
        {
            Atom* n_neighbor = (*it);
            if(n_neighbor->GetId().compare(target->GetId()) != 0 && n_neighbor->GetName().at(0) != 'C')
                pattern << n_neighbor->GetName().at(0);
            else if(n_neighbor->GetId().compare(target->GetId()) != 0 && n_neighbor->GetName().at(0) == 'C' && C != NULL)
                pattern << n_neighbor->GetName().at(0);
            else if(n_neighbor->GetId().compare(target->GetId()) != 0 && n_neighbor->GetName().at(0) == 'C' && C == NULL)
                C = n_neighbor;
        }
        if(C != NULL)
        {
            pattern_atoms.push_back(C);
            pattern << "-C";
            AtomVector c_neighbors = C->GetNode()->GetNodeNeighbors();
            for(AtomVector::iterator it = c_neighbors.begin(); it != c_neighbors.end(); it++)
            {
                Atom* c_neighbor = (*it);
                if(c_neighbor->GetId().compare(N_or_O->GetId()) != 0)
                    pattern << c_neighbor->GetName().at(0);
            }
        }
    }
    //            cout << "CheckxC_NxO_C: " << pattern.str() << endl;
    if(NxO == 'N')
    {
        if(pattern.str().compare("xCH-N-CHHH") == 0 ||  pattern.str().compare("xCH-NH-CHHH") == 0 || pattern.str().compare("xCH-NH-C") == 0 || pattern.str().compare("xCH-N-C") == 0 ||
                pattern.str().compare("xCHH-N-CHHH") == 0 || pattern.str().compare("xCHH-NH-CHHH") == 0 || pattern.str().compare("xCHH-NH-C") == 0 || pattern.str().compare("xCHH-N-C") == 0 ||
                pattern.str().compare("xC-NH-CHHH") == 0 || pattern.str().compare("xC-NH-C") == 0 || pattern.str().compare("xC-N-CHHH") == 0 || pattern.str().compare("xC-N-C") == 0)
            return "xC-N-CH3";
        else
            return "";
    }
    else if(NxO == 'O')
    {
        if(pattern.str().compare("xCH-O-CHHH") == 0 ||  pattern.str().compare("xCH-OH-CHHH") == 0 || pattern.str().compare("xCH-OH-C") == 0 || pattern.str().compare("xCH-O-C") == 0 ||
                pattern.str().compare("xCHH-O-CHHH") == 0 || pattern.str().compare("xCHH-OH-CHHH") == 0 || pattern.str().compare("xCHH-OH-C") == 0 || pattern.str().compare("xCHH-O-C") == 0 ||
                pattern.str().compare("xC-OH-CHHH") == 0 || pattern.str().compare("xC-OH-C") == 0 || pattern.str().compare("xC-O-CHHH") == 0 || pattern.str().compare("xC-O-C") == 0)
            return "xC-O-CH3";
        else
            return "";
    }
    else
        return "";
}

string Assembly::CheckxCOO(Atom *target, string cycle_atoms_str, AtomVector& pattern_atoms)
{
    stringstream pattern;
    pattern << "xC";
    Atom* O1 = NULL;
    Atom* O2 = NULL;
    AtomVector t_neighbors = target->GetNode()->GetNodeNeighbors();
    for(AtomVector::iterator it = t_neighbors.begin(); it != t_neighbors.end(); it++)
    {
        Atom* t_neighbor = (*it);
        if(cycle_atoms_str.find(t_neighbor->GetId()) == string::npos)
        {
            if(t_neighbor->GetName().at(0) != 'O')
                pattern << t_neighbor->GetName().at(0);
            else if(t_neighbor->GetName().at(0) == 'O' && O1 == NULL)
                O1 = t_neighbor;
            else if(t_neighbor->GetName().at(0) == 'O' && O2 == NULL)
                O2 = t_neighbor;
        }
    }
    if(O1 != NULL && O2 != NULL)
    {
        pattern << "OO";
        AtomVector o1_neighbors = O1->GetNode()->GetNodeNeighbors();
        for(AtomVector::iterator it = o1_neighbors.begin(); it != o1_neighbors.end(); it++)
        {
            Atom* o1_neighbor = (*it);
            if(o1_neighbor->GetId().compare(target->GetId()) != 0)
                pattern << o1_neighbor->GetName().at(0);
        }
        AtomVector o2_neighbors = O2->GetNode()->GetNodeNeighbors();
        for(AtomVector::iterator it = o2_neighbors.begin(); it != o2_neighbors.end(); it++)
        {
            Atom* o2_neighbor = (*it);
            if(o2_neighbor->GetId().compare(target->GetId()) != 0)
                pattern << o2_neighbor->GetName().at(0);
        }
    }
    //        cout << "CheckxCOO: " << pattern.str() << endl;
    if(pattern.str().compare("xCOO") == 0 || pattern.str().compare("xCHOO") == 0)
        return "xC-(O,O)";
    else if(pattern.str().compare("xCOOH") == 0 || pattern.str().compare("xCHOOH") == 0)
        return "xC-(O,OH)";
    else
        return "";
}

void Assembly::Ionizing(string ion_name, string lib_file, string parameter_file, int ion_count)
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

void Assembly::Solvation(double extension, double closeness, string lib_file)
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

double Assembly::GetTotalCharge()
{
    double charge = 0;
    AtomVector all_atoms_of_assembly = this->GetAllAtomsOfAssembly();
    for(AtomVector::iterator it = all_atoms_of_assembly.begin(); it != all_atoms_of_assembly.end(); it++)
    {
        Atom* atom = *it;
        if(atom->MolecularDynamicAtom::GetCharge() != dNotSet)
            charge += atom->MolecularDynamicAtom::GetCharge();

    }
    return charge;
}

double Assembly::GetRadius()
{
    double radius = -INFINITY;
    Coordinate* geometric_center = new Coordinate();
    this->GetCenterOfGeometry(geometric_center);
    AtomVector all_atoms_of_assembly = this->GetAllAtomsOfAssembly();
    for(AtomVector::iterator it = all_atoms_of_assembly.begin(); it != all_atoms_of_assembly.end(); it++)
    {
        Atom* atom = *it;
        double dist = sqrt((geometric_center->GetX() - atom->GetCoordinates().at(model_index_)->GetX()) *
                           (geometric_center->GetX() - atom->GetCoordinates().at(model_index_)->GetX()) +
                           (geometric_center->GetY() - atom->GetCoordinates().at(model_index_)->GetY()) *
                           (geometric_center->GetY() - atom->GetCoordinates().at(model_index_)->GetY()) +
                           (geometric_center->GetZ() - atom->GetCoordinates().at(model_index_)->GetZ()) *
                           (geometric_center->GetZ() - atom->GetCoordinates().at(model_index_)->GetZ()));
        double atom_radius = atom->MolecularDynamicAtom::GetRadius();
        if(atom_radius == dNotSet)
            atom_radius = MINIMUM_RADIUS;
        double dist_to_edge = 0;
        if(atom_radius != dNotSet)
            dist_to_edge = dist + atom_radius;
        else
            dist_to_edge = dist;
        if(dist_to_edge > radius)
            radius = dist_to_edge;
    }
    return radius;
}

double Assembly::GetTotalMass()
{
    double mass = 0.0;
    AtomVector all_atoms_of_assembly = this->GetAllAtomsOfAssembly();
    for(AtomVector::iterator it = all_atoms_of_assembly.begin(); it != all_atoms_of_assembly.end(); it++)
    {
        Atom* atom = *it;
        mass += atom->MolecularDynamicAtom::GetMass();
    }
    return mass;
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

void Assembly::GenerateCompleteSugarName(Monosaccharide *mono)
{
    stringstream in_bracket;
    stringstream head;
    stringstream tail;
    bool minus_one = false;
    if(mono->derivatives_map_.find("-1") != mono->derivatives_map_.end())
        minus_one = true;
    for(map<string, string>::iterator it1 = mono->derivatives_map_.begin(); it1 != mono->derivatives_map_.end(); it1++)
    {
        string key = (*it1).first;
        string value = (*it1).second;
        string long_name_pattern = "";
        string cond_name_pattern = "";
        string long_name_pattern_at_minus_one = "";
        string long_name_pattern_at_plus_one = "";
        string pattern = "";
        if(value.compare("xCH-N") == 0)
        {
            long_name_pattern = "-osamine";
            cond_name_pattern = "N";
            pattern = "CH-N";
            AddModificationRuleOneInfo(key, pattern, mono, long_name_pattern, cond_name_pattern, head, tail, minus_one, in_bracket);
        }
        if(value.compare("xC-N-C=OCH3") == 0)
        {
            long_name_pattern = "N-acetyl-";
            cond_name_pattern = "NAc";
            pattern = "C-N-C=OCH3";
            AddModificationRuleOneInfo(key, pattern, mono, long_name_pattern, cond_name_pattern, head, tail, minus_one, in_bracket);
        }
        if(value.compare("xC-N-C=OCH2OH") == 0)
        {
            long_name_pattern = "N-glycolyl-";
            cond_name_pattern = "NGc";
            pattern = "C-N-C=OCH2OH";
            AddModificationRuleOneInfo(key, pattern, mono, long_name_pattern, cond_name_pattern, head, tail, minus_one, in_bracket);
        }
        if(value.compare("xC-N-SO3") == 0)
        {
            long_name_pattern = "N-sulfo-";
            cond_name_pattern = "NS";
            pattern = "C-N-SO3";
            AddModificationRuleOneInfo(key, pattern, mono, long_name_pattern, cond_name_pattern, head, tail, minus_one, in_bracket);
        }
        if(value.compare("xC-N-PO3") == 0)
        {
            long_name_pattern = "N-phospho-";
            cond_name_pattern = "NP";
            pattern = "C-N-PO3";
            AddModificationRuleOneInfo(key, pattern, mono, long_name_pattern, cond_name_pattern, head, tail, minus_one, in_bracket);
        }
        if(value.compare("xC-N-CH3") == 0)
        {
            long_name_pattern = "N-methyl-";
            cond_name_pattern = "NMe";
            pattern = "C-N-CH3";
            AddModificationRuleOneInfo(key, pattern, mono, long_name_pattern, cond_name_pattern, head, tail, minus_one, in_bracket);
        }
        if(value.compare("xC-O-C=OCH3") == 0)
        {
            long_name_pattern = "-acetyl-";
            cond_name_pattern = "Ac";
            pattern = "C-O-C=OCH3";
            AddDerivativeRuleInfo(key, pattern, mono, long_name_pattern, cond_name_pattern, head, minus_one, in_bracket);
        }
        if(value.compare("xC-O-C=OCH2OH") == 0)
        {
            long_name_pattern = "-glycolyl-";
            cond_name_pattern = "Gc";
            pattern = "C-O-C=OCH2OH";
            AddDerivativeRuleInfo(key, pattern, mono, long_name_pattern, cond_name_pattern, head, minus_one, in_bracket);
        }
        if(value.compare("xC-O-SO3") == 0)
        {
            long_name_pattern = "-sulfo-";
            cond_name_pattern = "S";
            pattern = "C-O-SO3";
            AddDerivativeRuleInfo(key, pattern, mono, long_name_pattern, cond_name_pattern, head, minus_one, in_bracket);
        }
        if(value.compare("xC-O-PO3") == 0)
        {
            long_name_pattern = "-phospho-";
            cond_name_pattern = "P";
            pattern = "C-O-PO3";
            AddDerivativeRuleInfo(key, pattern, mono, long_name_pattern, cond_name_pattern, head, minus_one, in_bracket);
        }
        if(value.compare("xC-O-CH3") == 0)
        {
            long_name_pattern = "-methyl-";
            cond_name_pattern = "Me";
            pattern = "C-O-CH3";
            AddDerivativeRuleInfo(key, pattern, mono, long_name_pattern, cond_name_pattern, head, minus_one, in_bracket);
        }
        if(value.compare("xC-(O,OH)") == 0)
        {
            long_name_pattern_at_minus_one = "-ulosonic acid";
            long_name_pattern_at_plus_one = "-uronic acid";
            cond_name_pattern = "AH";
            pattern = "C-(O,OH)";
            AddModificationRuleTwoInfo(key, pattern, mono, long_name_pattern_at_minus_one, long_name_pattern_at_plus_one, cond_name_pattern, tail, minus_one, in_bracket);
        }
        if(value.compare("xC-(O,O)") == 0)
        {
            long_name_pattern_at_minus_one = "-ulosonate";
            long_name_pattern_at_plus_one = "-uronate";
            cond_name_pattern = "A";
            pattern = "C-(O,O)";
            AddModificationRuleTwoInfo(key, pattern, mono, long_name_pattern_at_minus_one, long_name_pattern_at_plus_one, cond_name_pattern, tail, minus_one, in_bracket);
        }
    }
    if(in_bracket.str().size() != 0)
    {
        stringstream short_name;
        string sn;
        if(mono->sugar_name_.monosaccharide_short_name_.compare("") != 0)
            sn = mono->sugar_name_.monosaccharide_short_name_;
        else
            sn = mono->sugar_name_.monosaccharide_stereochemistry_short_name_;

        ///moving a, b or x to after the bracket: short-name + [...] + a/b/x and removing ", " from the end of bracket stream
        int condensed_name_size = sn.size();
        string condensed_name = sn;
        string new_name_part1 = condensed_name.substr(0, (condensed_name_size - 1));///short_name
        char new_name_part2 = condensed_name.at(condensed_name_size - 1);///a/b/x
        short_name << new_name_part1 << "[" << in_bracket.str().substr(0, in_bracket.str().size() - 1) << "]" << new_name_part2;

        mono->sugar_name_.monosaccharide_short_name_ = short_name.str();
    }
    else if(mono->sugar_name_.monosaccharide_stereochemistry_short_name_.compare("") != 0 && mono->sugar_name_.monosaccharide_short_name_.compare("") == 0)
    {
        mono->sugar_name_.monosaccharide_short_name_ = mono->sugar_name_.monosaccharide_stereochemistry_short_name_;
    }
    stringstream long_name;
    if(mono->sugar_name_.monosaccharide_stereochemistry_name_.compare("") != 0)
    {
        long_name << head.str() << mono->sugar_name_.monosaccharide_stereochemistry_name_ << tail.str();
        mono->sugar_name_.monosaccharide_name_ = long_name.str();
    }
}
void Assembly::AddModificationRuleOneInfo(string key, string pattern, Monosaccharide* mono, string long_name_pattern, string cond_name_pattern, stringstream& head,
                                          stringstream& tail, bool minus_one, stringstream& in_bracket)
{
    stringstream ss;
    ss << pattern;
    if(key.compare("a") == 0)
    {
        ss << " is at warning position: anomeric";
        gmml::log(__LINE__, __FILE__,  gmml::WAR, ss.str());
        Note* der_mod_note = new Note();
        der_mod_note->type_ = Glycan::WARNING;
        der_mod_note->category_ = Glycan::DER_MOD;
        stringstream note;
        note << mono->sugar_name_.monosaccharide_short_name_ << ": " << ss.str();
        der_mod_note->description_ = note.str();
        this->AddNote(der_mod_note);
        cout << ss.str() << endl;
    }
    else if(key.compare("2") == 0 && mono->sugar_name_.ring_type_.compare("P") == 0 &&
            find(mono->chemical_code_->right_down_.begin(), mono->chemical_code_->right_down_.end(), "-1") == mono->chemical_code_->right_down_.end() &&
            find(mono->chemical_code_->right_up_.begin(), mono->chemical_code_->right_up_.end(), "-1") == mono->chemical_code_->right_up_.end() &&
            mono->sugar_name_.monosaccharide_stereochemistry_name_.compare("") != 0)
    {
        if(long_name_pattern.compare("-osamine") == 0)
            tail << long_name_pattern;
        else
            head << long_name_pattern;
        stringstream short_name;
        if(mono->sugar_name_.monosaccharide_stereochemistry_short_name_.compare("") != 0)
        {
            ///moving a, b or x to after the N expression: short-name + Condensed name pattern + a/b/x
            int stereo_condensed_name_size = mono->sugar_name_.monosaccharide_stereochemistry_short_name_.size();
            string stereo_condensed_name = mono->sugar_name_.monosaccharide_stereochemistry_short_name_;
            string new_name_part1 = stereo_condensed_name.substr(0, (stereo_condensed_name_size - 1));///short_name
            char new_name_part2 = stereo_condensed_name.at(stereo_condensed_name_size - 1);///a/b/x
            short_name << new_name_part1 << cond_name_pattern << new_name_part2;

            mono->sugar_name_.monosaccharide_short_name_ = short_name.str();
        }
    }
    else if(mono->sugar_name_.ring_type_.compare("F") == 0 && key.compare("4") == 0)
    {
        if(!minus_one)
            ss << " is at error position: 4";
        else
            ss << " is at error position: 5";
        gmml::log(__LINE__, __FILE__,  gmml::ERR, ss.str());
        Note* der_mod_note = new Note();
        der_mod_note->type_ = Glycan::ERROR;
        der_mod_note->category_ = Glycan::DER_MOD;
        stringstream note;
        note << mono->sugar_name_.monosaccharide_short_name_ << ": " << ss.str();
        der_mod_note->description_ = note.str();
        this->AddNote(der_mod_note);
        cout << ss.str() << endl;
    }
    else if(mono->sugar_name_.ring_type_.compare("P") == 0 && key.compare("5") == 0)
    {
        if(!minus_one)
            ss << " is at error position: 5";
        else
            ss << " is at error position: 6";
        gmml::log(__LINE__, __FILE__,  gmml::ERR, ss.str());
        Note* der_mod_note = new Note();
        der_mod_note->type_ = Glycan::ERROR;
        der_mod_note->category_ = Glycan::DER_MOD;
        stringstream note;
        note << mono->sugar_name_.monosaccharide_short_name_ << ": " << ss.str();
        der_mod_note->description_ = note.str();
        this->AddNote(der_mod_note);
        cout << ss.str() << endl;
    }
    else
    {
        if(!minus_one)
        {
            if( key.compare("+1") == 0 || key.compare("+2") == 0 || key.compare("+3") == 0)
                in_bracket << mono->cycle_atoms_.size() - 1 + ConvertString<int>(key) << cond_name_pattern << ",";
            else
                in_bracket << key << cond_name_pattern << ",";
        }
        else
        {
            if(key.compare("-1") == 0)
                in_bracket << "1" << cond_name_pattern << ",";
            else if( key.compare("+1") == 0 || key.compare("+2") == 0 || key.compare("+3") == 0)
                in_bracket << mono->cycle_atoms_.size() + ConvertString<int>(key) << cond_name_pattern << ",";
            else
                in_bracket << ConvertString<int>(key) + 1 << cond_name_pattern << ",";
        }
    }
}
void Assembly::AddDerivativeRuleInfo(string key, string pattern, Monosaccharide *mono, string long_name_pattern, string cond_name_pattern, stringstream &head,
                                     bool minus_one, stringstream &in_bracket)
{
    stringstream ss;
    ss << pattern;
    if(mono->sugar_name_.monosaccharide_stereochemistry_name_.compare("") != 0)
    {
        if(!minus_one)
        {
            if(key.compare("+1") == 0 || key.compare("+2") == 0 || key.compare("+3") == 0)
                head << mono->cycle_atoms_.size() - 1 + ConvertString<int>(key) << long_name_pattern;
            //            else if(key.compare("a") == 0)
            //                head << "2" << long_name_pattern;
            else if(key.compare("a") != 0)
                head << ConvertString<int>(key) << long_name_pattern;
        }
        else
        {
            if(key.compare("-1") == 0)
                head << "1" << long_name_pattern;
            else if(key.compare("+1") == 0 || key.compare("+2") == 0 || key.compare("+3") == 0)
                head << mono->cycle_atoms_.size() + ConvertString<int>(key) << long_name_pattern;
            //            else if(key.compare("a") == 0)
            //                head << "2" << long_name_pattern;
            else if(key.compare("a") != 0)
                head << ConvertString<int>(key) + 1 << long_name_pattern;
        }
    }
    if(mono->sugar_name_.monosaccharide_stereochemistry_short_name_.compare("") != 0)
    {
        if(mono->sugar_name_.ring_type_.compare("F") == 0 && key.compare("4") == 0)
        {
            if(!minus_one)
                ss << " is at error position: 4";
            else
                ss << " is at error position: 5";
            gmml::log(__LINE__, __FILE__,  gmml::ERR, ss.str());
            Note* der_mod_note = new Note();
            der_mod_note->type_ = Glycan::ERROR;
            der_mod_note->category_ = Glycan::DER_MOD;
            stringstream note;
            note << mono->sugar_name_.monosaccharide_short_name_ << ": " << ss.str();
            der_mod_note->description_ = note.str();
            this->AddNote(der_mod_note);
            cout << ss.str() << endl;
        }
        else if(mono->sugar_name_.ring_type_.compare("P") == 0 && key.compare("5") == 0)
        {
            if(!minus_one)
                ss << " is at error position: 5";
            else
                ss << " is at error position: 6";
            gmml::log(__LINE__, __FILE__,  gmml::ERR, ss.str());
            Note* der_mod_note = new Note();
            der_mod_note->type_ = Glycan::ERROR;
            der_mod_note->category_ = Glycan::DER_MOD;
            stringstream note;
            note << mono->sugar_name_.monosaccharide_short_name_ << ": " << ss.str();
            der_mod_note->description_ = note.str();
            this->AddNote(der_mod_note);
            cout << ss.str() << endl;
        }
        else
        {
            if(!minus_one)
            {
                if(key.compare("+1") == 0 || key.compare("+2") == 0 || key.compare("+3") == 0)
                    in_bracket << mono->cycle_atoms_.size() - 1 + ConvertString<int>(key) << cond_name_pattern << ",";
                //                else if(key.compare("a") == 0)
                //                    in_bracket << "2" << cond_name_pattern << ",";
                else if(key.compare("a") != 0)
                    in_bracket << ConvertString<int>(key) << cond_name_pattern << ",";
            }
            else
            {
                if(key.compare("-1") == 0)
                    in_bracket << "1" << cond_name_pattern << ",";
                else if( key.compare("+1") == 0 || key.compare("+2") == 0 || key.compare("+3") == 0)
                    in_bracket << mono->cycle_atoms_.size() + ConvertString<int>(key) << cond_name_pattern << ",";
                //                else if(key.compare("a") == 0)
                //                    in_bracket << "2" << cond_name_pattern << ",";
                else if(key.compare("a") != 0)
                    in_bracket << ConvertString<int>(key) + 1 << cond_name_pattern << ",";
            }
        }
    }
}
void Assembly::AddModificationRuleTwoInfo(string key, string pattern, Monosaccharide *mono, string long_name_pattern_at_minus_one, string long_name_pattern_at_plus_one,
                                          string cond_name_pattern, stringstream &tail, bool minus_one, stringstream &in_bracket)
{
    Note* der_mod_note = new Note();
    stringstream ss;
    ss << pattern;
    if((key.compare("-1") == 0 || key.compare("+2") == 0 || key.compare("+3") == 0) && mono->sugar_name_.monosaccharide_stereochemistry_name_.compare("") != 0)
    {
        if(key.compare("-1") == 0)
        {
            ss << " is at warning position: " << ConvertString<int>(key) + 1;
            der_mod_note->type_ = Glycan::WARNING;
        }
        if(!minus_one)
        {
            ss << " is at error position: " << key;
            der_mod_note->type_ = Glycan::ERROR;
        }
        else
        {
            ss << " is at error position: " << ConvertString<int>(key) + 1;
            der_mod_note->type_ = Glycan::ERROR;
        }
        cout << ss.str() << endl;
        gmml::log(__LINE__, __FILE__,  gmml::ERR, ss.str());
        tail << long_name_pattern_at_minus_one;

        if(mono->sugar_name_.monosaccharide_stereochemistry_short_name_.compare("") != 0)
        {
            if(!minus_one)
            {
                if( key.compare("+2") == 0 || key.compare("+3") == 0)
                    in_bracket << mono->cycle_atoms_.size() - 1 + ConvertString<int>(key) << cond_name_pattern << ",";
            }
            else
            {
                if(key.compare("-1") == 0)
                    in_bracket << "1" << cond_name_pattern << ",";
                else if( key.compare("+2") == 0 || key.compare("+3") == 0)
                    in_bracket << mono->cycle_atoms_.size() + ConvertString<int>(key) << cond_name_pattern << ",";
            }
        }
    }
    else if(key.compare("+1") == 0 && mono->sugar_name_.monosaccharide_stereochemistry_name_.compare("") != 0)
    {
        tail << long_name_pattern_at_plus_one;
        stringstream short_name;
        if(mono->sugar_name_.monosaccharide_stereochemistry_short_name_.compare("") != 0)
        {
            ///moving a, b or x to after the AH expression: short-name + AH + a/b/x
            int stereo_condensed_name_size = mono->sugar_name_.monosaccharide_stereochemistry_short_name_.size();
            string stereo_condensed_name = mono->sugar_name_.monosaccharide_stereochemistry_short_name_;
            string new_name_part1 = stereo_condensed_name.substr(0, (stereo_condensed_name_size - 1));///short_name
            char new_name_part2 = stereo_condensed_name.at(stereo_condensed_name_size - 1);///a/b/x
            short_name << new_name_part1 << cond_name_pattern << new_name_part2;

            mono->sugar_name_.monosaccharide_short_name_ = short_name.str();
        }
    }
    else
    {
        if(!minus_one)
        {
            if(key.compare("a") != 0)
                ss << " is at warning position: " << key;
            else
                ss << " is at warning position: 1";
        }
        else
        {
            if(key.compare("a") != 0)
                ss << " is at warning position: " << ConvertString<int>(key) + 1;
            else
                ss << " is at warning position: 2";
        }
        der_mod_note->type_ = Glycan::WARNING;
    }
    der_mod_note->category_ = Glycan::DER_MOD;
    stringstream note;
    note << mono->sugar_name_.monosaccharide_short_name_ << ": " << ss.str();
    der_mod_note->description_ = note.str();
    this->AddNote(der_mod_note);
    cout << ss.str() << endl;
    gmml::log(__LINE__, __FILE__,  gmml::WAR, ss.str());
}

void Assembly::UpdateComplexSugarChemicalCode(Monosaccharide *mono)
{
    if(mono->side_atoms_.at(0).at(0) != NULL)
    {
        if(mono->derivatives_map_["-1"].compare("xC-(O,O)") == 0)
        {
            vector<string>::iterator index_it;
            if((index_it = find(mono->chemical_code_->right_up_.begin(), mono->chemical_code_->right_up_.end(), "-1")) != mono->chemical_code_->right_up_.end())
                (*index_it) = "-1A";
            else if((index_it = find(mono->chemical_code_->right_down_.begin(), mono->chemical_code_->right_down_.end(), "-1")) != mono->chemical_code_->right_down_.end())
                (*index_it) = "-1A";
        }
    }
    if(mono->cycle_atoms_.at(2) != NULL)
    {
        if(mono->derivatives_map_["-1"].compare("xC-O-C=OCH3") == 0)
        {
            vector<string>::iterator index_it;
            if((index_it = find(mono->chemical_code_->left_up_.begin(), mono->chemical_code_->left_up_.end(), "3")) != mono->chemical_code_->left_up_.end())
                (*index_it) = "3Ac";
            else if((index_it = find(mono->chemical_code_->left_down_.begin(), mono->chemical_code_->left_down_.end(), "3")) != mono->chemical_code_->left_down_.end())
                (*index_it) = "3Ac";
        }
    }
    if(mono->cycle_atoms_.at(3) != NULL)
    {
        if(mono->derivatives_map_["4"].compare("xCH-N") == 0)
        {
            vector<string>::iterator index_it;
            if((index_it = find(mono->chemical_code_->left_up_.begin(), mono->chemical_code_->left_up_.end(), "4")) != mono->chemical_code_->left_up_.end())
                (*index_it) = "4N";
            else if((index_it = find(mono->chemical_code_->left_down_.begin(), mono->chemical_code_->left_down_.end(), "4")) != mono->chemical_code_->left_down_.end())
                (*index_it) = "4N";
        }
        if(mono->derivatives_map_["4"].compare("xC-N-C=OCH3") == 0)
        {
            vector<string>::iterator index_it;
            if((index_it = find(mono->chemical_code_->left_up_.begin(), mono->chemical_code_->left_up_.end(), "4")) != mono->chemical_code_->left_up_.end())
                (*index_it) = "4NAc";
            else if((index_it = find(mono->chemical_code_->left_down_.begin(), mono->chemical_code_->left_down_.end(), "4")) != mono->chemical_code_->left_down_.end())
                (*index_it) = "4NAc";
        }
        if(mono->derivatives_map_["4"].compare("xC-N-C=OCH2OH") == 0)
        {
            vector<string>::iterator index_it;
            if((index_it = find(mono->chemical_code_->left_up_.begin(), mono->chemical_code_->left_up_.end(), "4")) != mono->chemical_code_->left_up_.end())
                (*index_it) = "4NGc";
            else if((index_it = find(mono->chemical_code_->left_down_.begin(), mono->chemical_code_->left_down_.end(), "4")) != mono->chemical_code_->left_down_.end())
                (*index_it) = "4NGc";
        }
    }
}

vector<Oligosaccharide*> Assembly::ExtractOligosaccharides(vector<Monosaccharide*> monos, ResidueNameMap dataset_residue_names,
                                                           int& number_of_covalent_links, int& number_of_probable_non_covalent_complexes)
{
    string terminal_residue_name = "";
    ResidueNameMap common_terminal_residues = gmml::InitializeCommonTerminalResidueMap();
    map<Monosaccharide*, vector<Monosaccharide*> > monos_table = map<Monosaccharide*, vector<Monosaccharide*> >();
    map<Monosaccharide*, vector<string> > monos_table_linkages = map<Monosaccharide*, vector<string> >();

    ///Iterating on list of monos to check if there is a connection to another mono in the list
    for(vector<Monosaccharide*>::iterator it = monos.begin(); it != monos.end(); it++)
    {
        Monosaccharide* mono1 = (*it);

        monos_table[mono1] = vector<Monosaccharide*>();
        monos_table_linkages[mono1] = vector<string>();

        for(vector<AtomVector>::iterator it1 = mono1->side_atoms_.begin(); it1 != mono1->side_atoms_.end(); it1++) ///iterate on side atoms
        {
            int index = distance(mono1->side_atoms_.begin(), it1);
            AtomVector sides = (*it1);
            map<Atom*, Atom*> target_parent_map = map<Atom*, Atom*>(); /// A map of target atom to it's parent atom. Target atom is a non ring oxygen or nitrogen

            if(it1 == mono1->side_atoms_.begin())///side atoms of anomeric
            {
                if(sides.at(1) != NULL)
                    target_parent_map[sides.at(1)] = mono1->cycle_atoms_.at(0);
            }
            else if(it1 == mono1->side_atoms_.end() - 1) ///side atoms of last carbon of the ring
            {
                for(AtomVector::iterator last_c_side_it = sides.begin(); last_c_side_it != sides.end(); last_c_side_it++)
                {
                    Atom* side_of_last_carbon = (*last_c_side_it);
                    if(side_of_last_carbon != NULL)
                    {
                        AtomVector last_c_side_neighbors = side_of_last_carbon->GetNode()->GetNodeNeighbors();
                        for(AtomVector::iterator it2 = last_c_side_neighbors.begin(); it2 != last_c_side_neighbors.end(); it2++)
                        {
                            if((*it2)->GetId().at(0) == 'O' || (*it2)->GetId().at(0) == 'N')
                            {
                                target_parent_map[(*it2)] = side_of_last_carbon;
                                break;
                            }
                        }
                    }
                }
            }
            else
            {
                if(sides.at(1) != NULL)
                    target_parent_map[sides.at(1)] = mono1->cycle_atoms_.at(index);///index 1 of each side is for non-carbon side atoms in the vector<AtomVector> structure
            }
            ///Examine neighbors of each target atom to check if they can be found in other monos side/ring atoms
            for(map<Atom*, Atom*>::iterator map_it = target_parent_map.begin(); map_it != target_parent_map.end(); map_it++)
            {
                bool found_in_other_mono = false;
                Atom* target = (*map_it).first;
                Atom* target_parent = (*map_it).second;
                AtomVector t_neighbors = target->GetNode()->GetNodeNeighbors();
                for(AtomVector::iterator it2 = t_neighbors.begin(); it2 != t_neighbors.end(); it2++)
                {
                    Atom* t_neighbor = (*it2);
                    if(t_neighbor->GetId().compare(target_parent->GetId()) != 0)///making sure neighbor is not the parent of target atom
                    {
                        for(vector<Monosaccharide*>::iterator it3 = monos.begin(); it3 != monos.end(); it3++)
                        {
                            if(it3 != it)///Cheking monos other than the current mono
                            {
                                Monosaccharide* mono2 = (*it3);
                                AtomVector mono2_sides = mono2->side_atoms_.at(mono2->side_atoms_.size() - 1); ///side of last ring carbon

                                bool found_in_side = false;
                                for(AtomVector::iterator mono2_last_c_side_it = mono2_sides.begin(); mono2_last_c_side_it != mono2_sides.end(); mono2_last_c_side_it++)
                                {
                                    Atom* mono2_last_c_side = (*mono2_last_c_side_it);
                                    if(mono2_last_c_side != NULL)
                                    {
                                        if(t_neighbor->GetId().compare(mono2_last_c_side->GetId()) == 0) ///target atom has been attached to another cycle's side atom
                                            found_in_side = true;
                                    }
                                }
                                if(found_in_side || mono2->cycle_atoms_str_.find(t_neighbor->GetId()) != string::npos) //if target's neighbor found in another mono's side or ring atoms
                                {
                                    found_in_other_mono = true;
                                    monos_table[mono1].push_back(mono2);

                                    string mono1_carbon = target_parent->GetId();
                                    string mono1_name = "";
                                    string mono2_carbon = t_neighbor->GetId();
                                    string mono2_name = "";
                                    if(mono1->sugar_name_.monosaccharide_short_name_.compare("") != 0)
                                        mono1_name = mono1->sugar_name_.monosaccharide_short_name_;
                                    else
                                        mono1_name = mono1->sugar_name_.monosaccharide_stereochemistry_short_name_;

                                    if(mono2->sugar_name_.monosaccharide_short_name_.compare("") != 0)
                                        mono2_name = mono2->sugar_name_.monosaccharide_short_name_;
                                    else
                                        mono2_name = mono2->sugar_name_.monosaccharide_stereochemistry_short_name_;

                                    stringstream linkage;
                                    linkage << mono1_carbon << "-" << target->GetId() << "-" << mono2_carbon;
                                    monos_table_linkages[mono1].push_back(linkage.str());
                                    break;
                                }
                            }
                        }
                    }
                }
                if(found_in_other_mono)
                    break;
            }
        }
    }
    vector<int> visited_monos = vector<int>();
    vector<Oligosaccharide*> oligosaccharides = vector<Oligosaccharide*>();

    vector<string> checked_linkages = vector<string>();
    for(map<Monosaccharide*, vector<Monosaccharide*> >::iterator it = monos_table.begin(); it != monos_table.end(); it++)
    {
        Monosaccharide* key = (*it).first;
        vector<Monosaccharide*> values = (*it).second;

        vector<string> visited_linkages = vector<string>();
        if(find(visited_monos.begin(), visited_monos.end(), key->mono_id) == visited_monos.end())///if the mono is not visited
        {
            bool isRoot = false;
            stringstream anomeric_linkage;
            anomeric_linkage << key->cycle_atoms_.at(0)->GetId() << "-";
            vector<string> mono_linkages = monos_table_linkages[key];
            AtomVector terminal_atoms = AtomVector();
            if(values.size() == 0) ///mono is not attached to any other mono
            {
                Atom* anomeric_o = NULL;
                if(key->side_atoms_.at(0).at(1) != NULL)
                    anomeric_o = key->side_atoms_.at(0).at(1);
                if(anomeric_o != NULL)
                {
                    if(dataset_residue_names.find(anomeric_o->GetResidue()->GetName()) != dataset_residue_names.end() ||
                            common_terminal_residues.find(anomeric_o->GetResidue()->GetName()) != common_terminal_residues.end())///check if there is any terminal
                    {
                        terminal_residue_name = anomeric_o->GetResidue()->GetName();
                    }
                    else
                    {
                        terminal_residue_name = CheckTerminals(anomeric_o, terminal_atoms);
                    }
                }
                isRoot = true;
            }
            else if (values.size() == 1 && ///mono is attached to one other mono
                     (monos_table_linkages[values.at(0)].size() == 1))///the other mono is only attached to this mono
            {
                ///CHECKING LINKAGE ISSUES, e.g. C1-O3-C4 is an issue
                CheckLinkageNote(key, values.at(0), mono_linkages.at(0), checked_linkages);
                stringstream other_mono_anomeric_linkage_as_right_side;
                other_mono_anomeric_linkage_as_right_side << "-" << values.at(0)->cycle_atoms_.at(0)->GetId();///atom id on the right side of the linkage c-o-c

                Atom* anomeric_o = NULL;
                Atom* o_neighbor_1 = NULL;
                Atom* o_neighbor_2 = NULL;
                AtomVector o_neighbors = AtomVector();
                if(key->side_atoms_.at(0).at(1) != NULL)///Getting the information of anomeric oxygen's neighbors is needed for choosing the root
                {
                    anomeric_o = key->side_atoms_.at(0).at(1);
                    o_neighbors = anomeric_o->GetNode()->GetNodeNeighbors();
                    if(o_neighbors.size() > 1)
                    {
                        o_neighbor_1 = o_neighbors.at(0);
                        o_neighbor_2 = o_neighbors.at(1);
                    }
                }
                if(anomeric_o != NULL)
                {
                    ///RULE1: anomeric to anomeric linkage
                    if(((mono_linkages.at(0)).find(anomeric_linkage.str()) != string::npos) && ///this mono is attached to other mono through anomeric
                            (mono_linkages.at(0).find(other_mono_anomeric_linkage_as_right_side.str()) != string::npos))///the other mono is only attached to this mono through anomeric
                        isRoot = true;
                    ///RULE2: Directed graph
                    else if(((mono_linkages.at(0)).find(anomeric_linkage.str()) == string::npos) && ///this mono is not attached to other mono through anomeric
                            (mono_linkages.at(0).find(other_mono_anomeric_linkage_as_right_side.str()) != string::npos)) ///the other mono is attached to this mono through anomeric
                    {
                        isRoot = true;
                        terminal_residue_name = CheckTerminals(anomeric_o, terminal_atoms);
                    }
                    ///RULE3: Terminal
                    else if(o_neighbors.size() == 1) ///anomeric oxygen is not attached to anything else except the carbon of the ring
                    {
                        isRoot = true;
                        terminal_residue_name = CheckTerminals(anomeric_o, terminal_atoms);
                    }
                    else if(o_neighbors.size() == 2 && (((o_neighbor_1->GetDescription().find("Het;") != string::npos) && (o_neighbor_2->GetDescription().find("Het;") == string::npos)) ||
                                                        ((o_neighbor_2->GetDescription().find("Het;") != string::npos) && (o_neighbor_1->GetDescription().find("Het;") == string::npos))) )
                    {
                        ///anomeric oxygen is attached to protein
                        isRoot = true;
                        terminal_residue_name = CheckTerminals(anomeric_o, terminal_atoms);
                        number_of_covalent_links++;
                        if(terminal_residue_name.compare("NLN") != 0 && terminal_residue_name.compare("OLS") != 0 && terminal_residue_name.compare("OLT") != 0)
                        {
                            stringstream ss;
                            ss << "Root anomeric atom is attached to a non-standard " << terminal_residue_name << " protein residue!";
                            gmml::log(__LINE__, __FILE__, gmml::WAR, ss.str());
                            cout << ss.str() << endl;
                            terminal_residue_name = "";
                        }
                    }
                    else if(dataset_residue_names.find(anomeric_o->GetResidue()->GetName()) != dataset_residue_names.end() ||
                            common_terminal_residues.find(anomeric_o->GetResidue()->GetName()) != common_terminal_residues.end())///mono is attached to a terminal through anomeric oxygen
                    {
                        terminal_residue_name = anomeric_o->GetResidue()->GetName();
                        isRoot = true;
                    }
                    else if((terminal_residue_name = CheckTerminals(anomeric_o, terminal_atoms)).compare("") != 0)
                        isRoot = true;
                }
                ///RULE2: Directed graph
                else if((mono_linkages.at(0).find(other_mono_anomeric_linkage_as_right_side.str()) != string::npos)) ///this mono doesn't have anomeric oxygen and the other mono is attached to this mono through anomeric
                {
                    isRoot = true;
                    terminal_residue_name = CheckTerminals(anomeric_o, terminal_atoms);
                    number_of_probable_non_covalent_complexes++;
                }
            }
            else
            {
                Atom* anomeric_o = NULL;
                Atom* o_neighbor_1 = NULL;
                Atom* o_neighbor_2 = NULL;
                AtomVector o_neighbors = AtomVector();
                if(key->side_atoms_.at(0).at(1) != NULL)///Getting the information of anomeric oxygen's neighbors is needed for choosing the root
                {
                    anomeric_o = key->side_atoms_.at(0).at(1);
                    o_neighbors = anomeric_o->GetNode()->GetNodeNeighbors();
                    if(o_neighbors.size() > 1)
                    {
                        o_neighbor_1 = o_neighbors.at(0);
                        o_neighbor_2 = o_neighbors.at(1);
                    }
                }
                if(anomeric_o != NULL)
                {
                    ///RULE1: anomeric to anomeric linkage
                    for(int i = 0; i < values.size(); i++)
                    {
                        CheckLinkageNote(key, values.at(i), mono_linkages.at(i), checked_linkages);
                        stringstream other_mono_anomeric_linkage_as_right_side;
                        other_mono_anomeric_linkage_as_right_side << "-" << values.at(i)->cycle_atoms_.at(0)->GetId();///atom id on the right side of the linkage c-o-c

                        if(((mono_linkages.at(i)).find(anomeric_linkage.str()) != string::npos) && ///this mono is attached to another mono through anomeric
                                (mono_linkages.at(i).find(other_mono_anomeric_linkage_as_right_side.str()) != string::npos))///the other mono is attached to this mono through anomeric
                        {
                            isRoot = true;
                            break;
                        }
                    }
                    if(!isRoot) ///RULE2: Directed graph
                    {
                        for(int i = 0; i < values.size(); i++)
                        {
                            CheckLinkageNote(key, values.at(i), mono_linkages.at(i), checked_linkages);
                            stringstream other_mono_anomeric_linkage_as_right_side;
                            other_mono_anomeric_linkage_as_right_side << "-" << values.at(i)->cycle_atoms_.at(0)->GetId();///atom id on the right side of the linkage c-o-c
                            if(((mono_linkages.at(i)).find(anomeric_linkage.str()) != string::npos)) ///this mono is attached to other mono through anomeric
                            {
                                isRoot = false;
                                break;
                            }
                            else if((mono_linkages.at(i).find(other_mono_anomeric_linkage_as_right_side.str()) != string::npos)) ///the other mono is attached to this mono through anomeric
                            {
                                isRoot = true;
                                terminal_residue_name = CheckTerminals(anomeric_o, terminal_atoms);
                            }
                        }
                    }
                    else if(!isRoot)///RULE3: Terminal
                    {
                        if(o_neighbors.size() == 1) ///anomeric oxygen is not attached to anything else, except the carbon of the ring
                        {
                            isRoot = true;
                            terminal_residue_name = CheckTerminals(anomeric_o, terminal_atoms);
                        }
                        else if(o_neighbors.size() == 2 && (((o_neighbor_1->GetDescription().find("Het;") != string::npos) && (o_neighbor_2->GetDescription().find("Het;") == string::npos)) ||
                                                            ((o_neighbor_2->GetDescription().find("Het;") != string::npos) && (o_neighbor_1->GetDescription().find("Het;") == string::npos))) )
                        {
                            ///anomeric oxygen is attached to protein
                            isRoot = true;
                            terminal_residue_name = CheckTerminals(anomeric_o, terminal_atoms);
                            number_of_covalent_links++;
                            if(terminal_residue_name.compare("NLN") != 0 && terminal_residue_name.compare("OLS") != 0 && terminal_residue_name.compare("OLT") != 0)
                            {
                                stringstream ss;
                                ss << "Root anomeric atom is attached to a non-standard " << terminal_residue_name << " protein residue!";
                                gmml::log(__LINE__, __FILE__, gmml::WAR, ss.str());
                                cout << ss.str() << endl;
                                terminal_residue_name = "";
                            }
                        }
                        else if(dataset_residue_names.find(anomeric_o->GetResidue()->GetName()) != dataset_residue_names.end() ||
                                common_terminal_residues.find(anomeric_o->GetResidue()->GetName()) != common_terminal_residues.end())///mono is attached to a terminal through anomeric oxygen
                        {
                            terminal_residue_name = anomeric_o->GetResidue()->GetName();
                            isRoot = true;
                        }
                        else if((terminal_residue_name = CheckTerminals(anomeric_o, terminal_atoms)).compare("") != 0)
                            isRoot = true;
                    }
                }
                //this mono doesn't have anomeric oxygen
                else ///RULE2: Directed graph
                {
                    for(int i = 0; i < values.size(); i++)
                    {
                        CheckLinkageNote(key, values.at(i), mono_linkages.at(i), checked_linkages);
                        vector<string> other_mono_linkage = monos_table_linkages[values.at(i)];
                        stringstream other_mono_anomeric_linkage;
                        other_mono_anomeric_linkage << values.at(i)->cycle_atoms_.at(0)->GetId() << "-";///atom id on the left side of the linkage c-o-c
                        if((other_mono_linkage.at(0).find(other_mono_anomeric_linkage.str()) != string::npos)) ///the other mono is attached to this mono through anomeric
                        {
                            isRoot = true;
                            terminal_residue_name = CheckTerminals(anomeric_o, terminal_atoms);
                            number_of_probable_non_covalent_complexes++;
                            break;
                        }
                    }
                }
            }
            if(isRoot)
            {
                Oligosaccharide* oligo = new Oligosaccharide();
                BuildOligosaccharideTreeStructure(key, values, oligo, visited_monos, monos_table, monos_table_linkages, visited_linkages);
                oligo->terminal_ = terminal_residue_name;
                oligosaccharides.push_back(oligo);
            }
        }
    }

    for(map<Monosaccharide*, vector<Monosaccharide*> >::iterator it = monos_table.begin(); it != monos_table.end(); it++)
    {
        Monosaccharide* key = (*it).first;
        vector<Monosaccharide*> values = (*it).second;
        if(values.size() > 1)
        {
            vector<string> visited_linkages = vector<string>();
            if(find(visited_monos.begin(), visited_monos.end(), key->mono_id) == visited_monos.end())///if the mono is not visited
            {
                vector<string> mono_linkages = monos_table_linkages[key];
                stringstream anomeric_linkage;
                anomeric_linkage << key->cycle_atoms_.at(0)->GetId() << "-";

                for(vector<string>::iterator it1 = mono_linkages.begin(); it1 != mono_linkages.end(); it1++)
                {
                    if((*it1).find(anomeric_linkage.str()) != string::npos)///mono is attached to another mono through anomeric
                    {
                        Oligosaccharide* oligo = new Oligosaccharide();
                        BuildOligosaccharideTreeStructure(key, values, oligo, visited_monos, monos_table, monos_table_linkages, visited_linkages);
                        oligosaccharides.push_back(oligo);
                        break;
                    }
                }
            }
        }
    }
    return oligosaccharides;
}

string Assembly::CheckOMETerminal(Atom* target, AtomVector& terminal_atoms)
{
    terminal_atoms = AtomVector();
    AtomVector atoms_1 = AtomVector();
    AtomVector atoms_2 = AtomVector();
    stringstream pattern;
    pattern << "O";
    atoms_1.push_back(target);
    atoms_2.push_back(target);
    Atom* C1 = NULL;
    Atom* C2 = NULL;
    AtomVector o_neighbors = target->GetNode()->GetNodeNeighbors();
    for(AtomVector::iterator it = o_neighbors.begin(); it != o_neighbors.end(); it++)
    {
        Atom* o_neighbor = (*it);

        if(o_neighbor->GetName().at(0) == 'C' && C1 == NULL && C2 == NULL)
            C1 = o_neighbor;
        else if(o_neighbor->GetName().at(0) == 'C' && C1 != NULL && C2 == NULL)
            C2 = o_neighbor;
    }
    stringstream temp;
    if(C1 != NULL)
    {
        temp << pattern.str() << "-C";
        atoms_1.push_back(C1);
        AtomVector c1_neighbors = C1->GetNode()->GetNodeNeighbors();
        for(AtomVector::iterator it = c1_neighbors.begin(); it != c1_neighbors.end(); it++)
        {
            Atom* c1_neighbor = (*it);
            if(c1_neighbor->GetId().compare(target->GetId()) != 0)
            {
                temp << c1_neighbor->GetName().at(0);
                atoms_1.push_back(c1_neighbor);
            }
        }
    }
    if(temp.str().compare("O-C") == 0 || temp.str().compare("O-CHHH") == 0)
    {
        terminal_atoms = atoms_1;
        return "OME";
    }
    if(C2 != NULL)
    {
        pattern << "-C";
        atoms_2.push_back(C2);
        AtomVector c2_neighbors = C2->GetNode()->GetNodeNeighbors();
        for(AtomVector::iterator it = c2_neighbors.begin(); it != c2_neighbors.end(); it++)
        {
            Atom* c2_neighbor = (*it);
            if(c2_neighbor->GetId().compare(target->GetId()) != 0)
            {
                pattern << c2_neighbor->GetName().at(0);
                atoms_2.push_back(c2_neighbor);
            }
        }
        if(pattern.str().compare("O-C") == 0 || pattern.str().compare("O-CHHH") == 0)
        {
            terminal_atoms = atoms_2;
            return "OME";
        }
    }
    return "";
}
string Assembly::CheckROHTerminal(Atom* target, AtomVector& terminal_atoms)
{
    terminal_atoms = AtomVector();
    AtomVector o_neighbors = target->GetNode()->GetNodeNeighbors();
    if(o_neighbors.size() == 1)
    {
        terminal_atoms.push_back(target);
        return "ROH";
    }
    else if (o_neighbors.size() > 1)
    {
        if((o_neighbors.at(0)->GetName().at(0) == 'H' && o_neighbors.at(1)->GetName().at(0) != 'H'))
        {
            terminal_atoms.push_back(target);
            terminal_atoms.push_back(o_neighbors.at(0));
            return "ROH";
        }
        else if(o_neighbors.at(1)->GetName().at(0) == 'H' && o_neighbors.at(0)->GetName().at(0) != 'H')
        {
            terminal_atoms.push_back(target);
            terminal_atoms.push_back(o_neighbors.at(1));
            return "ROH";
        }
    }
    return "";
}
string Assembly::CheckTBTTerminal(Atom *target, AtomVector& terminal_atoms)
{
    terminal_atoms = AtomVector();
    AtomVector atoms_1 = AtomVector();
    AtomVector atoms_2 = AtomVector();
    stringstream pattern;
    pattern << "O";
    atoms_1.push_back(target);
    atoms_2.push_back(target);
    Atom* C1 = NULL;
    Atom* C2 = NULL;
    AtomVector o_neighbors = target->GetNode()->GetNodeNeighbors();
    for(AtomVector::iterator it = o_neighbors.begin(); it != o_neighbors.end(); it++)
    {
        Atom* o_neighbor = (*it);
        if(o_neighbor->GetName().at(0) == 'C' && C1 == NULL && C2 == NULL)
            C1 = o_neighbor;
        else if(o_neighbor->GetName().at(0) == 'C' && C1 != NULL && C2 == NULL)
            C2 = o_neighbor;
    }
    stringstream temp;
    if(C1 != NULL)
    {
        Atom* C1C1 = NULL;
        Atom* C1C2 = NULL;
        Atom* C1C3 = NULL;
        temp << pattern.str() << "-" << "C";
        atoms_1.push_back(C1);
        AtomVector c1_neighbors = C1->GetNode()->GetNodeNeighbors();
        for(AtomVector::iterator it = c1_neighbors.begin(); it != c1_neighbors.end(); it++)
        {
            Atom* c1_neighbor = (*it);
            if(c1_neighbor->GetId().compare(target->GetId()) != 0 && c1_neighbor->GetName().at(0) == 'C' && C1C1 == NULL)
                C1C1 = c1_neighbor;
            else if(c1_neighbor->GetId().compare(target->GetId()) != 0 && c1_neighbor->GetName().at(0) == 'C' && C1C1 != NULL && C1C2 == NULL)
                C1C2 = c1_neighbor;
            else if(c1_neighbor->GetId().compare(target->GetId()) != 0 && c1_neighbor->GetName().at(0) == 'C' && C1C1 != NULL && C1C2 != NULL && C1C3 == NULL)
                C1C3 = c1_neighbor;
        }
        if(C1C1 != NULL && C1C2 != NULL && C1C3 != NULL)
        {
            temp << "C";
            atoms_1.push_back(C1C1);
            AtomVector c1c1_neighbors = C1C1->GetNode()->GetNodeNeighbors();
            for(AtomVector::iterator it = c1c1_neighbors.begin(); it != c1c1_neighbors.end(); it++)
            {
                Atom* c1c1_neighbor = (*it);
                if(c1c1_neighbor->GetId().compare(C1->GetId()) != 0)
                {
                    temp << c1c1_neighbor->GetName().at(0);
                    atoms_1.push_back(c1c1_neighbor);
                }
            }
            temp << "C";
            atoms_1.push_back(C1C2);
            AtomVector c1c2_neighbors = C1C2->GetNode()->GetNodeNeighbors();
            for(AtomVector::iterator it = c1c2_neighbors.begin(); it != c1c2_neighbors.end(); it++)
            {
                Atom* c1c2_neighbor = (*it);
                if(c1c2_neighbor->GetId().compare(C1->GetId()) != 0)
                {
                    temp << c1c2_neighbor->GetName().at(0);
                    atoms_1.push_back(c1c2_neighbor);
                }
            }
            temp << "C";
            atoms_1.push_back(C1C3);
            AtomVector c1c3_neighbors = C1C3->GetNode()->GetNodeNeighbors();
            for(AtomVector::iterator it = c1c3_neighbors.begin(); it != c1c3_neighbors.end(); it++)
            {
                Atom* c1c3_neighbor = (*it);
                if(c1c3_neighbor->GetId().compare(C1->GetId()) != 0)
                {
                    temp << c1c3_neighbor->GetName().at(0);
                    atoms_1.push_back(c1c3_neighbor);
                }
            }
        }
        if(temp.str().compare("O-CCCC") == 0 || temp.str().compare("O-CCHHHCHHHCHHH") == 0 )
        {
            terminal_atoms = atoms_1;
            return "TBT";
        }
    }
    if(C2 != NULL)
    {
        Atom* C2C1 = NULL;
        Atom* C2C2 = NULL;
        Atom* C2C3 = NULL;
        pattern << "-" << "C";
        atoms_2.push_back(C2);
        AtomVector c2_neighbors = C2->GetNode()->GetNodeNeighbors();
        for(AtomVector::iterator it = c2_neighbors.begin(); it != c2_neighbors.end(); it++)
        {
            Atom* c2_neighbor = (*it);
            if(c2_neighbor->GetId().compare(target->GetId()) != 0 && c2_neighbor->GetName().at(0) == 'C' && C2C1 == NULL)
                C2C1 = c2_neighbor;
            else if(c2_neighbor->GetId().compare(target->GetId()) != 0 && c2_neighbor->GetName().at(0) == 'C' && C2C1 != NULL && C2C2 == NULL)
                C2C2 = c2_neighbor;
            else if(c2_neighbor->GetId().compare(target->GetId()) != 0 && c2_neighbor->GetName().at(0) == 'C' && C2C1 != NULL && C2C2 != NULL && C2C3 == NULL)
                C2C3 = c2_neighbor;
        }
        if(C2C1 != NULL && C2C2 != NULL && C2C3 != NULL)
        {
            pattern << "C";
            atoms_2.push_back(C2C1);
            AtomVector c2c1_neighbors = C2C1->GetNode()->GetNodeNeighbors();
            for(AtomVector::iterator it = c2c1_neighbors.begin(); it != c2c1_neighbors.end(); it++)
            {
                Atom* c2c1_neighbor = (*it);
                if(c2c1_neighbor->GetId().compare(C2->GetId()) != 0)
                {
                    pattern << c2c1_neighbor->GetName().at(0);
                    atoms_2.push_back(c2c1_neighbor);
                }
            }
            pattern << "C";
            atoms_2.push_back(C2C2);
            AtomVector c2c2_neighbors = C2C2->GetNode()->GetNodeNeighbors();
            for(AtomVector::iterator it = c2c2_neighbors.begin(); it != c2c2_neighbors.end(); it++)
            {
                Atom* c2c2_neighbor = (*it);
                if(c2c2_neighbor->GetId().compare(C2->GetId()) != 0)
                {
                    pattern << c2c2_neighbor->GetName().at(0);
                    atoms_2.push_back(c2c2_neighbor);
                }
            }
            pattern << "C";
            atoms_2.push_back(C2C3);
            AtomVector c2c3_neighbors = C2C3->GetNode()->GetNodeNeighbors();
            for(AtomVector::iterator it = c2c3_neighbors.begin(); it != c2c3_neighbors.end(); it++)
            {
                Atom* c2c3_neighbor = (*it);
                if(c2c3_neighbor->GetId().compare(C2->GetId()) != 0)
                {
                    pattern << c2c3_neighbor->GetName().at(0);
                    atoms_2.push_back(c2c3_neighbor);
                }
            }
        }
        if(pattern.str().compare("O-CCCC") == 0 || pattern.str().compare("O-CCHHHCHHHCHHH") == 0 )
        {
            terminal_atoms = atoms_2;
            return "TBT";
        }
    }
    return "";
}
string Assembly::CheckTerminals(Atom* target, AtomVector& terminal_atoms)
{
    if(target != NULL)
    {
        AtomVector o_neighbors = target->GetNode()->GetNodeNeighbors();
        if(CheckROHTerminal(target, terminal_atoms).compare("") != 0)
            return "ROH";
        else if(CheckOMETerminal(target, terminal_atoms).compare("") != 0)
            return "OME";
        else if(CheckTBTTerminal(target, terminal_atoms).compare("") != 0)
            return "TBT";
        else if(o_neighbors.size() == 2)
        {
            Atom* target_o_neighbor = NULL;
            if(o_neighbors.at(0)->GetDescription().find("Het;") != string::npos && o_neighbors.at(1)->GetDescription().find("Het;") == string::npos)
                target_o_neighbor = o_neighbors.at(1);
            else if(o_neighbors.at(0)->GetDescription().find("Het;") == string::npos && o_neighbors.at(1)->GetDescription().find("Het;") != string::npos)
                target_o_neighbor = o_neighbors.at(0);

            if(target_o_neighbor != NULL)
            {
                ResidueVector residues = this->GetAllResiduesOfAssembly();
                Residue* target_residue = NULL;
                for(ResidueVector::iterator it = residues.begin(); it != residues.end(); it++)
                {
                    Residue* residue = *it;
                    if(residue->GetId().compare(target_o_neighbor->GetResidue()->GetId()) == 0)
                    {
                        target_residue = residue;
                        break;
                    }
                }
                if(target_residue != NULL)
                    terminal_atoms = target_residue->GetAtoms();


                AmberGlycamMap amber_glycam = AmberGlycamLookup(target_o_neighbor->GetResidue()->GetName());
                AmberGlycamMap glycam_amber = GlycamAmberLookup(target_o_neighbor->GetResidue()->GetName());

                if(amber_glycam.amber_name_.compare("") != 0)
                    return amber_glycam.amber_name_;
                else if(glycam_amber.glycam_name_.compare("") != 0)
                    return glycam_amber.amber_name_;
                else
                    return target_o_neighbor->GetResidue()->GetName();
            }
            else
                return "";
        }
        else
            return "";
    }
    else
        return "";
}
void Assembly::CheckLinkageNote(Monosaccharide* mono1, Monosaccharide* mono2, string linkage, vector<string>& checked_linkages)
{
    if(find(checked_linkages.begin(), checked_linkages.end(), linkage) == checked_linkages.end())///If this linkage hasn't been checked before by calling the function on other side of the linkage
    {
        vector<string> linkage_tokens = Split(linkage, "-");
        stringstream reverse_linkage;
        reverse_linkage << linkage_tokens.at(2) << "-" << linkage_tokens.at(1) << "-" << linkage_tokens.at(0);
        checked_linkages.push_back(linkage);
        checked_linkages.push_back(reverse_linkage.str());

        int left_c_index = ConvertString<int>(Split(Split(linkage_tokens.at(0), "_").at(0), "C*,\'").at(0));
        int right_c_index = ConvertString<int>(Split(Split(linkage_tokens.at(2), "_").at(0), "C*,\'").at(0));
        int glycosidic_o_index = ConvertString<int>(Split(Split(linkage_tokens.at(1), "_").at(0), "ON*,\'").at(0));
        if(left_c_index != glycosidic_o_index && right_c_index != glycosidic_o_index)
        {
            Note* linkage_note = new Note();
            linkage_note->type_ = Glycan::ERROR;
            linkage_note->category_ = Glycan::GLYCOSIDIC;
            stringstream n;
            n << mono1->sugar_name_.monosaccharide_short_name_ << ": Glycosidic oxygen/nitrogen index does not conform to carbon index in the linkage to "
              << mono2->sugar_name_.monosaccharide_short_name_ << ". " << Split(linkage_tokens.at(0), "_").at(0) << "-" << Split(linkage_tokens.at(1), "_").at(0)
              << "-" << Split(linkage_tokens.at(2), "_").at(0);
            linkage_note->description_ = n.str();
            this->AddNote(linkage_note);
        }
    }
}

void Assembly::BuildOligosaccharideTreeStructure(Monosaccharide *key, vector<Monosaccharide*> values, Oligosaccharide *oligo,
                                                 vector<int>& visited_monos, map<Monosaccharide*, vector<Monosaccharide*> > monos_table,
                                                 map<Monosaccharide*, vector<string> > monos_table_linkages, vector<string>& visited_linkages)
{
    oligo->root_ = key;
    if(values.size() == 0)
    {
        oligo->child_oligos_ = vector<Oligosaccharide*>();
        oligo->child_oligos_linkages_ = vector<string>();
        visited_monos.push_back(key->mono_id);
        return;
    }
    else
    {
        for(vector<Monosaccharide*>::iterator it = values.begin(); it != values.end(); it++)
        {
            Monosaccharide* value_mono = (*it);
            if(find(visited_monos.begin(), visited_monos.end(), value_mono->mono_id) == visited_monos.end())
            {
                int it_index = distance(values.begin(), it);
                vector<string> key_mono_linkages = monos_table_linkages[key];
                string link = key_mono_linkages.at(it_index);
                stringstream reverse_link;
                reverse_link << Split(link, "-").at(2) << "-" << Split(link, "-").at(1) << "-" << Split(link, "-").at(0);
                if(find(visited_linkages.begin(), visited_linkages.end(), link) == visited_linkages.end() &&
                        find(visited_linkages.begin(), visited_linkages.end(), reverse_link.str()) == visited_linkages.end())
                {
                    //                    cout << "key id " << key->mono_id  << ", value id " << value_mono->mono_id << endl;
                    Oligosaccharide* child_oligo = new Oligosaccharide();
                    vector<Monosaccharide*> value_mono_values = monos_table[value_mono];
                    visited_linkages.push_back(link);
                    //                    cout << "call " << value_mono->mono_id << endl;
                    BuildOligosaccharideTreeStructure(value_mono, value_mono_values, child_oligo, visited_monos, monos_table, monos_table_linkages, visited_linkages);
                    oligo->child_oligos_.push_back(child_oligo);
                    oligo->child_oligos_linkages_.push_back(link);
                }
            }
        }
        visited_monos.push_back(key->mono_id);
        return;
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



//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void Assembly::Print(ostream &out)
{
    out << "===================== " << name_ << " ============================" << endl;
    out << "Source file: " << source_file_ << endl;
    if(assemblies_.size() != 0)
    {
        for(AssemblyVector::iterator it = assemblies_.begin(); it != assemblies_.end(); it++)
        {
            Assembly* assembly = (*it);
            assembly->Print(out);
        }
    }
    else
    {
        for(ResidueVector::iterator it = residues_.begin(); it != residues_.end(); it++)
        {
            Residue* residue = (*it);
            residue->Print(out);
        }
    }
}

void Assembly::PrettyPrintHet(ostream &out)
{
    out << "===================== " << "PDB" << " ============================" << endl;
    out << "PDB file name: " << source_file_ << endl;
    if(assemblies_.size() != 0)
    {
        for(AssemblyVector::iterator it = assemblies_.begin(); it != assemblies_.end(); it++)
        {
            Assembly* assembly = (*it);
            assembly->PrettyPrintHet(out);
        }
    }
    else
    {
        for(ResidueVector::iterator it = residues_.begin(); it != residues_.end(); it++)
        {
            Residue* residue = (*it);
            string name = residue->GetName();
            if(name.compare("HOH") != 0)
                residue->PrettyPrintHet(out);
        }
    }
}

void Assembly::PrintHetResidues(ostream &out)
{
    if(assemblies_.size() != 0)
    {
        for(AssemblyVector::iterator it = assemblies_.begin(); it != assemblies_.end(); it++)
        {
            Assembly* assembly = (*it);
            assembly->PrintHetResidues(out);
        }
    }
    else
    {
        for(ResidueVector::iterator it = residues_.begin(); it != residues_.end(); it++)
        {
            Residue* residue = (*it);
            string name = residue->GetName();
            if(name.compare("HOH") != 0)
                residue->PrintHetResidues(out);
        }
    }
}
void Assembly::PrintHetAtoms(ostream &out)
{
    if(assemblies_.size() != 0)
    {
        for(AssemblyVector::iterator it = assemblies_.begin(); it != assemblies_.end(); it++)
        {
            Assembly* assembly = (*it);
            assembly->PrintHetAtoms(out);
        }
    }
    else
    {
        for(ResidueVector::iterator it = residues_.begin(); it != residues_.end(); it++)
        {
            Residue* residue = (*it);
            string name = residue->GetName();
            if(name.compare("HOH") != 0)
                residue->PrintHetAtoms(out);
        }
    }
}

void Assembly::WriteHetResidues(string file_name)
{
    ofstream out_file;
    out_file.open(file_name.c_str());

    for(ResidueVector::iterator it = residues_.begin(); it != residues_.end(); it++)
    {
        Residue* residue = (*it);
        string name = residue->GetName();
        if(name.compare("HOH") != 0)
            residue->WriteHetResidues(out_file);
    }
    out_file.close();
}

void Assembly::WriteHetAtoms(string file_name)
{
    ofstream out_file;
    out_file.open(file_name.c_str());

    for(ResidueVector::iterator it = residues_.begin(); it != residues_.end(); it++)
    {
        Residue* residue = (*it);
        string name = residue->GetName();
        if(name.compare("HOH") != 0)
            residue->WriteHetAtoms(out_file);
    }
}

