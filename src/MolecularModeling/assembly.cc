#include <math.h>
#include <fstream>

#include "../../includes/MolecularModeling/assembly.hpp"
#include "../../includes/MolecularModeling/residue.hpp"
#include "../../includes/MolecularModeling/atom.hpp"
#include "../../includes/MolecularModeling/atomnode.hpp"

#include "../../includes/FileSet/TopologyFileSpace/topologyfile.hpp"
#include "../../includes/FileSet/TopologyFileSpace/topologyassembly.hpp"
#include "../../includes/FileSet/TopologyFileSpace/topologyresidue.hpp"
#include "../../includes/FileSet/TopologyFileSpace/topologyatom.hpp"
#include "../../includes/FileSet/TopologyFileSpace/topologybond.hpp"
#include "../../includes/FileSet/TopologyFileSpace/topologybondtype.hpp"
#include "../../includes/FileSet/TopologyFileSpace/topologyangle.hpp"
#include "../../includes/FileSet/TopologyFileSpace/topologyangletype.hpp"
#include "../../includes/FileSet/TopologyFileSpace/topologydihedral.hpp"
#include "../../includes/FileSet/TopologyFileSpace/topologydihedraltype.hpp"
#include "../../includes/FileSet/TopologyFileSpace/topologyatompair.hpp"
#include "../../includes/FileSet/CoordinateFileSpace/coordinatefile.hpp"
#include "../../includes/ParameterSet/PrepFileSpace/prepfile.hpp"
#include "../../includes/ParameterSet/PrepFileSpace/prepfileresidue.hpp"
#include "../../includes/ParameterSet/PrepFileSpace/prepfileatom.hpp"
#include "../../includes/FileSet/PdbFileSpace/pdbfile.hpp"
#include "../../includes/FileSet/PdbFileSpace/pdbtitlecard.hpp"
#include "../../includes/FileSet/PdbFileSpace/pdbmodelcard.hpp"
#include "../../includes/FileSet/PdbFileSpace/pdbmodel.hpp"
#include "../../includes/FileSet/PdbFileSpace/pdbmodelresidueset.hpp"
#include "../../includes/FileSet/PdbFileSpace/pdbatomcard.hpp"
#include "../../includes/FileSet/PdbFileSpace/pdbheterogenatomcard.hpp"
#include "../../includes/FileSet/PdbFileSpace/pdbatom.hpp"
#include "../../includes/FileSet/PdbFileSpace/pdbconnectcard.hpp"
#include "../../includes/FileSet/PdbFileSpace/pdbfileprocessingexception.hpp"
#include "../../includes/FileSet/PdbqtFileSpace/pdbqtfile.hpp"
#include "../../includes/FileSet/PdbqtFileSpace/pdbqtatom.hpp"
#include "../../includes/FileSet/PdbqtFileSpace/pdbqtmodel.hpp"
#include "../../includes/FileSet/PdbqtFileSpace/pdbqtmodelcard.hpp"
#include "../../includes/FileSet/PdbqtFileSpace/pdbqtatomcard.hpp"
#include "../../includes/FileSet/PdbqtFileSpace/pdbqtmodelresidueset.hpp"
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
#include "../../includes/Geometry/grid.hpp"
#include "../../includes/Geometry/cell.hpp"

//#include "raptor2/raptor.h"
//#include "raptor2/raptor2.h"
//#include "rasqal/rasqal.h"
//#include "redland.h"

using namespace std;
using namespace MolecularModeling;
using namespace TopologyFileSpace;
using namespace CoordinateFileSpace;
using namespace PrepFileSpace;
using namespace PdbFileSpace;
using namespace PdbqtFileSpace;
using namespace ParameterFileSpace;
using namespace Geometry;
using namespace LibraryFileSpace;
using namespace gmml;
using namespace Glycam;

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

//////////////////////////////////////////////////////////
//                       FUNCTIONS                      //
//////////////////////////////////////////////////////////
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
                    new_atom->AddCoordinate(new Geometry::Coordinate(atom->GetAtomOrthogonalCoordinate()));
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
                            PdbAtomCard::PdbAtomMap atom_map = atom_card->GetAtoms();
                            PdbAtom* matching_atom = atom_map[atom->GetAtomSerialNumber()];
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
                                Geometry::Coordinate* coordinate = new Geometry::Coordinate(matching_atom->GetAtomOrthogonalCoordinate());
                                new_atom->AddCoordinate(coordinate);
                                new_atom->SetDescription("Atom;");
                            }
                        }
                        else if(card_index.at(0).compare("HETATOM") == 0)
                        {
                            PdbModelResidueSet::HeterogenAtomCardVector heterogen_atom_cards = residue_set->GetHeterogenAtoms();
                            PdbHeterogenAtomCard* heterogen_atom_card = heterogen_atom_cards.at(gmml::ConvertString<int>(card_index.at(1)));
                            PdbAtomCard::PdbAtomMap heterogen_atom_map = heterogen_atom_card->GetHeterogenAtoms();
                            PdbAtom* matching_heterogen_atom = heterogen_atom_map[atom->GetAtomSerialNumber()];
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
                                Geometry::Coordinate* coordinate = new Geometry::Coordinate(matching_heterogen_atom->GetAtomOrthogonalCoordinate());
                                new_atom->AddCoordinate(coordinate);
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
                    new_atom->AddCoordinate(new Geometry::Coordinate(atom->GetAtomOrthogonalCoordinate()));
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
                            PdbAtomCard::PdbAtomMap atom_map = atom_card->GetAtoms();
                            PdbAtom* matching_atom = atom_map[atom->GetAtomSerialNumber()];
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
                                Geometry::Coordinate* coordinate = new Geometry::Coordinate(matching_atom->GetAtomOrthogonalCoordinate());
                                new_atom->AddCoordinate(coordinate);
                                new_atom->SetDescription("Atom;");
                            }
                        }
                        else if(card_index.at(0).compare("HETATOM") == 0)
                        {
                            PdbModelResidueSet::HeterogenAtomCardVector heterogen_atom_cards = residue_set->GetHeterogenAtoms();
                            PdbHeterogenAtomCard* heterogen_atom_card = heterogen_atom_cards.at(gmml::ConvertString<int>(card_index.at(1)));
                            PdbAtomCard::PdbAtomMap heterogen_atom_map = heterogen_atom_card->GetHeterogenAtoms();
                            PdbAtom* matching_heterogen_atom = heterogen_atom_map[atom->GetAtomSerialNumber()];
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
                                Geometry::Coordinate* coordinate = new Geometry::Coordinate(matching_heterogen_atom->GetAtomOrthogonalCoordinate());
                                new_atom->AddCoordinate(coordinate);
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
                    new_atom->AddCoordinate(new Geometry::Coordinate(atom->GetAtomOrthogonalCoordinate()));
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
                            Geometry::Coordinate* coordinate = new Geometry::Coordinate(matching_atom->GetAtomOrthogonalCoordinate());
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
                    new_atom->AddCoordinate(new Geometry::Coordinate(atom->GetAtomOrthogonalCoordinate()));
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
                            Geometry::Coordinate* coordinate = new Geometry::Coordinate(matching_atom->GetAtomOrthogonalCoordinate());
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
    TopologyAssembly::TopologyResidueMap topology_residues = topology_file->GetAssembly()->GetResidues();
    for(TopologyAssembly::TopologyResidueMap::iterator it = topology_residues.begin(); it != topology_residues.end(); it++)
    {
        Residue* assembly_residue = new Residue();
        assembly_residue->SetAssembly(this);
        TopologyResidue* topology_residue = (*it).second;
        string residue_name = topology_residue->GetResidueName();
        assembly_residue->SetName(residue_name);
        stringstream id;
        id << residue_name << "_" << gmml::BLANK_SPACE << "_" << topology_residue->GetIndex() << "_" << gmml::BLANK_SPACE << "_"
           << gmml::BLANK_SPACE << "_" << id_;
        assembly_residue->SetId(id.str());

        TopologyResidue::TopologyAtomMap topology_atoms = topology_residue->GetAtoms();
        int serial_number = 0;
        for(TopologyResidue::TopologyAtomMap::iterator it1 = topology_atoms.begin(); it1 != topology_atoms.end(); it1++)
        {
            serial_number++;
            Atom* assembly_atom = new Atom();
            string atom_name = (*it1).first;
            assembly_atom->SetName(atom_name);
            stringstream atom_id;
            atom_id << atom_name << "_" << serial_number << "_" << id.str();
            assembly_atom->SetId(atom_id.str());
            TopologyAtom* topology_atom = (*it1).second;
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
    TopologyAssembly::TopologyResidueMap topology_residues = topology_file->GetAssembly()->GetResidues();
    for(TopologyAssembly::TopologyResidueMap::iterator it = topology_residues.begin(); it != topology_residues.end(); it++)
    {
        Residue* assembly_residue = new Residue();
        assembly_residue->SetAssembly(this);
        TopologyResidue* topology_residue = (*it).second;
        string residue_name = topology_residue->GetResidueName();
        assembly_residue->SetName(residue_name);
        stringstream id;
        id << residue_name << "_" << gmml::BLANK_SPACE << "_" << topology_residue->GetIndex()
           << "_" << gmml::BLANK_SPACE << "_" << gmml::BLANK_SPACE << "_" << id_;
        assembly_residue->SetId(id.str());

        TopologyResidue::TopologyAtomMap topology_atoms = topology_residue->GetAtoms();
        int serial_number = 0;
        for(TopologyResidue::TopologyAtomMap::iterator it1 = topology_atoms.begin(); it1 != topology_atoms.end(); it1++)
        {
            serial_number++;
            Atom* assembly_atom = new Atom();
            string atom_name = (*it1).first;
            assembly_atom->SetName(atom_name);
            stringstream atom_id;
            atom_id << atom_name << "_" << serial_number << "_" << id.str();
            assembly_atom->SetId(atom_id.str());
            TopologyAtom* topology_atom = (*it1).second;
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
        int serial_number = 0;
        for(LibraryFileResidue::AtomMap::iterator it1 = library_atoms.begin(); it1 != library_atoms.end(); it1++)
        {
            serial_number++;
            Atom* assembly_atom = new Atom();
            LibraryFileAtom* library_atom = (*it1).second;
            string atom_name = library_atom->GetName();
            assembly_atom->SetName(atom_name);
            stringstream atom_id;
            atom_id << atom_name << "_" << serial_number << "_" << id.str();
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
        int serial_number = 0;
        for(LibraryFileResidue::AtomMap::iterator it1 = library_atoms.begin(); it1 != library_atoms.end(); it1++)
        {
            serial_number++;
            Atom* assembly_atom = new Atom();
            LibraryFileAtom* library_atom = (*it1).second;
            string atom_name = library_atom->GetName();
            assembly_atom->SetName(atom_name);
            stringstream atom_id;
            atom_id << atom_name << "_" << serial_number << "_" << id.str();
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
    TopologyAssembly::TopologyResidueMap topology_residues = topology_file->GetAssembly()->GetResidues();
    for(TopologyAssembly::TopologyResidueMap::iterator it = topology_residues.begin(); it != topology_residues.end(); it++)
    {
        Residue* assembly_residue = new Residue();
        assembly_residue->SetAssembly(this);
        TopologyResidue* topology_residue = (*it).second;
        string residue_name = topology_residue->GetResidueName();
        assembly_residue->SetName(residue_name);
        stringstream id;
        id << residue_name << "_" << gmml::BLANK_SPACE << "_" << topology_residue->GetIndex() << "_" << gmml::BLANK_SPACE << "_"
           << gmml::BLANK_SPACE << "_" << id_;
        assembly_residue->SetId(id.str());

        TopologyResidue::TopologyAtomMap topology_atoms = topology_residue->GetAtoms();
        int serial_number = 0;
        for(TopologyResidue::TopologyAtomMap::iterator it1 = topology_atoms.begin(); it1 != topology_atoms.end(); it1++)
        {
            serial_number++;
            Atom* assembly_atom = new Atom();
            string atom_name = (*it1).first;
            assembly_atom->SetName(atom_name);
            stringstream atom_id;
            atom_id << atom_name << "_" << serial_number << "_" << id.str();
            assembly_atom->SetId(atom_id.str());
            TopologyAtom* topology_atom = (*it1).second;

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
            vector<Geometry::Coordinate*> coord_file_coordinates = coordinate_file->GetCoordinates();
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
    TopologyAssembly::TopologyResidueMap topology_residues = topology_file->GetAssembly()->GetResidues();
    for(TopologyAssembly::TopologyResidueMap::iterator it = topology_residues.begin(); it != topology_residues.end(); it++)
    {
        Residue* assembly_residue = new Residue();
        assembly_residue->SetAssembly(this);
        TopologyResidue* topology_residue = (*it).second;
        string residue_name = topology_residue->GetResidueName();
        assembly_residue->SetName(residue_name);
        stringstream id;
        id << residue_name << "_" << gmml::BLANK_SPACE << "_" << topology_residue->GetIndex() << "_" << gmml::BLANK_SPACE << "_"
           << gmml::BLANK_SPACE << "_" << id_;
        assembly_residue->SetId(id.str());

        TopologyResidue::TopologyAtomMap topology_atoms = topology_residue->GetAtoms();
        int serial_number = 0;
        for(TopologyResidue::TopologyAtomMap::iterator it1 = topology_atoms.begin(); it1 != topology_atoms.end(); it1++)
        {
            serial_number++;
            Atom* assembly_atom = new Atom();
            string atom_name = (*it1).first;
            assembly_atom->SetName(atom_name);
            stringstream atom_id;
            atom_id << atom_name << "_" << serial_number << id.str();
            assembly_atom->SetId(atom_id.str());
            TopologyAtom* topology_atom = (*it1).second;

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

            vector<Geometry::Coordinate*> coord_file_coordinates = coordinate_file->GetCoordinates();
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
        int serial_number = 0;
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
        int serial_number = 0;
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

PdbFile* Assembly::BuildPdbFileStructureFromAssembly()
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

    ExtractPdbModelCardFromAssembly(residue_set, serial_number, sequence_number, model_index_);

    model->SetModelResidueSet(residue_set);
    models[1] = model;
    model_card->SetModels(models);
    pdb_file->SetModels(model_card);

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

    return pdbqt_file;
}

void Assembly::ExtractPdbModelCardFromAssembly(PdbModelResidueSet* residue_set, int &serial_number, int &sequence_number, int model_number)
{
    for(AssemblyVector::iterator it = this->assemblies_.begin(); it != this->assemblies_.end(); it++)
    {
        Assembly* assembly = (*it);
        AssemblyVector assemblies = assembly->GetAssemblies();
        for(AssemblyVector::iterator it1 = assemblies.begin(); it1 != assemblies.end(); it1++)
        {
            ExtractPdbModelCardFromAssembly(residue_set, serial_number, sequence_number, model_number);
        }
        PdbAtomCard* atom_card = new PdbAtomCard();
        PdbHeterogenAtomCard* het_atom_card = new PdbHeterogenAtomCard();
        PdbAtomCard::PdbAtomMap atom_map = PdbAtomCard::PdbAtomMap();
        PdbHeterogenAtomCard::PdbHeterogenAtomMap het_atom_map = PdbHeterogenAtomCard::PdbHeterogenAtomMap();
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
                PdbAtom* pdb_atom = new PdbAtom(serial_number, atom->GetName(), ' ', atom->GetResidue()->GetName(), ' ', sequence_number, ' ',
                                                *((atom->GetCoordinates()).at(model_number)), dNotSet, dNotSet, atom->GetElementSymbol(), "");//ss.str());
                if(find(dscr.begin(), dscr.end(), "Atom") != dscr.end())
                {

                    atom_map[serial_number] = pdb_atom;
                    serial_number++;
                }
                else if(find(dscr.begin(), dscr.end(), "Het") != dscr.end())
                {
                    het_atom_map[serial_number] = pdb_atom;
                    serial_number++;
                }
                else
                {
                    atom_map[serial_number] = pdb_atom;
                    serial_number++;
                }
            }
            sequence_number++;
        }
        atom_card->SetAtoms(atom_map);
        het_atom_card->SetHeterogenAtoms(het_atom_map);
        residue_set->AddAtom(atom_card);
        residue_set->AddHeterogenAtom(het_atom_card);
    }
    PdbAtomCard* atom_card = new PdbAtomCard();
    PdbHeterogenAtomCard* het_atom_card = new PdbHeterogenAtomCard();
    PdbAtomCard::PdbAtomMap atom_map = PdbAtomCard::PdbAtomMap();
    PdbHeterogenAtomCard::PdbHeterogenAtomMap het_atom_map = PdbHeterogenAtomCard::PdbHeterogenAtomMap();
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
            PdbAtom* pdb_atom = new PdbAtom(serial_number, atom->GetName(), ' ', atom->GetResidue()->GetName(), ' ', sequence_number, ' ',
                                            *((atom->GetCoordinates()).at(model_number)), dNotSet, dNotSet, atom->GetElementSymbol(), "");//ss.str());

            if(find(dscr.begin(), dscr.end(), "Atom") != dscr.end())
            {

                atom_map[serial_number] = pdb_atom;
                serial_number++;
            }
            else if(find(dscr.begin(), dscr.end(), "Het") != dscr.end())
            {
                het_atom_map[serial_number] = pdb_atom;
                serial_number++;
            }
            else
            {
                atom_map[serial_number] = pdb_atom;
                serial_number++;
            }
        }
        sequence_number++;
    }
    atom_card->SetAtoms(atom_map);
    het_atom_card->SetHeterogenAtoms(het_atom_map);
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

    topology_file->SetNumberOfAtoms(this->CountNumberOfAtoms());
    topology_file->SetNumberOfTypes(this->CountNumberOfAtomTypes());
    topology_file->SetNumberOfBondsIncludingHydrogen(this->CountNumberOfBondsIncludingHydrogen());
    topology_file->SetNumberOfBondsExcludingHydrogen(this->CountNumberOfBondsExcludingHydrogen());
    topology_file->SetNumberOfAnglesIncludingHydrogen(this->CountNumberOfAnglesIncludingHydrogen());
    topology_file->SetNumberOfAnglesExcludingHydrogen(this->CountNumberOfAnglesExcludingHydrogen());
    topology_file->SetNumberOfDihedralsIncludingHydrogen(this->CountNumberOfDihedralsIncludingHydrogen(parameter_file_path));
    topology_file->SetNumberOfDihedralsExcludingHydrogen(this->CountNumberOfDihedralsExcludingHydrogen(parameter_file_path));
    //    topology_file->SetNumberOfHydrogenParameters();
    //    topology_file->SetNumberOfParameters();
    topology_file->SetNumberOfExcludedAtoms(this->CountNumberOfExcludedAtoms());   // Does not match
    topology_file->SetNumberOfResidues(this->CountNumberOfResidues());
    topology_file->SetTotalNumberOfBonds(this->CountNumberOfBondsExcludingHydrogen());
    topology_file->SetTotalNumberOfAngles(this->CountNumberOfAnglesExcludingHydrogen());
    topology_file->SetTotalNumberOfDihedrals(this->CountNumberOfDihedralsExcludingHydrogen(parameter_file_path));
    topology_file->SetNumberOfBondTypes(this->CountNumberOfBondTypes(parameter_file_path));
    topology_file->SetNumberOfAngleTypes(this->CountNumberOfAngleTypes(parameter_file_path));
    topology_file->SetNumberOfDihedralTypes(this->CountNumberOfDihedralTypes(parameter_file_path));
//        topology_file->SetNumberOfAtomTypesInParameterFile();
    //    topology_file->SetNumberOfDistinctHydrogenBonds();
    //    topology_file->SetPerturbationOption();
    //    topology_file->SetNumberOfBondsPerturbed();
    //    topology_file->SetNumberOfAnglesPerturbed();
    //    topology_file->SetNumberOfDihedralsPerturbed();
    //    topology_file->SetNumberOfBondsGroupPerturbed();
    //    topology_file->SetNumberOfAnglesGroupPerturbed();
    //    topology_file->SetNumberOfDihedralsGroupPerturbed();
    //    topology_file->SetStandardPeriodicBoxOption();
    topology_file->SetNumberOfAtomsInLargestResidue(this->CountMaxNumberOfAtomsInLargestResidue());
    //    topology_file->SetCapOption();
    //    topology_file->SetNumberOfExtraPoints();
    //    topology_file->SetNumberOfBeads();
    TopologyAssembly* topology_assembly = new TopologyAssembly();
    ResidueVector assembly_residues = this->GetAllResiduesOfAssembly();
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
        topology_residue->SetResidueName(assembly_residue->GetName());
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
                        topology_atom->AddExcludedAtom(key2.str());
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
                                topology_atom->AddExcludedAtom(key3.str());
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
                                        topology_atom->AddExcludedAtom(key4.str());
                                    }
                                }
                            }
                        }
                    }
                }
            }
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
            stringstream ss;
            ss << atom_pair_type.at(0) << "-" << atom_pair_type.at(1) << " bond type does not exist in the parameter files";
            cout << ss.str() << endl;
            gmml::log(__LINE__, __FILE__, gmml::ERR, ss.str());
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
    atom_pair_name.push_back(assembly_atom->GetName());
    atom_pair_name.push_back(neighbor->GetName());
    reverse_atom_pair_name.push_back(neighbor->GetName());
    reverse_atom_pair_name.push_back(assembly_atom->GetName());
    vector<string> residue_names = vector<string>();;
    vector<string> reverse_residue_names = vector<string>();;
    residue_names.push_back(assembly_atom->GetResidue()->GetName());
    residue_names.push_back(neighbor->GetResidue()->GetName());
    reverse_residue_names.push_back(neighbor->GetResidue()->GetName());
    reverse_residue_names.push_back(assembly_atom->GetResidue()->GetName());
    vector<string> bond = vector<string>();;
    vector<string> reverse_bond = vector<string>();;
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
            stringstream ss;
            ss << atom_pair_type.at(0) << "-" << atom_pair_type.at(1) << " bond type does not exist in the parameter files";
            cout << ss.str() << endl;
            gmml::log(__LINE__, __FILE__, gmml::ERR, ss.str());
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
            stringstream ss;
            ss << angle_type.at(0) << "-" << angle_type.at(1) << "-" << angle_type.at(2) << " angle type does not exist in the parameter files";
            cout << ss.str() << endl;
            gmml::log(__LINE__, __FILE__, gmml::ERR, ss.str());
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
    angle_atom_names.push_back(assembly_atom->GetName());
    angle_atom_names.push_back(neighbor->GetName());
    angle_atom_names.push_back(neighbor_of_neighbor->GetName());
    reverse_angle_atom_names.push_back(neighbor_of_neighbor->GetName());
    reverse_angle_atom_names.push_back(neighbor->GetName());
    reverse_angle_atom_names.push_back(assembly_atom->GetName());

    vector<string> residue_names = vector<string>();
    vector<string> reverse_residue_names = vector<string>();
    residue_names.push_back(assembly_atom->GetResidue()->GetName());
    residue_names.push_back(neighbor->GetResidue()->GetName());
    residue_names.push_back(neighbor_of_neighbor->GetResidue()->GetName());
    reverse_residue_names.push_back(neighbor_of_neighbor->GetResidue()->GetName());
    reverse_residue_names.push_back(neighbor->GetResidue()->GetName());
    reverse_residue_names.push_back(assembly_atom->GetResidue()->GetName());
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
            stringstream ss;
            ss << angle_type.at(0) << "-" << angle_type.at(1) << "-" << angle_type.at(2) << " angle type does not exist in the parameter files";
            cout << ss.str() << endl;
            gmml::log(__LINE__, __FILE__, gmml::ERR, ss.str());
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
        stringstream ss;
        ss << all_atom_type_permutations.at(0).at(0) << "-" << all_atom_type_permutations.at(0).at(1) << "-" << all_atom_type_permutations.at(0).at(2) << "-"
           << all_atom_type_permutations.at(0).at(3) << " dihedral type (or any other permutation of it) does not exist in the parameter files";
//        cout << ss.str() << endl;
        gmml::log(__LINE__, __FILE__, gmml::ERR, ss.str());
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
            stringstream ss;
            ss << all_improper_dihedrals_atom_type_permutations.at(0).at(0) << "-" << all_improper_dihedrals_atom_type_permutations.at(0).at(1) << "-"
               << all_improper_dihedrals_atom_type_permutations.at(0).at(2) << "-" << all_improper_dihedrals_atom_type_permutations.at(0).at(3)
               << " improer dihedral type (or any other permutation of it) does not exist in the parameter files";
//            cout << ss.str() << endl;
            gmml::log(__LINE__, __FILE__, gmml::ERR, ss.str());
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
            dihedral_atom_names.push_back(assembly_atom->GetName());
            dihedral_atom_names.push_back(neighbor->GetName());
            dihedral_atom_names.push_back(neighbor_of_neighbor->GetName());
            dihedral_atom_names.push_back(neighbor_of_neighbor_of_neighbor->GetName());
            reverse_dihedral_atom_names.push_back(neighbor_of_neighbor_of_neighbor->GetName());
            reverse_dihedral_atom_names.push_back(neighbor_of_neighbor->GetName());
            reverse_dihedral_atom_names.push_back(neighbor->GetName());
            reverse_dihedral_atom_names.push_back(assembly_atom->GetName());

            vector<string> residue_names = vector<string>();
            vector<string> reverse_residue_names = vector<string>();
            residue_names.push_back(assembly_atom->GetResidue()->GetName());
            residue_names.push_back(neighbor->GetResidue()->GetName());
            residue_names.push_back(neighbor_of_neighbor->GetResidue()->GetName());
            residue_names.push_back(neighbor_of_neighbor_of_neighbor->GetResidue()->GetName());
            reverse_residue_names.push_back(neighbor_of_neighbor_of_neighbor->GetResidue()->GetName());
            reverse_residue_names.push_back(neighbor_of_neighbor->GetResidue()->GetName());
            reverse_residue_names.push_back(neighbor->GetResidue()->GetName());
            reverse_residue_names.push_back(assembly_atom->GetResidue()->GetName());

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
                dihedral_atom_names1.push_back(neighbor1->GetName());
                dihedral_atom_names1.push_back(neighbor2->GetName());
                dihedral_atom_names1.push_back(assembly_atom->GetName());
                dihedral_atom_names1.push_back(neighbor3->GetName());
                vector<string> dihedral_atom_names2 = vector<string>();
                dihedral_atom_names2.push_back(neighbor1->GetName());
                dihedral_atom_names2.push_back(assembly_atom->GetName());
                dihedral_atom_names2.push_back(neighbor3->GetName());
                dihedral_atom_names2.push_back(neighbor2->GetName());
                vector<string> dihedral_atom_names3 = vector<string>();
                dihedral_atom_names3.push_back(neighbor1->GetName());
                dihedral_atom_names3.push_back(neighbor3->GetName());
                dihedral_atom_names3.push_back(assembly_atom->GetName());
                dihedral_atom_names3.push_back(neighbor2->GetName());

                vector<string> reverse_dihedral_atom_names1 = vector<string>();
                reverse_dihedral_atom_names1.push_back(neighbor3->GetName());
                reverse_dihedral_atom_names1.push_back(assembly_atom->GetName());
                reverse_dihedral_atom_names1.push_back(neighbor2->GetName());
                reverse_dihedral_atom_names1.push_back(neighbor1->GetName());
                vector<string> reverse_dihedral_atom_names2 = vector<string>();
                reverse_dihedral_atom_names2.push_back(neighbor2->GetName());
                reverse_dihedral_atom_names2.push_back(neighbor3->GetName());
                reverse_dihedral_atom_names2.push_back(assembly_atom->GetName());
                reverse_dihedral_atom_names2.push_back(neighbor1->GetName());
                vector<string> reverse_dihedral_atom_names3 = vector<string>();
                reverse_dihedral_atom_names3.push_back(neighbor2->GetName());
                reverse_dihedral_atom_names3.push_back(assembly_atom->GetName());
                reverse_dihedral_atom_names3.push_back(neighbor3->GetName());
                reverse_dihedral_atom_names3.push_back(neighbor1->GetName());

                vector<string> residue_names1 = vector<string>();
                residue_names1.push_back(neighbor1->GetResidue()->GetName());
                residue_names1.push_back(neighbor2->GetResidue()->GetName());
                residue_names1.push_back(assembly_atom->GetResidue()->GetName());
                residue_names1.push_back(neighbor3->GetResidue()->GetName());
                vector<string> residue_names2 = vector<string>();
                residue_names2.push_back(neighbor1->GetName());
                residue_names2.push_back(assembly_atom->GetResidue()->GetName());
                residue_names2.push_back(neighbor3->GetResidue()->GetName());
                residue_names2.push_back(neighbor2->GetResidue()->GetName());
                vector<string> residue_names3 = vector<string>();
                residue_names3.push_back(neighbor1->GetResidue()->GetName());
                residue_names3.push_back(neighbor3->GetResidue()->GetName());
                residue_names3.push_back(assembly_atom->GetResidue()->GetName());
                residue_names3.push_back(neighbor2->GetResidue()->GetName());

                vector<string> reverse_residue_names1 = vector<string>();
                reverse_residue_names1.push_back(neighbor3->GetResidue()->GetName());
                reverse_residue_names1.push_back(assembly_atom->GetResidue()->GetName());
                reverse_residue_names1.push_back(neighbor2->GetResidue()->GetName());
                reverse_residue_names1.push_back(neighbor1->GetResidue()->GetName());
                vector<string> reverse_residue_names2 = vector<string>();
                reverse_residue_names2.push_back(neighbor2->GetResidue()->GetName());
                reverse_residue_names2.push_back(neighbor3->GetResidue()->GetName());
                reverse_residue_names2.push_back(assembly_atom->GetResidue()->GetName());
                reverse_residue_names2.push_back(neighbor1->GetName());
                vector<string> reverse_residue_names3 = vector<string>();
                reverse_residue_names3.push_back(neighbor2->GetResidue()->GetName());
                reverse_residue_names3.push_back(assembly_atom->GetResidue()->GetName());
                reverse_residue_names3.push_back(neighbor3->GetResidue()->GetName());
                reverse_residue_names3.push_back(neighbor1->GetResidue()->GetName());

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
                        if(permutation_index % 6 == 2)
                        {
                            topology_dihedral->SetResidueNames(residue_names2);
                            topology_dihedral->SetDihedrals(dihedral_atom_names2);
                        }
                        if(permutation_index % 6 == 4)
                        {
                            topology_dihedral->SetResidueNames(residue_names3);
                            topology_dihedral->SetDihedrals(dihedral_atom_names3);
                        }
                        if(permutation_index % 6 == 1)
                        {
                            topology_dihedral->SetResidueNames(reverse_residue_names1);
                            topology_dihedral->SetDihedrals(reverse_dihedral_atom_names1);
                        }
                        if(permutation_index % 6 == 3)
                        {
                            topology_dihedral->SetResidueNames(reverse_residue_names2);
                            topology_dihedral->SetDihedrals(reverse_dihedral_atom_names2);
                        }
                        if(permutation_index % 6 == 5)
                        {
                            topology_dihedral->SetResidueNames(reverse_residue_names3);
                            topology_dihedral->SetDihedrals(reverse_dihedral_atom_names3);
                        }

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

void Assembly::BuildStructureByDistance(double cutoff, int model_index)
{
    cout << "Building structure by distance ..." << endl;
    gmml::log(__LINE__, __FILE__, gmml::INF, "Building structure by distance ...");
    model_index_ = model_index;
    AtomVector all_atoms_of_assembly = this->GetAllAtomsOfAssembly();
    int i = 0;
    for(AtomVector::iterator it = all_atoms_of_assembly.begin(); it != all_atoms_of_assembly.end() - 1; it++)
    {
        Atom* atom = (*it);
        AtomNode* atom_node;
        if(atom->GetNode() == NULL)
        {
            atom_node = new AtomNode();
            atom_node->SetAtom(atom);
        }
        else
            atom_node = atom->GetNode();
        atom_node->SetId(i);
        i++;
        for(AtomVector::iterator it1 = it + 1; it1 != all_atoms_of_assembly.end(); it1++)
        {
//            if(it != it1)
            {
                Atom* neighbor_atom = (*it1);
                if((atom->GetCoordinates().at(model_index)->Distance(*(neighbor_atom->GetCoordinates().at(model_index)))) < cutoff)
                {
                    AtomNode* neighbor_node;
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
                }
            }
        }
        atom->SetNode(atom_node);
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
                   << ":" << gmml::Split(atom_1->GetId(), "_").at(0) << "-"
                   << gmml::Split(atom_2->GetId(), "_").at(2) << "(" << gmml::Split(atom_2->GetId(), "_").at(4) << ")"
                   << ":" << gmml::Split(atom_2->GetId(), "_").at(0);
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
                        vector<string> atom_id_tokens = gmml::Split(assembly_atom_id, "_");
                        stringstream sss;
                        sss << atom_id_tokens.at(2) << ":" << atom_id_tokens.at(0);
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

int Assembly::CountNumberOfBondsIncludingHydrogen()
{
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
            if((atom_name.substr(0,1).compare("H") == 0 ||
                (atom_name.substr(1,1).compare("H") == 0 && isdigit(ConvertString<char>(atom_name.substr(0,1))))))
            {
                counter += node_neighbors.size();
            }
            else
            {
                for(AtomVector::iterator it1 = node_neighbors.begin(); it1 != node_neighbors.end(); it1++)
                {
                    Atom* node_neighbor = (*it1);
                    string node_neighbor_name = node_neighbor->GetName();
                    if((node_neighbor_name.substr(0,1).compare("H") == 0 ||
                        (node_neighbor_name.substr(1,1).compare("H") == 0 && isdigit(ConvertString<char>(node_neighbor_name.substr(0,1))))))
                    {
                        counter++;
                    }
                }
            }
        }
    }
    return counter/2;
}

int Assembly::CountNumberOfBondsExcludingHydrogen()
{
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

            if((atom_name.substr(0,1).compare("H") == 0 ||
                (atom_name.substr(1,1).compare("H") == 0 && isdigit(ConvertString<char>(atom_name.substr(0,1))))))
            {}
            else
            {
                for(AtomVector::iterator it1 = node_neighbors.begin(); it1 != node_neighbors.end(); it1++)
                {
                    Atom* node_neighbor = (*it1);
                    string node_neighbor_name = node_neighbor->GetName();
                    if((node_neighbor_name.substr(0,1).compare("H") == 0 ||
                        (node_neighbor_name.substr(1,1).compare("H") == 0 && isdigit(ConvertString<char>(node_neighbor_name.substr(0,1))))))
                    {}
                    else
                    {
                        counter++;
                    }
                }
            }
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

int Assembly::CountNumberOfAnglesIncludingHydrogen()
{
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
                            counter++;
                        }
                    }
                }
            }
        }
    }
    return counter/2;
}

int Assembly::CountNumberOfAnglesExcludingHydrogen()
{
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
                            counter++;
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
    for(AtomVector::iterator it = atoms.begin(); it != atoms.end(); it++)
    {
        Atom* atom = (*it);
        stringstream ss;
        ss << atom->GetId();
        AtomNode* node = atom->GetNode();
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
                    excluded_atom_list.push_back(first_order_interaction.str());
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
                            excluded_atom_list.push_back(second_order_interaction.str());
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
                                    excluded_atom_list.push_back(third_order_interaction.str());
                            }
                        }
                    }
                }
            }
        }
    }
    return excluded_atom_list.size();
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

void Assembly::ExtractSugars(vector<string> amino_lib_files)
{
    ResidueNameMap dataset_residue_names = GetAllResidueNamesFromMultipleLibFilesMap(amino_lib_files);

    CycleMap cycles = DetectCyclesByExhaustiveRingPerception();

    //    CycleMap cycles = DetectCyclesByDFS();


    cout << endl << "All detected cycles" << endl;
    gmml::log(__LINE__, __FILE__,  gmml::INF,"All detected cycles");
    for(CycleMap::iterator it = cycles.begin(); it != cycles.end(); it++)
    {
        string cycle_atoms_str = (*it).first;
        cout << cycle_atoms_str << endl;
        gmml::log(__LINE__, __FILE__,  gmml::INF, cycle_atoms_str);
    }

    RemoveFusedCycles(cycles);
    FilterAllCarbonCycles(cycles);
    cout << endl << "Cycles after discarding rings that are all-carbon" << endl;
    gmml::log(__LINE__, __FILE__,  gmml::INF, "Cycles after discarding rings that are all-carbon");
    CycleMap sorted_cycles = CycleMap();
    vector<string> anomeric_carbons_status = vector<string>();
    for(CycleMap::iterator it = cycles.begin(); it != cycles.end(); it++)///find anomeric and sort cycles
    {
        string cycle_atoms_str = (*it).first;
        AtomVector cycle_atoms = (*it).second;
        cout << cycle_atoms_str << endl;
        gmml::log(__LINE__, __FILE__,  gmml::INF, cycle_atoms_str);
        Atom* anomeric = FindAnomericCarbon(anomeric_carbons_status, cycle_atoms, cycle_atoms_str);
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

    vector<Monosaccharide*> monos = vector<Monosaccharide*>();
    int mono_id = 0;
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

        mono->cycle_atoms_str_ = cycle_atoms_str;
        mono->cycle_atoms_ = cycle;
        vector<string> orientations = GetSideGroupOrientations(mono, cycle_atoms_str);

        cout << "Side group atoms: " << endl;
        gmml::log(__LINE__, __FILE__,  gmml::INF, "Side group atoms: ");
        for(vector<AtomVector>::iterator it1 = mono->side_atoms_.begin(); it1 != mono->side_atoms_.end(); it1++)
        {
            AtomVector sides = (*it1);
            stringstream side_atoms;
            if(it1 == mono->side_atoms_.begin())///side atoms of anomeric carbon
            {
                if(sides.at(0) != NULL && sides.at(1) != NULL)
                {
                    side_atoms << "[1] -> " << sides.at(0)->GetId() << ", " << sides.at(1)->GetId();
                    cout << side_atoms.str() << endl;
                    gmml::log(__LINE__, __FILE__,  gmml::INF, side_atoms.str());
                }
                else if(sides.at(1) != NULL)
                {
                    side_atoms << "[1] -> " << sides.at(1)->GetId();
                    cout << side_atoms.str() << endl;
                    gmml::log(__LINE__, __FILE__,  gmml::INF, side_atoms.str());

                }
                else if(sides.at(0) != NULL)
                {
                    side_atoms << "[1] -> " << sides.at(0)->GetId();
                    cout << side_atoms.str() << endl;
                    gmml::log(__LINE__, __FILE__,  gmml::INF, side_atoms.str());

                }
            }
            else if(it1 == mono->side_atoms_.end() - 1)//side atoms of last carbon of the ring
            {
                side_atoms << "[" << mono->cycle_atoms_.size() - 1 << "] -> ";
                if(sides.at(0) != NULL)
                {
                    side_atoms << sides.at(0)->GetId();
                    cout << side_atoms.str() << endl;
                    gmml::log(__LINE__, __FILE__,  gmml::INF, side_atoms.str());
                }
            }
            else if(sides.at(1) != NULL)
            {
                int cycle_atom_index = distance(mono->side_atoms_.begin(), it1);
                side_atoms << "[" << cycle_atom_index + 1 << "] -> " << sides.at(1)->GetId();
                cout << side_atoms.str() << endl;
                gmml::log(__LINE__, __FILE__,  gmml::INF, side_atoms.str());
            }
        }
        stringstream anomeric_status;
        anomeric_status << mono->anomeric_status_ << mono->cycle_atoms_.at(0)->GetId();
        cout << endl << anomeric_status.str() << endl;
        gmml::log(__LINE__, __FILE__,  gmml::INF, anomeric_status.str());


        ChemicalCode* code = BuildChemicalCode(orientations);

        if(code != NULL)
        {
            mono->chemical_code_ = code;
        }
        cout << endl << "Stereo chemistry chemical code:"  << endl;
        gmml::log(__LINE__, __FILE__,  gmml::INF, "Stereo chemistry chemical code:");

        code->Print(cout);
        cout << endl;
        string code_str = code->toString();
        gmml::log(__LINE__, __FILE__,  gmml::INF, code_str);

        ///check for +2 and +3 and update the side atoms
        AtomVector plus_sides = ExtractAdditionalSideAtoms(mono);
        mono->sugar_name_ = SugarStereoChemistryNameLookup(code_str);
        if(plus_sides.size() <= 1)
        {

            ExtractDerivatives(mono, cycle_atoms_str);
            for(map<string, string>::iterator it1 = mono->derivatives_map_.begin(); it1 != mono->derivatives_map_.end(); it1++)
            {
                string key = (*it1).first;
                string value = (*it1).second;
                stringstream derivatives;
                derivatives << "Carbon at position " << key << " is attached to " << value;
                cout << derivatives.str() << endl;
                gmml::log(__LINE__, __FILE__,  gmml::INF, derivatives.str());
            }

            GenerateCompleteSugarName(mono);
        }
        else
        {
            ExtractDerivatives(mono, mono->cycle_atoms_str_);
            for(map<string, string>::iterator it1 = mono->derivatives_map_.begin(); it1 != mono->derivatives_map_.end(); it1++)
            {
                string key = (*it1).first;
                string value = (*it1).second;
                stringstream derivatives;
                derivatives << "Carbon at position " << key << " is attached to " << value;
                cout << derivatives.str() << endl;
                gmml::log(__LINE__, __FILE__,  gmml::INF, derivatives.str());
            }

            if(plus_sides.size() == 3)
            {
                vector<string>::iterator index_it;
                if((index_it = find(mono->chemical_code_->right_up_.begin(), mono->chemical_code_->right_up_.end(), "+1")) != mono->chemical_code_->right_up_.end())
                {
                    //check R or S
                    stringstream plus_one;
                    string orientation = CalculateRSOrientations(mono->cycle_atoms_.at(mono->cycle_atoms_.size() - 2), plus_sides.at(0), plus_sides.at(1));
                    plus_one << "+1" << orientation;
                    (*index_it) = plus_one.str();
                    stringstream plus_two;
                    orientation = CalculateRSOrientations(plus_sides.at(0), plus_sides.at(1), plus_sides.at(2));
                    plus_two << "+2" << orientation;
                    mono->chemical_code_->right_up_.push_back(plus_two.str());
                    mono->chemical_code_->right_up_.push_back("+3");
                }
                else if((index_it = find(mono->chemical_code_->right_down_.begin(), mono->chemical_code_->right_down_.end(), "+1")) != mono->chemical_code_->right_down_.end())
                {
                    //check R or S
                    stringstream plus_one;
                    string orientation = CalculateRSOrientations(mono->cycle_atoms_.at(mono->cycle_atoms_.size() - 2), plus_sides.at(0), plus_sides.at(1));
                    plus_one << "+1" << orientation;
                    (*index_it) = plus_one.str();
                    stringstream plus_two;
                    orientation = CalculateRSOrientations(plus_sides.at(0), plus_sides.at(1), plus_sides.at(2));
                    plus_two << "+2" << orientation;
                    mono->chemical_code_->right_down_.push_back(plus_two.str());
                    mono->chemical_code_->right_down_.push_back("+3");
                }
                //update chemical code
                UpdateComplexSugarChemicalCode(mono);

                cout << "Complex structure side group atoms: " << endl;
                gmml::log(__LINE__, __FILE__,  gmml::INF, "Complex structure side group atoms: ");
                for(vector<AtomVector>::iterator it1 = mono->side_atoms_.begin(); it1 != mono->side_atoms_.end(); it1++)
                {
                    stringstream complex_structure_side;
                    AtomVector sides = (*it1);
                    if(it1 == mono->side_atoms_.begin())///side atoms of anomeric carbon
                    {
                        if(sides.at(0) != NULL && sides.at(1) != NULL)
                        {
                            complex_structure_side << "[1] -> " << sides.at(0)->GetId() << ", " << sides.at(1)->GetId();
                            cout << complex_structure_side.str() << endl;
                            gmml::log(__LINE__, __FILE__,  gmml::INF, complex_structure_side.str());

                        }
                        else if(sides.at(1) != NULL)
                        {
                            complex_structure_side << "[1] -> " << sides.at(1)->GetId() ;
                            cout << complex_structure_side.str() << endl;
                            gmml::log(__LINE__, __FILE__,  gmml::INF, complex_structure_side.str());
                        }
                        else if(sides.at(0) != NULL)
                        {
                            complex_structure_side << "[1] -> " << sides.at(0)->GetId();
                            cout << complex_structure_side.str() << endl;
                            gmml::log(__LINE__, __FILE__,  gmml::INF, complex_structure_side.str());
                        }
                    }
                    else if(it1 == mono->side_atoms_.end() - 1)//side atoms of last carbon of the ring
                    {
                        complex_structure_side << "[" << mono->cycle_atoms_.size() - 1 << "]";
                        for(int i = 0; i < plus_sides.size() ; i++)
                            complex_structure_side << " -> " << sides.at(i)->GetId();
                        cout << complex_structure_side.str() << endl;
                        gmml::log(__LINE__, __FILE__,  gmml::INF, complex_structure_side.str());
                    }
                    else if(sides.at(1) != NULL)
                    {
                        int cycle_atom_index = distance(mono->side_atoms_.begin(), it1);                        
                        complex_structure_side << "[" << cycle_atom_index + 1 << "] -> " << sides.at(1)->GetId();
                        cout << complex_structure_side.str() << endl;
                        gmml::log(__LINE__, __FILE__,  gmml::INF, complex_structure_side.str());
                    }
                }
                cout << endl << "Complex sugar chemical code:" << endl;
                gmml::log(__LINE__, __FILE__,  gmml::INF, "Complex sugar chemical code:");
                gmml::log(__LINE__, __FILE__,  gmml::INF, mono->chemical_code_->toString());
                mono->chemical_code_->Print(cout);
                //lookup in complex map
                mono->sugar_name_ = ComplexSugarNameLookup(mono->chemical_code_->toString());

            }
            else if(plus_sides.size() == 2)
            {
                vector<string>::iterator index_it;
                if((index_it = find(mono->chemical_code_->right_up_.begin(), mono->chemical_code_->right_up_.end(), "+1")) != mono->chemical_code_->right_up_.end())
                {
                    //check R or S
                    stringstream plus_one;
                    string orientation = CalculateRSOrientations(mono->cycle_atoms_.at(mono->cycle_atoms_.size() - 2), plus_sides.at(0), plus_sides.at(1));
                    plus_one << "+1" << orientation;
                    (*index_it) = plus_one.str();
                    mono->chemical_code_->right_up_.push_back("+2");
                }
                else if((index_it = find(mono->chemical_code_->right_down_.begin(), mono->chemical_code_->right_down_.end(), "+1")) != mono->chemical_code_->right_down_.end())
                {
                    //check R or S
                    stringstream plus_one;
                    string orientation = CalculateRSOrientations(mono->cycle_atoms_.at(mono->cycle_atoms_.size() - 2), plus_sides.at(0), plus_sides.at(1));
                    plus_one << "+1" << orientation;
                    (*index_it) = plus_one.str();
                    mono->chemical_code_->right_down_.push_back("+2");
                }

                //update chemical code
                UpdateComplexSugarChemicalCode(mono);

                cout << "Complex structure side group atoms: " << endl;
                gmml::log(__LINE__, __FILE__,  gmml::INF, "Complex structure side group atoms: ");
                for(vector<AtomVector>::iterator it1 = mono->side_atoms_.begin(); it1 != mono->side_atoms_.end(); it1++)
                {
                    stringstream complex_sugar_side;
                    AtomVector sides = (*it1);
                    cout << "side atom size" << sides.size() << endl;
                    if(it1 == mono->side_atoms_.begin())///side atoms of anomeric carbon
                    {
                        if(sides.at(0) != NULL && sides.at(1) != NULL)
                        {
                            complex_sugar_side << "[1] -> " << sides.at(0)->GetId() << ", " << sides.at(1)->GetId();
                            cout << complex_sugar_side.str() << endl;
                            gmml::log(__LINE__, __FILE__,  gmml::INF, complex_sugar_side.str());
                        }
                        else if(sides.at(1) != NULL)
                        {
                            complex_sugar_side << "[1] -> " << sides.at(1)->GetId();
                            cout << complex_sugar_side.str() << endl;
                            gmml::log(__LINE__, __FILE__,  gmml::INF, complex_sugar_side.str());
                        }
                        else if(sides.at(0) != NULL)
                        {
                            complex_sugar_side << "[1] -> " << sides.at(0)->GetId();
                            cout << complex_sugar_side.str() << endl;
                            gmml::log(__LINE__, __FILE__,  gmml::INF, complex_sugar_side.str());
                        }
                    }
                    else if(it1 == mono->side_atoms_.end() - 1)//side atoms of last carbon of the ring
                    {
                        complex_sugar_side << "[" << mono->cycle_atoms_.size() - 1 << "]";
                        for(int i = 0; i < plus_sides.size() ; i++)
                        {
                            complex_sugar_side << " -> " << sides.at(i)->GetId();
                        }
                        cout << complex_sugar_side.str() << endl;
                        gmml::log(__LINE__, __FILE__,  gmml::INF, complex_sugar_side.str());
                    }
                    else if(sides.at(1) != NULL)
                    {
                        int cycle_atom_index = distance(mono->side_atoms_.begin(), it1);
                        complex_sugar_side << "[" << cycle_atom_index + 1 << "] -> " << sides.at(1)->GetId();
                        cout << complex_sugar_side.str() << endl;
                        gmml::log(__LINE__, __FILE__,  gmml::INF, complex_sugar_side.str());
                    }
                }

                cout << endl << "Complex sugar chemical code:" << endl;
                gmml::log(__LINE__, __FILE__,  gmml::INF, "Complex sugar chemical code:");
                gmml::log(__LINE__, __FILE__,  gmml::INF, mono->chemical_code_->toString());
                mono->chemical_code_->Print(cout);
                //lookup in complex map
                mono->sugar_name_ = ComplexSugarNameLookup(mono->chemical_code_->toString());
                //generate complete name
                GenerateCompleteSugarName(mono);
            }

        }
        cout << endl;
        if(mono->sugar_name_.monosaccharide_stereochemistry_name_.compare("") == 0 && mono->sugar_name_.monosaccharide_name_.compare("") == 0)
        {
            mono->sugar_name_.monosaccharide_stereochemistry_name_ = "Unknown";
            mono->sugar_name_.monosaccharide_stereochemistry_short_name_ = "Unknown";
            mono->sugar_name_.monosaccharide_name_ = "Unknown";
            mono->sugar_name_.monosaccharide_short_name_ = "Unkown";
        }
        stringstream stereo;
        stereo << "Stereochemistry name: " << mono->sugar_name_.monosaccharide_stereochemistry_name_;
        cout << stereo.str() << endl;
        gmml::log(__LINE__, __FILE__,  gmml::INF, stereo.str());

        stringstream stereo_short;
        stereo_short << "Stereochemistry short name: " << mono->sugar_name_.monosaccharide_stereochemistry_short_name_;
        cout << stereo_short.str() << endl;
        gmml::log(__LINE__, __FILE__,  gmml::INF, stereo_short.str());

        stringstream mono_name;
        mono_name << "Complete name: " << mono->sugar_name_.monosaccharide_name_;
        cout << mono_name.str() << endl;
        gmml::log(__LINE__, __FILE__,  gmml::INF, mono_name.str());

        stringstream mono_short;
        mono_short << "Short name: " << mono->sugar_name_.monosaccharide_short_name_;
        cout << mono_short.str() << endl;
        gmml::log(__LINE__, __FILE__,  gmml::INF, mono_short.str());

        cout << "-------------------------------------------------------------------------------------------------------------------------------------------" << endl;

        mono_id++;
        mono->mono_id = mono_id;
        monos.push_back(mono);
    }
    cout << endl << "Oligosaccharides:" << endl;
    gmml::log(__LINE__, __FILE__,  gmml::INF, "Oligosaccharides:");
    string terminal_residue_name = "";
    vector<Oligosaccharide*> oligosaccharides = ExtractOligosaccharides(monos, dataset_residue_names, terminal_residue_name);
    for(vector<Oligosaccharide*>::iterator it = oligosaccharides.begin(); it != oligosaccharides.end(); it++)
        (*it)->Print(terminal_residue_name, cout);
}

//////////populate ontology

Assembly::CycleMap Assembly::DetectCyclesByExhaustiveRingPerception()
{
    CycleMap cycles = CycleMap();
    AtomVector atoms = GetAllAtomsOfAssemblyExceptProteinWaterResiduesAtoms();
    map<string, Atom*> IdAtom = map<string, Atom*>();

    for(AtomVector::iterator it = atoms.begin(); it != atoms.end(); it++)
    {
        Atom* atom = (*it);
        IdAtom[atom->GetId()] = atom;
    }
    ///Pruning the graph (filter out atoms with less than 2 neighbors)
    PruneGraph(atoms);
    vector<string> path_graph_edges = vector<string> ();
    vector<string> path_graph_labels = vector<string> ();

    ///Converting the molecular graph into a path graph
    ConvertIntoPathGraph(path_graph_edges, path_graph_labels, atoms);

    vector<string> reduced_path_graph_edges = path_graph_edges;
    vector<string> reduced_path_graph_labels = path_graph_labels;

    ///Reducing the path graph
    vector<string> cycless = vector<string>();
    int neighbor_counter = 2;
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
        //            cout << "================================" << (*common_atom_it)->GetId() << "=====================================" << endl;
        atoms.erase(common_atom_it);

        path_graph_edges = reduced_path_graph_edges;
        path_graph_labels = reduced_path_graph_labels;
        //            cout << "---------------------------------------------------------------------" << endl;
        //            for(int i = 0; i < reduced_path_graph_edges.size(); i ++)
        //            {
        //                cout << reduced_path_graph_edges.at(i) << "-------------->" << reduced_path_graph_labels.at(i) << endl;
        //            }
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
}

void Assembly::ReducePathGraph(vector<string> path_graph_edges, vector<string> path_graph_labels, vector<string>& reduced_path_graph_edges, vector<string>& reduced_path_graph_labels, string common_atom, vector<string>& cycles)
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
            {///if the edges a != b and b != c

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
                    cycles.push_back(new_label.str());
                    //                    cout << "CYCLE FOUND " << endl;
                    //                    for(int i = 0; i < path_graph_edges.size(); i++)
                    //                    {
                    //                        string edge = path_graph_edges.at(i);
                    //                        if(edge.find(new_edge_atoms.at(1)) != string::npos)
                    //                            to_be_deleted_edges.push_back(i);
                    //                    }
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

Atom* Assembly::FindAnomericCarbon(vector<string>& anomeric_carbons_status, AtomVector cycle, string cycle_atoms_str)
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
//                    cout << "Anomeric carbon is: " << anomeric_carbon->GetId() << endl;
                    anomeric_carbons_status.push_back("Anomeric carbon: ");
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
//                cout << "Anomeric carbon probably is: " << o_neighbor1->GetId() << endl;
                anomeric_carbons_status.push_back("Anomeric carbon probably is: ");
                return o_neighbor1;
            }
            if(ConvertString<int>(ss2.str()) < ConvertString<int>(ss1.str()))
            {
//                cout << "Anomeric carbon probably is: " << o_neighbor2->GetId() << endl;
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
//                cout << "Anomeric carbon probably is: " << o_neighbor2->GetId() << endl;
                anomeric_carbons_status.push_back("Anomeric carbon probably is: ");
                return o_neighbor2;
            }
            else if(!neighbor2_is_anomeric)
            {
//                cout << "Anomeric carbon probably is: " << o_neighbor1->GetName() << endl;
                anomeric_carbons_status.push_back("Anomeric carbon probably is: ");
                return o_neighbor1;
            }

//            cout << "Not enough information to detect the anomeric carbon, it has been chosen randomely: " << o_neighbor1->GetId() << endl;
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
                Atom* a0 = cycle.at(0);
                if(a0->GetName().substr(0,1).compare("O") == 0)///a0 is oxygen so the vector is in reverse order
                {
                    for(AtomVector::iterator it1 = it - 1; it1 != cycle.begin(); it1--)///atoms before the anomeric atom in reverse order
                    {
                        Atom* a = (*it1);
                        sorted_cycle.push_back(a);
                        sorted_cycle_stream << a->GetId() << "-";
                    }
                    sorted_cycle.push_back((*cycle.begin()));
                    sorted_cycle_stream << (*cycle.begin())->GetId();
                }
                else
                {
                    for(AtomVector::iterator it1 = cycle.begin(); it1 != it; it1++)///atoms before the anomeric atom from beginning of vector
                    {
                        Atom* a = (*it1);
                        sorted_cycle.push_back(a);
                        if(it1 == it -1)
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
                    for(AtomVector::iterator it1 = it; it1 != cycle.begin(); it1--) ///atoms befor anomeric atom to down to beginning of the vector
                    {
                        Atom* a_before = (*it1);
                        sorted_cycle.push_back(a_before);
                        sorted_cycle_stream << a_before->GetId() << "-";
                    }
                    sorted_cycle.push_back((*cycle.begin()));
                    sorted_cycle_stream << (*cycle.begin())->GetId() << "-";
                    for(AtomVector::iterator it2 = cycle.end() - 1; it2 != it; it2--)///atoms from end of the vector down to anomeric atom
                    {
                        Atom* atom_after = (*it2);
                        sorted_cycle.push_back(atom_after);
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
                        if(it1 == cycle.end())
                            sorted_cycle_stream << atom_after->GetId();
                        else
                            sorted_cycle_stream << atom_after->GetId() << "-";
                    }
                    for(AtomVector::iterator it2 = cycle.begin(); it2 != it; it2++)///atoms befor the anomeric atom from beginning of vector
                    {
                        Atom* atom_before = (*it2);
                        sorted_cycle.push_back(atom_before);
                        sorted_cycle_stream << atom_before->GetId() << "-";
                    }
                }
            }
        }
    }
    return sorted_cycle;
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

void Assembly::ExtractDerivatives(Monosaccharide * mono, string cycle_atoms_str)
{
    for(AtomVector::iterator it = mono->cycle_atoms_.begin(); it != mono->cycle_atoms_.end() - 1; it++) ///iterate on cycle atoms except the oxygen in the ring
    {
        int index = distance(mono->cycle_atoms_.begin(), it);
        Atom* target = (*it);
        //        cout << "Target ring: " << target->GetName() << endl;
        string key = "";
        string value = "";
        AtomVector t_neighbors = target->GetNode()->GetNodeNeighbors();
        for(AtomVector::iterator it1 = t_neighbors.begin(); it1 != t_neighbors.end(); it1++)
        {
            Atom* t_neighbor = (*it1);
            if(t_neighbor->GetName().at(0) == 'N' && cycle_atoms_str.find(t_neighbor->GetId()) == string::npos)///check formulas with nitrogen
            {
                if((value = CheckxC_N(target, cycle_atoms_str)).compare("") != 0)///xCH-N
                {
                    //                    cout << "expected pattern: " << "xCH-N" << endl;
                    break;
                }
                if((value = CheckxC_NxO_CO_C(target, cycle_atoms_str, 'N')).compare("") != 0)///xC-N-C=OCH3
                {
                    //                    cout << "expected pattern : " << "xC-N-C=OCH3" << endl;
                    break;
                }
                if((value = CheckxC_NxO_CO_CO(target, cycle_atoms_str, 'N')).compare("") != 0)///xC-N-C=OCH2OH
                {
                    //                    cout << "expected pattern : " << "xC-N-C=OCH2OH" << endl;
                    break;
                }
                if((value = CheckxC_NxO_SO3(target, cycle_atoms_str, 'N')).compare("") != 0)///xC-N-SO3
                {
                    //                    cout << "expected pattern : " << "xC-N-SO3" << endl;
                    break;
                }
                if((value = CheckxC_NxO_PO3(target, cycle_atoms_str, 'N')).compare("") != 0)///xC-N-PO3
                {
                    //                    cout << "expected pattern : " << "xC-N-PO3" << endl;
                    break;
                }
                if((value = CheckxC_NxO_C(target, cycle_atoms_str, 'N')).compare("") != 0)///xC-N-CH3
                {
                    //                    cout << "expected pattern : " << "xC-N-CH3" << endl;
                    break;
                }
            }
            if(t_neighbor->GetName().at(0) == 'O' && cycle_atoms_str.find(t_neighbor->GetId()) == string::npos)///check formulas with oxygen
            {
                if((value = CheckxC_NxO_CO_C(target, cycle_atoms_str, 'O')).compare("") != 0)///xC-O-C=OCH3
                {
                    //                    cout << "expected pattern : " << "xC-O-C=OCH3" << endl;
                    break;
                }
                if((value = CheckxC_NxO_CO_CO(target, cycle_atoms_str, 'O')).compare("") != 0)///xC-O-C=OCH2OH
                {
                    //                    cout << "expected pattern : " << "xC-O-C=OCH2OH" << endl;
                    break;
                }
                if((value = CheckxC_NxO_SO3(target, cycle_atoms_str, 'O')).compare("") != 0)///xC-O-SO3
                {
                    //                    cout << "expected pattern : " << "xC-O-SO3" << endl;
                    break;
                }
                if((value = CheckxC_NxO_PO3(target, cycle_atoms_str, 'O')).compare("") != 0)///xC-O-PO3
                {
                    //                    cout << "expected pattern : " << "xC-O-PO3" << endl;
                    break;
                }
                if((value = CheckxC_NxO_C(target, cycle_atoms_str, 'O')).compare("") != 0)///xC-O-CH3
                {
                    //                    cout << "expected pattern : " << "xC-O-CH3" << endl;
                    break;
                }
                if((value = CheckxCOO(target, cycle_atoms_str)).compare("") != 0)///xC-(O,O) and xC-(O,OH)
                {
                    //                    cout << "expected pattern : " << "xC-(O,O) and xC-(O,OH)" << endl;
                    break;
                }
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
    //    cout << "SDIES! " << endl;
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
            //            cout << "Target side: " << target->GetName() << endl;
            AtomVector t_neighbors = target->GetNode()->GetNodeNeighbors();
            for(AtomVector::iterator it1 = t_neighbors.begin(); it1 != t_neighbors.end(); it1++)
            {
                Atom* t_neighbor = (*it1);
                if(t_neighbor->GetName().at(0) == 'N' && cycle_atoms_str.find(t_neighbor->GetId()) == string::npos)///check formulas with nitrogen
                {
                    if((value = CheckxC_N(target, cycle_atoms_str)).compare("") != 0)///xCH-N
                    {
                        //                        cout << "expected pattern : " << "xCH-N" << endl;
                        break;
                    }
                    if((value = CheckxC_NxO_CO_C(target, cycle_atoms_str, 'N')).compare("") != 0)///xC-N-C=OCH3
                    {
                        //                        cout << "expected pattern : " << "xC-N-C=OCH3" << endl;
                        break;
                    }
                    if((value = CheckxC_NxO_CO_CO(target, cycle_atoms_str, 'N')).compare("") != 0)///xC-N-C=OCH2OH
                    {
                        //                        cout << "expected pattern : " << "xC-N-C=OCH2OH" << endl;
                        break;
                    }
                    if((value = CheckxC_NxO_SO3(target, cycle_atoms_str, 'N')).compare("") != 0)///xC-N-SO3
                    {
                        //                        cout << "expected pattern : " << "xC-N-SO3" << endl;
                        break;
                    }
                    if((value = CheckxC_NxO_PO3(target, cycle_atoms_str, 'N')).compare("") != 0)///xC-N-PO3
                    {
                        //                        cout << "expected pattern : " << "xC-N-PO3" << endl;
                        break;
                    }
                    if((value = CheckxC_NxO_C(target, cycle_atoms_str, 'N')).compare("") != 0)///xC-N-CH3
                    {
                        //                        cout << "expected pattern : " << "xC-N-CH3" << endl;
                        break;
                    }
                }
                if(t_neighbor->GetName().at(0) == 'O' && cycle_atoms_str.find(t_neighbor->GetId()) == string::npos)///check formulas with oxygen
                {
                    if((value = CheckxC_NxO_CO_C(target, cycle_atoms_str, 'O')).compare("") != 0)///xC-O-C=OCH3
                    {
                        //                        cout << "expected pattern : " << "xC-O-C=OCH3" << endl;
                        break;
                    }
                    if((value = CheckxC_NxO_CO_CO(target, cycle_atoms_str, 'O')).compare("") != 0)///xC-O-C=OCH2OH
                    {
                        //                        cout << "expected pattern : " << "xC-O-C=OCH2OH" << endl;
                        break;
                    }
                    if((value = CheckxC_NxO_SO3(target, cycle_atoms_str, 'O')).compare("") != 0)///xC-O-SO3
                    {
                        //                        cout << "expected pattern : " << "xC-O-SO3" << endl;
                        //                        cout << "value : " << value << endl;
                        break;
                    }
                    if((value = CheckxC_NxO_PO3(target, cycle_atoms_str, 'O')).compare("") != 0)///xC-O-PO3
                    {
                        //                        cout << "expected pattern : " << "xC-O-PO3" << endl;
                        break;
                    }
                    if((value = CheckxC_NxO_C(target, cycle_atoms_str, 'O')).compare("") != 0)///xC-O-CH3
                    {
                        //                        cout << "expected pattern : " << "xC-O-CH3" << endl;
                        break;
                    }
                    if((value = CheckxCOO(target, cycle_atoms_str)).compare("") != 0)///xC-(O,O) and xC-(O,OH)
                    {
                        //                        cout << "expected pattern : " << "xC-(O,O) and xC-(O,OH)" << endl;
                        break;
                    }
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

string Assembly::CheckxC_N(Atom* target, string cycle_atoms_str)
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
    //    cout << "CheckxC_N:" << pattern.str() << endl;
    if(pattern.str().compare("xCH-NHH") == 0 || pattern.str().compare("xC-N") == 0 || pattern.str().compare("xCHH-NHH") == 0 )
        return "xCH-N";
    else
        return "";
}

string Assembly::CheckxC_NxO_CO_C(Atom *target, string cycle_atoms_str, char NxO)
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
//        cout << "CheckxC_NxO_CO_C:" << pattern.str() << endl;
    if(NxO == 'N')
    {
        if(pattern.str().compare("xCH-NH-CO-CHHH") == 0 || pattern.str().compare("xC-N-CO-C") == 0 || pattern.str().compare("xCHH-NH-CO-CHHH") == 0 )
            return "xC-N-C=OCH3";
        else
            return "";
    }
    else if(NxO == 'O')
    {
        if(pattern.str().compare("xCH-OH-CO-CHHH") == 0 || pattern.str().compare("xC-O-CO-C") == 0 || pattern.str().compare("xCHH-OH-CO-CHHH") == 0 )
            return "xC-O-C=OCH3";
        else
            return "";
    }
    else
        return "";
}

string Assembly::CheckxC_NxO_CO_CO(Atom *target, string cycle_atoms_str, char NxO)
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
    //    cout << "CheckxC_NxO_CO_CO: " << pattern.str() << endl;
    if(NxO == 'N')
    {
        if(pattern.str().compare("xCH-NH-CO-CHH-OH") == 0 || pattern.str().compare("xC-N-CO-C-O") == 0 || pattern.str().compare("xCHH-NH-CO-CHH-OH") == 0 )
            return "xC-N-C=OCH2OH";
        else
            return "";
    }
    else if(NxO == 'O')
    {
        if(pattern.str().compare("xCH-OH-CO-CHH-OH") == 0 || pattern.str().compare("xC-O-CO-C-O") == 0 || pattern.str().compare("xCHH-OH-CO-CHH-OH") == 0 )
            return "xC-O-C=OCH2OH";
        else
            return "";
    }
    else
        return "";
}

string Assembly::CheckxC_NxO_SO3(Atom *target, string cycle_atoms_str, char NxO)
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
    //    cout << "CheckxC_NxO_SO3: " << pattern.str() << endl;
    if(NxO == 'N')
    {
        if(pattern.str().compare("xCH-NH-SOOOH") == 0 || pattern.str().compare("xCHH-NH-SOOOH") == 0 || pattern.str().compare("xC-N-SOOO") == 0 ||
                pattern.str().compare("xCH-NH-SOOO") == 0 || pattern.str().compare("xCHH-NH-SOOO") == 0)
            return "xC-N-SO3";
        else
            return "";
    }
    else if(NxO == 'O')
    {
        if(pattern.str().compare("xCH-O-SOOOH") == 0 || pattern.str().compare("xCHH-O-SOOOH") == 0 || pattern.str().compare("xC-O-SOOO") == 0 ||
                pattern.str().compare("xCH-O-SOOO") == 0 || pattern.str().compare("xCHH-O-SOOO") == 0)
            return "xC-O-SO3";
        else
            return "";
    }
    else
        return "";
}

string Assembly::CheckxC_NxO_PO3(Atom *target, string cycle_atoms_str, char NxO)
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
    //    cout << "CheckxC_NxO_PO3: " << pattern.str() << endl;
    if(NxO == 'N')
    {
        if(pattern.str().compare("xCH-NH-POOOH") == 0 || pattern.str().compare("xCHH-NH-POOOH") == 0  || pattern.str().compare("xC-N-POOO") == 0 ||
                pattern.str().compare("xCH-NH-POOO") == 0 || pattern.str().compare("xCHH-NH-POOO") == 0)
            return "xC-N-PO3";
        else
            return "";
    }
    else if(NxO = 'O')
    {
        if(pattern.str().compare("xCH-O-POOOH") == 0 || pattern.str().compare("xCHH-O-POOOH") == 0  || pattern.str().compare("xC-O-POOO") == 0 ||
                pattern.str().compare("xCH-O-POOO") == 0 || pattern.str().compare("xCHH-O-POOO") == 0)
            return "xC-O-PO3";
        else
            return "";
    }
    else
        return "";
}

string Assembly::CheckxC_NxO_C(Atom *target, string cycle_atoms_str, char NxO)
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
//        cout << "CheckxC_NxO_C: " << pattern.str() << endl;
    if(NxO == 'N')
    {///xCH-N-CHHH??
        if(pattern.str().compare("xCH-N-CHHH") == 0 || pattern.str().compare("xCH-NH-CHHH") == 0 || pattern.str().compare("xCHH-NH-CHHH") == 0 || pattern.str().compare("xC-N-C") == 0 )
            return "xC-N-CH3";
        else
            return "";
    }
    else if(NxO == 'O')
    {///xCH-O-CHHH??
        if(pattern.str().compare("xCH-O-CHHH") == 0 || pattern.str().compare("xCH-OH-CHHH") == 0 || pattern.str().compare("xCHH-OH-CHHH") == 0 || pattern.str().compare("xC-O-C") == 0 )
            return "xC-O-CH3";
        else
            return "";
    }
    else
        return "";
}

string Assembly::CheckxCOO(Atom *target, string cycle_atoms_str)
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
    //    cout << "CheckxCOO: " << pattern.str() << endl;
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
        cout << "Neutralizing ......." << endl;
        LibraryFile* lib = new LibraryFile(lib_file);
        ParameterFile* param = new ParameterFile(parameter_file, gmml::IONICMOD);
        double charge = this->GetTotalCharge();
        if(fabs(charge) < CHARGE_TOLERANCE)
        {
            cout << "The assembly has 0 charge and is neutral." << endl;
            return;
        }
        else
            cout << "Total charge of the assembly is " << charge << endl;
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
                cout << "The ion has 0 charge" << endl;
                return;
            }
            else if(ion_charge > 0 && charge > 0)
            {
                cout << "The assembly and the given have positive charge, neutralizing process is aborted." << endl;
                return;
            }
            else if(ion_charge < 0 && charge < 0)
            {
                cout << "The assembly and the given have negative charge, neutralizing process is aborted." << endl;
                return;
            }
            else
            {
                int number_of_neutralizing_ion = (int)(fabs(charge) + gmml::CHARGE_TOLERANCE) / (int)(fabs(ion_charge) + gmml::CHARGE_TOLERANCE);
                cout << "The assembly will be neutralized by " << number_of_neutralizing_ion << " ion(s)" << endl;

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
            cout << "The ion has not been found in the library file." << endl;
        }
    }
    else if (ion_count > 0)
    {
        cout << "Ionizing ......." << endl;
        LibraryFile* lib = new LibraryFile(lib_file);
        ParameterFile* param = new ParameterFile(parameter_file, gmml::IONICMOD);
        double charge = this->GetTotalCharge();
        cout << "Total charge of the assembly is " << charge << endl;
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
                cout << "The ion has 0 charge" << endl;
                return;
            }
            else
            {
                cout << "The assembly will be charged by " << ion_count << " ion(s)" << endl;

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
            cout << "The ion has not been found in the library file." << endl;
        }
    }
    else
    {
        cout << "Please have a non-negative number as the number of ion(s) want to add" << endl;
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
//    center_of_geometry = new Coordinate();
    AtomVector all_atoms_of_assembly = this->GetAllAtomsOfAssembly();
    for(AtomVector::iterator it = all_atoms_of_assembly.begin(); it != all_atoms_of_assembly.end(); it++)
    {
        Atom* atom = *it;
        center_of_geometry->Translate(atom->GetCoordinates().at(model_index_)->GetX(),
                                    atom->GetCoordinates().at(model_index_)->GetY(),
                                    atom->GetCoordinates().at(model_index_)->GetZ());
    }
    center_of_geometry->operator /(Coordinate(center_of_geometry->GetX() / all_atoms_of_assembly.size(),
                                              center_of_geometry->GetY() / all_atoms_of_assembly.size(),
                                              center_of_geometry->GetZ() / all_atoms_of_assembly.size()));
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
            cout << "There is no information of the atom type/radius/charge of the atoms in the given library/parameter file" << endl;
            atom->MolecularDynamicAtom::SetRadius(DEFAULT_RADIUS);
            cout << "The default value has been set for " << atom->GetId() << endl;
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
}

void Assembly::GenerateCompleteSugarName(Monosaccharide *mono)
{
    stringstream in_bracket;
    stringstream head;
    stringstream tail;
    for(map<string, string>::iterator it1 = mono->derivatives_map_.begin(); it1 != mono->derivatives_map_.end(); it1++)
    {
        string key = (*it1).first;
        string value = (*it1).second;

        if(value.compare("xCH-N") == 0)
        {
            if(key.compare("a") == 0)
            {
                cout << "CH-N is at warning position: anomeric" << endl;
                gmml::log(__LINE__, __FILE__,  gmml::WAR, "CH-N is at warning position: anomeric");
            }
            else if(key.compare("2") == 0 && mono->sugar_name_.ring_type_.compare("P") == 0 &&
                    find(mono->chemical_code_->right_down_.begin(), mono->chemical_code_->right_down_.end(), "-1") == mono->chemical_code_->right_down_.end() &&
                    find(mono->chemical_code_->right_up_.begin(), mono->chemical_code_->right_up_.end(), "-1") == mono->chemical_code_->right_up_.end() &&
                    mono->sugar_name_.monosaccharide_stereochemistry_name_.compare("") != 0)
            {
                tail << "-osamine";
                stringstream short_name;
                if(mono->sugar_name_.monosaccharide_stereochemistry_short_name_.compare("") != 0)
                {
                    ///moving a, b or x to after the N expression: short-name + N + a/b/x
                    int stereo_condensed_name_size = mono->sugar_name_.monosaccharide_stereochemistry_short_name_.size();
                    string stereo_condensed_name = mono->sugar_name_.monosaccharide_stereochemistry_short_name_;
                    string new_name_part1 = stereo_condensed_name.substr(0, (stereo_condensed_name_size - 1));///short_name
                    char new_name_part2 = stereo_condensed_name.at(stereo_condensed_name_size - 1);///a/b/x
                    short_name << new_name_part1 << "N" << new_name_part2;

                    mono->sugar_name_.monosaccharide_short_name_ = short_name.str();
                }
            }
            else if(mono->sugar_name_.ring_type_.compare("F") == 0 && key.compare("4") == 0)
            {
                cout << "CH-N is at error position: 4" << endl;
                gmml::log(__LINE__, __FILE__,  gmml::ERR, "CH-N is at error position: 4");
            }
            else if(mono->sugar_name_.ring_type_.compare("P") == 0 && key.compare("5") == 0)
            {
                cout << "CH-N is at error position: 5" << endl;
                gmml::log(__LINE__, __FILE__,  gmml::ERR, "CH-N is at error position: 5");
            }
            else
            {
                if(key.compare("-1") == 0)
                    in_bracket << "0N, ";
                else if( key.compare("+1") == 0 || key.compare("+2") == 0 || key.compare("+3") == 0)
                    in_bracket << mono->cycle_atoms_.size() - 1 + ConvertString<int>(key) << "N, ";
                else
                    in_bracket << key << "N, ";
            }
        }

        if(value.compare("xC-N-C=OCH3") == 0)
        {
            if(key.compare("a") == 0)
            {
                cout << "C-N-C=OCH3 is at warning position: anomeric" << endl;
                gmml::log(__LINE__, __FILE__,  gmml::WAR, "C-N-C=OCH3 is at warning position: anomeric");
            }
            else if(key.compare("2") == 0 && mono->sugar_name_.ring_type_.compare("P") == 0 &&
                    find(mono->chemical_code_->right_down_.begin(), mono->chemical_code_->right_down_.end(), "-1") == mono->chemical_code_->right_down_.end() &&
                    find(mono->chemical_code_->right_up_.begin(), mono->chemical_code_->right_up_.end(), "-1") == mono->chemical_code_->right_up_.end() &&
                    mono->sugar_name_.monosaccharide_stereochemistry_name_.compare("") != 0)
            {
                head << "N-acetyl-";
                stringstream short_name;
                if(mono->sugar_name_.monosaccharide_stereochemistry_short_name_.compare("") != 0)
                {
                    ///moving a, b or x to after the NAc expression: short-name + NAc + a/b/x
                    int stereo_condensed_name_size = mono->sugar_name_.monosaccharide_stereochemistry_short_name_.size();
                    string stereo_condensed_name = mono->sugar_name_.monosaccharide_stereochemistry_short_name_;
                    string new_name_part1 = stereo_condensed_name.substr(0, (stereo_condensed_name_size - 1));///short_name
                    char new_name_part2 = stereo_condensed_name.at(stereo_condensed_name_size - 1);///a/b/x
                    short_name << new_name_part1 << "NAc" << new_name_part2;

                    mono->sugar_name_.monosaccharide_short_name_ = short_name.str();
                }
            }

            else if(mono->sugar_name_.ring_type_.compare("F") == 0 && key.compare("4") == 0)
            {
                cout << "C-N-C=OCH3 is at error position: 4" << endl;
                gmml::log(__LINE__, __FILE__,  gmml::ERR, "C-N-C=OCH3 is at error position: 4");
            }
            else if(mono->sugar_name_.ring_type_.compare("P") == 0 && key.compare("5") == 0)
            {
                cout << "C-N-C=OCH3 is at error position: 5" << endl;
                gmml::log(__LINE__, __FILE__,  gmml::ERR, "C-N-C=OCH3 is at error position: 5");
            }
            else
            {
                if(key.compare("-1") == 0)
                    in_bracket << "0NAc, ";
                else if( key.compare("+1") == 0 || key.compare("+2") == 0 || key.compare("+3") == 0)
                    in_bracket << mono->cycle_atoms_.size() - 1 + ConvertString<int>(key) << "NAc, ";
                else
                    in_bracket << key << "NAc, ";
            }
        }
        if(value.compare("xC-N-C=OCH2OH") == 0)
        {
            if(key.compare("a") == 0)
            {
                cout << "C-N-C=OCH2OH is at warning position: anomeric" << endl;
                gmml::log(__LINE__, __FILE__,  gmml::WAR, "C-N-C=OCH2OH is at warning position: anomeric");
            }
            else if(key.compare("2") == 0 && mono->sugar_name_.ring_type_.compare("P") == 0 &&
                    find(mono->chemical_code_->right_down_.begin(), mono->chemical_code_->right_down_.end(), "-1") == mono->chemical_code_->right_down_.end() &&
                    find(mono->chemical_code_->right_up_.begin(), mono->chemical_code_->right_up_.end(), "-1") == mono->chemical_code_->right_up_.end() &&
                    mono->sugar_name_.monosaccharide_stereochemistry_name_.compare("") != 0)
            {
                head << "N-glycolyl-";
                stringstream short_name;
                if(mono->sugar_name_.monosaccharide_stereochemistry_short_name_.compare("") != 0)
                {
                    ///moving a, b or x to after the NGc expression: short-name + NGc + a/b/x
                    int stereo_condensed_name_size = mono->sugar_name_.monosaccharide_stereochemistry_short_name_.size();
                    string stereo_condensed_name = mono->sugar_name_.monosaccharide_stereochemistry_short_name_;
                    string new_name_part1 = stereo_condensed_name.substr(0, (stereo_condensed_name_size - 1));///short_name
                    char new_name_part2 = stereo_condensed_name.at(stereo_condensed_name_size - 1);///a/b/x
                    short_name << new_name_part1 << "NGc" << new_name_part2;

                    mono->sugar_name_.monosaccharide_short_name_ = short_name.str();
                }
            }

            else if(mono->sugar_name_.ring_type_.compare("F") == 0 && key.compare("4") == 0)
            {
                cout << "C-N-C=OCH2OH is at error position: 4" << endl;
                gmml::log(__LINE__, __FILE__,  gmml::ERR, "C-N-C=OCH2OH is at error position: 4");
            }
            else if(mono->sugar_name_.ring_type_.compare("P") == 0 && key.compare("5") == 0)
            {
                cout << "C-N-C=OCH2OH is at error position: 5" << endl;
                gmml::log(__LINE__, __FILE__,  gmml::ERR, "C-N-C=OCH2OH is at error position: 5");
            }
            else
            {
                if(key.compare("-1") == 0)
                    in_bracket << "0NGc, ";
                else if( key.compare("+1") == 0 || key.compare("+2") == 0 || key.compare("+3") == 0)
                    in_bracket << mono->cycle_atoms_.size() - 1 + ConvertString<int>(key) << "NGc, ";
                else
                    in_bracket << key << "NGc, ";
            }
        }
        if(value.compare("xC-N-SO3") == 0)
        {
            if(key.compare("a") == 0)
            {
                cout << "C-N-SO3 is at warning position: anomeric" << endl;
                gmml::log(__LINE__, __FILE__,  gmml::WAR, "C-N-SO3 is at warning position: anomeric");
            }
            else if(key.compare("2") == 0 && mono->sugar_name_.ring_type_.compare("P") == 0 &&
                    find(mono->chemical_code_->right_down_.begin(), mono->chemical_code_->right_down_.end(), "-1") == mono->chemical_code_->right_down_.end() &&
                    find(mono->chemical_code_->right_up_.begin(), mono->chemical_code_->right_up_.end(), "-1") == mono->chemical_code_->right_up_.end() &&
                    mono->sugar_name_.monosaccharide_stereochemistry_name_.compare("") != 0)
            {
                head << "N-sulfo-";
                stringstream short_name;
                if(mono->sugar_name_.monosaccharide_stereochemistry_short_name_.compare("") != 0)
                {
                    ///moving a, b or x to after the NS expression: short-name + NS + a/b/x
                    int stereo_condensed_name_size = mono->sugar_name_.monosaccharide_stereochemistry_short_name_.size();
                    string stereo_condensed_name = mono->sugar_name_.monosaccharide_stereochemistry_short_name_;
                    string new_name_part1 = stereo_condensed_name.substr(0, (stereo_condensed_name_size - 1));///short_name
                    char new_name_part2 = stereo_condensed_name.at(stereo_condensed_name_size - 1);///a/b/x
                    short_name << new_name_part1 << "NS" << new_name_part2;

                    mono->sugar_name_.monosaccharide_short_name_ = short_name.str();
                }
            }
            else if(mono->sugar_name_.ring_type_.compare("F") == 0 && key.compare("4") == 0)
            {
                cout << "C-N-SO3 is at error position: 4" << endl;
                gmml::log(__LINE__, __FILE__,  gmml::ERR, "C-N-SO3 is at error position: 4");
            }
            else if(mono->sugar_name_.ring_type_.compare("P") == 0 && key.compare("5") == 0)
            {
                cout << "C-N-SO3 is at error position: 5" << endl;
                gmml::log(__LINE__, __FILE__,  gmml::ERR, "C-N-SO3 is at error position: 5");
            }
            else
            {
                if(key.compare("-1") == 0)
                    in_bracket << "0NS, ";
                else if( key.compare("+1") == 0 || key.compare("+2") == 0 || key.compare("+3") == 0)
                    in_bracket << mono->cycle_atoms_.size() - 1 + ConvertString<int>(key) << "NS, ";
                else
                    in_bracket << key << "NS, ";
            }
        }
        if(value.compare("xC-N-PO3") == 0)
        {
            if(key.compare("a") == 0)
            {
                cout << "C-N-PO3 is at warning position: anomeric" << endl;
                gmml::log(__LINE__, __FILE__,  gmml::WAR, "C-N-PO3 is at warning position: anomeric");
            }
            else if(key.compare("2") == 0 && mono->sugar_name_.ring_type_.compare("P") == 0 &&
                    find(mono->chemical_code_->right_down_.begin(), mono->chemical_code_->right_down_.end(), "-1") == mono->chemical_code_->right_down_.end() &&
                    find(mono->chemical_code_->right_up_.begin(), mono->chemical_code_->right_up_.end(), "-1") == mono->chemical_code_->right_up_.end() &&
                    mono->sugar_name_.monosaccharide_stereochemistry_name_.compare("") != 0)
            {
                head << "N-phospho-";
                stringstream short_name;
                if(mono->sugar_name_.monosaccharide_stereochemistry_short_name_.compare("") != 0)
                {
                    ///moving a, b or x to after the NP expression: short-name + NP + a/b/x
                    int stereo_condensed_name_size = mono->sugar_name_.monosaccharide_stereochemistry_short_name_.size();
                    string stereo_condensed_name = mono->sugar_name_.monosaccharide_stereochemistry_short_name_;
                    string new_name_part1 = stereo_condensed_name.substr(0, (stereo_condensed_name_size - 1));///short_name
                    char new_name_part2 = stereo_condensed_name.at(stereo_condensed_name_size - 1);///a/b/x
                    short_name << new_name_part1 << "NP" << new_name_part2;

                    mono->sugar_name_.monosaccharide_short_name_ = short_name.str();
                }
            }

            else if(mono->sugar_name_.ring_type_.compare("F") == 0 && key.compare("4") == 0)
            {
                cout << "C-N-PO3 is at error position: 4" << endl;
                gmml::log(__LINE__, __FILE__,  gmml::ERR, "C-N-PO3 is at error position: 4");
            }
            else if(mono->sugar_name_.ring_type_.compare("P") == 0 && key.compare("5") == 0)
            {
                cout << "C-N-PO3 is at error position: 5" << endl;
                gmml::log(__LINE__, __FILE__,  gmml::ERR, "C-N-PO3 is at error position: 5");
            }
            else
            {
                if(key.compare("-1") == 0)
                    in_bracket << "0NP, ";
                else if( key.compare("+1") == 0 || key.compare("+2") == 0 || key.compare("+3") == 0)
                    in_bracket << mono->cycle_atoms_.size() - 1 + ConvertString<int>(key) << "NP, ";
                else
                    in_bracket << key << "NP, ";
            }
        }
        if(value.compare("xC-N-CH3") == 0)
        {
            if(key.compare("a") == 0)
            {
                cout << "C-N-CH3 is at warning position: anomeric" << endl;
                gmml::log(__LINE__, __FILE__,  gmml::WAR, "C-N-CH3 is at warning position: anomeric");
            }
            else if(key.compare("2") == 0 && mono->sugar_name_.ring_type_.compare("P") == 0 &&
                    find(mono->chemical_code_->right_down_.begin(), mono->chemical_code_->right_down_.end(), "-1") == mono->chemical_code_->right_down_.end() &&
                    find(mono->chemical_code_->right_up_.begin(), mono->chemical_code_->right_up_.end(), "-1") == mono->chemical_code_->right_up_.end() &&
                    mono->sugar_name_.monosaccharide_stereochemistry_name_.compare("") != 0)
            {
                head << "N-methyl-";
                stringstream short_name;
                if(mono->sugar_name_.monosaccharide_stereochemistry_short_name_.compare("") != 0)
                {
                    ///moving a, b or x to after the NMe expression: short-name + NMe + a/b/x
                    int stereo_condensed_name_size = mono->sugar_name_.monosaccharide_stereochemistry_short_name_.size();
                    string stereo_condensed_name = mono->sugar_name_.monosaccharide_stereochemistry_short_name_;
                    string new_name_part1 = stereo_condensed_name.substr(0, (stereo_condensed_name_size - 1));///short_name
                    char new_name_part2 = stereo_condensed_name.at(stereo_condensed_name_size - 1);///a/b/x
                    short_name << new_name_part1 << "NMe" << new_name_part2;

                    mono->sugar_name_.monosaccharide_short_name_ = short_name.str();
                }
            }

            else if(mono->sugar_name_.ring_type_.compare("F") == 0 && key.compare("4") == 0)
            {
                cout << "C-N-CH3 is at error position: 4" << endl;
                gmml::log(__LINE__, __FILE__,  gmml::ERR, "C-N-CH3 is at error position: 4");
            }
            else if(mono->sugar_name_.ring_type_.compare("P") == 0 && key.compare("5") == 0)
            {
                cout << "C-N-CH3 is at error position: 5" << endl;
                gmml::log(__LINE__, __FILE__,  gmml::ERR, "C-N-CH3 is at error position: 5");
            }
            else
            {
                if(key.compare("-1") == 0)
                    in_bracket << "0NMe, ";
                else if( key.compare("+1") == 0 || key.compare("+2") == 0 || key.compare("+3") == 0)
                    in_bracket << mono->cycle_atoms_.size() - 1 + ConvertString<int>(key) << "NMe, ";
                else
                    in_bracket << key << "NMe, ";
            }
        }

        if(value.compare("xC-O-C=OCH3") == 0)
        {
            if(mono->sugar_name_.monosaccharide_stereochemistry_name_.compare("") != 0)
            {
                if(key.compare("-1") == 0)
                    head << "0-acetyl-";
                else if(key.compare("+1") == 0 || key.compare("+2") == 0 || key.compare("+3") == 0)
                    head << mono->cycle_atoms_.size() - 1 + ConvertString<int>(key) << "-acetyl-";
                else if(key.compare("a") == 0)
                    head << "1-acetyl-";
                else
                    head << key << "-acetyl-";
            }

            if(mono->sugar_name_.monosaccharide_stereochemistry_short_name_.compare("") != 0)
            {
                if(mono->sugar_name_.ring_type_.compare("F") == 0 && key.compare("4") == 0)
                {
                    cout << "C-O-C=OCH3 is at error position: 4" << endl;
                    gmml::log(__LINE__, __FILE__,  gmml::ERR, "C-O-C=OCH3 is at error position: 4");
                }
                else if(mono->sugar_name_.ring_type_.compare("P") == 0 && key.compare("5") == 0)
                {
                    cout << "C-O-C=OCH3 is at error position: 5" << endl;
                    gmml::log(__LINE__, __FILE__,  gmml::ERR, "C-O-C=OCH3 is at error position: 5");
                }
                else
                {
                    if(key.compare("-1") == 0)
                        in_bracket << "0Ac, ";
                    else if( key.compare("+1") == 0 || key.compare("+2") == 0 || key.compare("+3") == 0)
                        in_bracket << mono->cycle_atoms_.size() - 1 + ConvertString<int>(key) << "Ac, ";
                    else if(key.compare("a") == 0)
                        in_bracket << "1Ac, ";
                    else
                        in_bracket << key << "Ac, ";
                }
            }
        }
        if(value.compare("xC-O-C=OCH2OH") == 0)
        {
            if(mono->sugar_name_.monosaccharide_stereochemistry_name_.compare("") != 0)
            {
                if(key.compare("-1") == 0)
                    head << "0-glycolyl-";
                if(key.compare("+1") == 0 || key.compare("+2") == 0 || key.compare("+3") == 0)
                    head << mono->cycle_atoms_.size() - 1 + ConvertString<int>(key) << "-glycolyl-";
                else
                    head << key << "-glycolyl-";
            }

            if(mono->sugar_name_.monosaccharide_stereochemistry_short_name_.compare("") != 0)
            {
                if(mono->sugar_name_.ring_type_.compare("F") == 0 && key.compare("4") == 0)
                {
                    cout << "C-O-C=OCH2OH is at error position: 4" << endl;
                    gmml::log(__LINE__, __FILE__,  gmml::ERR, "C-O-C=OCH2OH is at error position: 4");
                }
                else if(mono->sugar_name_.ring_type_.compare("P") == 0 && key.compare("5") == 0)
                {
                    cout << "C-O-C=OCH2OH is at error position: 5" << endl;
                    gmml::log(__LINE__, __FILE__,  gmml::ERR, "C-O-C=OCH2OH is at error position: 5");
                }
                else
                {
                    if(key.compare("-1") == 0)
                        in_bracket << "0Gc, ";
                    else if( key.compare("+1") == 0 || key.compare("+2") == 0 || key.compare("+3") == 0)
                        in_bracket << mono->cycle_atoms_.size() - 1 + ConvertString<int>(key) << "Gc, ";
                    else if(key.compare("a") == 0)
                        in_bracket << "1Gc, ";
                    else
                        in_bracket << key << "Gc, ";
                }
            }
        }
        if(value.compare("xC-O-SO3") == 0)
        {
            if(mono->sugar_name_.monosaccharide_stereochemistry_name_.compare("") != 0)
            {
                if(key.compare("-1") == 0)
                    head << "0-sulfo-";
                if(key.compare("+1") == 0 || key.compare("+2") == 0 || key.compare("+3") == 0)
                    head << mono->cycle_atoms_.size() - 1 + ConvertString<int>(key) << "-sulfo-";
                else
                    head << key << "-sulfo-";
            }

            if(mono->sugar_name_.monosaccharide_stereochemistry_short_name_.compare("") != 0)
            {
                if(mono->sugar_name_.ring_type_.compare("F") == 0 && key.compare("4") == 0)
                {
                    cout << "C-O-SO3 is at error position: 4" << endl;
                    gmml::log(__LINE__, __FILE__,  gmml::ERR, "C-O-SO3 is at error position: 4");
                }
                else if(mono->sugar_name_.ring_type_.compare("P") == 0 && key.compare("5") == 0)
                {
                    cout << "C-O-SO3 is at error position: 5" << endl;
                    gmml::log(__LINE__, __FILE__,  gmml::ERR, "C-O-SO3 is at error position: 5");
                }
                else
                {
                    if(key.compare("-1") == 0)
                        in_bracket << "0S, ";
                    else if( key.compare("+1") == 0 || key.compare("+2") == 0 || key.compare("+3") == 0)
                        in_bracket << mono->cycle_atoms_.size() - 1 + ConvertString<int>(key) << "S, ";
                    else if(key.compare("a") == 0)
                        in_bracket << "1S, ";
                    else
                        in_bracket << key << "S, ";
                }
            }
        }
        if(value.compare("xC-O-PO3") == 0)
        {
            if(mono->sugar_name_.monosaccharide_stereochemistry_name_.compare("") != 0)
            {
                if(key.compare("-1") == 0)
                    head << "0-phospho-";
                if(key.compare("+1") == 0 || key.compare("+2") == 0 || key.compare("+3") == 0)
                    head << mono->cycle_atoms_.size() - 1 + ConvertString<int>(key) << "-phospho-";
                else
                    head << key << "-phospho-";
            }

            if(mono->sugar_name_.monosaccharide_stereochemistry_short_name_.compare("") != 0)
            {
                if(mono->sugar_name_.ring_type_.compare("F") == 0 && key.compare("4") == 0)
                {
                    cout << "C-O-PO3 is at error position: 4" << endl;
                    gmml::log(__LINE__, __FILE__,  gmml::ERR, "C-O-PO3 is at error position: 4");
                }
                else if(mono->sugar_name_.ring_type_.compare("P") == 0 && key.compare("5") == 0)
                {
                    cout << "C-O-PO3 is at error position: 5" << endl;
                    gmml::log(__LINE__, __FILE__,  gmml::ERR, "C-O-SO3 is at error position: 5");
                }
                else
                {
                    if(key.compare("-1") == 0)
                        in_bracket << "0P, ";
                    else if( key.compare("+1") == 0 || key.compare("+2") == 0 || key.compare("+3") == 0)
                        in_bracket << mono->cycle_atoms_.size() - 1 + ConvertString<int>(key) << "P, ";
                    else if(key.compare("a") == 0)
                        in_bracket << "1P, ";
                    else
                        in_bracket << key << "P, ";
                }
            }
        }
        if(value.compare("xC-O-CH3") == 0)
        {
            if(mono->sugar_name_.monosaccharide_stereochemistry_name_.compare("") != 0)
            {
                if(key.compare("-1") == 0)
                    head << "0-methyl-";
                if(key.compare("+1") == 0 || key.compare("+2") == 0 || key.compare("+3") == 0)
                    head << mono->cycle_atoms_.size() - 1 + ConvertString<int>(key) << "-methyl-";
                else
                    head << key << "-methyl-";
            }

            if(mono->sugar_name_.monosaccharide_stereochemistry_short_name_.compare("") != 0)
            {
                if(mono->sugar_name_.ring_type_.compare("F") == 0 && key.compare("4") == 0)
                {
                    cout << "C-O-CH3 is at error position: 4" << endl;
                    gmml::log(__LINE__, __FILE__,  gmml::ERR, "C-O-CH3 is at error position: 4");
                }
                else if(mono->sugar_name_.ring_type_.compare("P") == 0 && key.compare("5") == 0)
                {
                    cout << "C-O-CH3 is at error position: 5" << endl;
                    gmml::log(__LINE__, __FILE__,  gmml::ERR, "C-O-CH3 is at error position: 5");
                }
                else
                {
                    if(key.compare("-1") == 0)
                        in_bracket << "0Me, ";
                    else if( key.compare("+1") == 0 || key.compare("+2") == 0 || key.compare("+3") == 0)
                        in_bracket << mono->cycle_atoms_.size() - 1 + ConvertString<int>(key) << "Me, ";
                    else if(key.compare("a") == 0)
                        in_bracket << "1Me, ";
                    else
                        in_bracket << key << "Me, ";
                }
            }
        }
        if(value.compare("xC-(O,OH)") == 0)
        {
            stringstream err_pos;
            if((key.compare("-1") == 0 || key.compare("+2") == 0 || key.compare("+3") == 0) && mono->sugar_name_.monosaccharide_stereochemistry_name_.compare("") != 0)
            {
                err_pos << "C-(O,OH) is at error position: " << key;
                cout << err_pos.str() << endl;
                gmml::log(__LINE__, __FILE__,  gmml::ERR, err_pos.str());

                tail << "-ulosonic acid";
                if(mono->sugar_name_.monosaccharide_stereochemistry_short_name_.compare("") != 0)
                {
                    if(key.compare("-1") == 0)
                        in_bracket << "0AH, ";
                    else if( key.compare("+2") == 0 || key.compare("+3") == 0)
                        in_bracket << mono->cycle_atoms_.size() - 1 + ConvertString<int>(key) << "AH, ";
                }
            }
            else if(key.compare("+1") == 0 && mono->sugar_name_.monosaccharide_stereochemistry_name_.compare("") != 0)
            {
                tail << "-uronic acid";
                stringstream short_name;
                if(mono->sugar_name_.monosaccharide_stereochemistry_short_name_.compare("") != 0)
                {
                    ///moving a, b or x to after the AH expression: short-name + AH + a/b/x
                    int stereo_condensed_name_size = mono->sugar_name_.monosaccharide_stereochemistry_short_name_.size();
                    string stereo_condensed_name = mono->sugar_name_.monosaccharide_stereochemistry_short_name_;
                    string new_name_part1 = stereo_condensed_name.substr(0, (stereo_condensed_name_size - 1));///short_name
                    char new_name_part2 = stereo_condensed_name.at(stereo_condensed_name_size - 1);///a/b/x
                    short_name << new_name_part1 << "AH" << new_name_part2;

                    mono->sugar_name_.monosaccharide_short_name_ = short_name.str();
                }
            }
            else
            {
                err_pos << "C-(O,OH) is at error position: " << key;
                cout << err_pos.str() << endl;
                gmml::log(__LINE__, __FILE__,  gmml::ERR, err_pos.str());
            }
        }
        if(value.compare("xC-(O,O)") == 0)
        {
            stringstream err_pos;
            if((key.compare("-1") == 0 || key.compare("+2") == 0 || key.compare("+3") == 0) && mono->sugar_name_.monosaccharide_stereochemistry_name_.compare("") != 0)
            {
                err_pos << "C-(O,O) is at error position: " << key;
                cout << err_pos.str() << endl;
                gmml::log(__LINE__, __FILE__,  gmml::ERR, err_pos.str());

                tail << "-ulosonate";
                if(mono->sugar_name_.monosaccharide_stereochemistry_short_name_.compare("") != 0)
                {
                    if(key.compare("-1") == 0)
                        in_bracket << "0A, ";
                    else if( key.compare("+2") == 0 || key.compare("+3") == 0)
                        in_bracket << mono->cycle_atoms_.size() - 1 + ConvertString<int>(key) << "A, ";
                }
            }
            else if(key.compare("+1") == 0 && mono->sugar_name_.monosaccharide_stereochemistry_name_.compare("") != 0)
            {
                tail << "-uronate";
                stringstream short_name;
                if(mono->sugar_name_.monosaccharide_stereochemistry_short_name_.compare("") != 0)
                {
                    ///moving a, b or x to after the A expression: short-name + A + a/b/x
                    int stereo_condensed_name_size = mono->sugar_name_.monosaccharide_stereochemistry_short_name_.size();
                    string stereo_condensed_name = mono->sugar_name_.monosaccharide_stereochemistry_short_name_;
                    string new_name_part1 = stereo_condensed_name.substr(0, (stereo_condensed_name_size - 1));///short_name
                    char new_name_part2 = stereo_condensed_name.at(stereo_condensed_name_size - 1);///a/b/x
                    short_name << new_name_part1 << "A" << new_name_part2;

                    mono->sugar_name_.monosaccharide_short_name_ = short_name.str();
                }
            }
            else
            {
                err_pos << "C-(O,O) is at error position: " << key;
                cout << err_pos.str() << endl;
                gmml::log(__LINE__, __FILE__,  gmml::ERR, err_pos.str());
            }
        }
    }
    if(in_bracket.str().size() != 0)
    {
        stringstream short_name;
        if(mono->sugar_name_.monosaccharide_short_name_.compare("") != 0)
            short_name << mono->sugar_name_.monosaccharide_short_name_ << "[" << in_bracket.str().substr(0, in_bracket.str().size() - 2) << "]";
        else
            short_name << mono->sugar_name_.monosaccharide_stereochemistry_short_name_ << "[" << in_bracket.str().substr(0, in_bracket.str().size() - 2) << "]";
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

vector<Oligosaccharide*> Assembly::ExtractOligosaccharides(vector<Monosaccharide*> monos, ResidueNameMap dataset_residue_names, string& terminal_residue_name)
{
    ResidueNameMap common_terminal_residues = gmml::InitializeCommonTerminalResidueMap();
    map<Monosaccharide*, vector<Monosaccharide*> > monos_table = map<Monosaccharide*, vector<Monosaccharide*> >();
    map<Monosaccharide*, vector<string> > monos_table_linkages = map<Monosaccharide*, vector<string> >();
    for(vector<Monosaccharide*>::iterator it = monos.begin(); it != monos.end(); it++)
    {
        Monosaccharide* mono1 = (*it);
        {
            monos_table[mono1] = vector<Monosaccharide*>();
            monos_table_linkages[mono1] = vector<string>();
            for(vector<AtomVector>::iterator it1 = mono1->side_atoms_.begin(); it1 != mono1->side_atoms_.end(); it1++) ///iterate on side atoms
            {
                int index = distance(mono1->side_atoms_.begin(), it1);
                int side_branch_last_carbon_index = 0;
                AtomVector sides = (*it1);
                Atom* target = NULL;
                Atom* target_parent = NULL;
                if(it1 == mono1->side_atoms_.begin())///side atoms of anomeric
                {
                    if(sides.at(1) != NULL)
                    {
                        target = sides.at(1);
                        target_parent = mono1->cycle_atoms_.at(0);
                    }
                }
                else if(it1 == mono1->side_atoms_.end() - 1) ///side atoms of last carbon of the ring
                {
                    if(sides.at(0) != NULL)
                    {
                        for(side_branch_last_carbon_index = sides.size() - 1; sides.at(side_branch_last_carbon_index) == NULL; side_branch_last_carbon_index-- ){}
                        Atom* last_c = sides.at(side_branch_last_carbon_index);
                        AtomVector last_c_neighbors = last_c->GetNode()->GetNodeNeighbors();
                        for(AtomVector::iterator it2 = last_c_neighbors.begin(); it2 != last_c_neighbors.end(); it2++)
                        {
                            if((*it2)->GetId().at(0) == 'O' || (*it2)->GetId().at(0) == 'N')
                            {
                                target = (*it2);
                                target_parent = last_c;
                                break;
                            }

                        }
                    }
                }
                else
                {
                    if(sides.at(1) != NULL)
                    {
                        target = sides.at(1);///index 1 of each side is for non-carbon side atoms in the vector<AtomVector> structure
                        target_parent = mono1->cycle_atoms_.at(index);
                    }
                }
                if(target != NULL)
                {
                    AtomVector t_neighbors = target->GetNode()->GetNodeNeighbors();
                    for(AtomVector::iterator it2 = t_neighbors.begin(); it2 != t_neighbors.end(); it2++)
                    {
                        Atom* t_neighbor = (*it2);
                        if(t_neighbor->GetId().compare(target_parent->GetId()) != 0)///neighbor is not the atom that we are coming from
                        {
                            for(vector<Monosaccharide*>::iterator it3 = monos.begin(); it3 != monos.end(); it3++)
                            {
                                if(it3 != it)
                                {
                                    Monosaccharide* mono2 = (*it3);

                                    int mono2_side_branch_last_carbon_index = 0;
                                    AtomVector mono2_sides = mono2->side_atoms_.at(mono2->side_atoms_.size() - 1);
                                    Atom* mono2_last_c = NULL;
                                    if(mono2_sides.at(0) != NULL)
                                    {
                                        for(mono2_side_branch_last_carbon_index = mono2_sides.size() - 1; mono2_sides.at(mono2_side_branch_last_carbon_index) == NULL; mono2_side_branch_last_carbon_index-- ){}
                                        mono2_last_c = mono2_sides.at(mono2_side_branch_last_carbon_index);
                                    }

                                    if(mono2->cycle_atoms_str_.find(t_neighbor->GetId()) != string::npos ///target atom has been found in another cycle
                                            || (mono2_last_c != NULL && t_neighbor->GetId().compare(mono2_last_c->GetId()) == 0)) ///target atom has been attached to another cycle's side atom
                                    {
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
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    vector<int> visited_monos = vector<int>();
    vector<Oligosaccharide*> oligosaccharides = vector<Oligosaccharide*>();

    for(map<Monosaccharide*, vector<Monosaccharide*> >::iterator it = monos_table.begin(); it != monos_table.end(); it++)
    {
        Monosaccharide* key = (*it).first;
        vector<Monosaccharide*> values = (*it).second;

//        cout << endl << key->mono_id << " >>> " ;
//        for(vector<Monosaccharide*>::iterator it1 = values.begin(); it1 != values.end(); it1++)
//        {
//            Monosaccharide* mono_val = (*it1);
//            cout << mono_val->mono_id << ", ";
//        }
//        cout << endl;
//        vector<string> link = monos_table_linkages[key];
//        for(vector<string>::iterator it1 = link.begin(); it1 != link.end(); it1++)
//        {
//            string l = (*it1);
//            cout << l << ", ";
//        }
//        cout << endl;

        vector<string> visited_linkages = vector<string>();
        if(find(visited_monos.begin(), visited_monos.end(), key->mono_id) == visited_monos.end())///if the mono is not visited
        {
            bool isRoot = false;
            stringstream anomeric_linkage;
            anomeric_linkage << key->cycle_atoms_.at(0)->GetId() << "-";
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
                }
                isRoot = true;
            }
            else if (values.size() == 1)///mono is attached to one other mono
            {
                vector<string> mono_linkage = monos_table_linkages[key];
                vector<string> mono2_linkage = monos_table_linkages[values.at(0)];
                stringstream mono2_anomeric_linkage;
                mono2_anomeric_linkage << values.at(0)->cycle_atoms_.at(0)->GetId() << "-";///atom id on the left side of the linkage c-o-c

                Atom* anomeric_o = NULL;
                if(key->side_atoms_.at(0).at(1) != NULL)
                    anomeric_o = key->side_atoms_.at(0).at(1);

                if(anomeric_o == NULL)
                    isRoot = true;
                else
                {
                    if(dataset_residue_names.find(anomeric_o->GetResidue()->GetName()) != dataset_residue_names.end() ||
                            common_terminal_residues.find(anomeric_o->GetResidue()->GetName()) != common_terminal_residues.end())///mono is attached to a terminal through anomeric oxygen
                    {
                        terminal_residue_name = anomeric_o->GetResidue()->GetName();
                        isRoot = true;
                    }
                    else if((mono_linkage.at(0)).find(anomeric_linkage.str()) == string::npos)///mono is not attached to other mono through anomeric
                        isRoot = true;
                    else if((monos_table_linkages[values.at(0)].size() == 1) && (mono2_linkage.at(0).find(mono2_anomeric_linkage.str()) != string::npos))///mono is attached to other mono through anomeric, the other mono is only attached to this mono through anomeric
                        isRoot = true;
                }
            }
            else if (values.size() > 1)
            {
                Atom* anomeric_o = NULL;
                if(key->side_atoms_.at(0).at(1) != NULL)
                    anomeric_o = key->side_atoms_.at(0).at(1);

                if(anomeric_o == NULL)
                    isRoot = true;
                else
                {
                    if(dataset_residue_names.find(anomeric_o->GetResidue()->GetName()) != dataset_residue_names.end() ||
                            common_terminal_residues.find(anomeric_o->GetResidue()->GetName()) != common_terminal_residues.end())///mono is attached to a terminal through anomeric oxygen
                    {
                        terminal_residue_name = anomeric_o->GetResidue()->GetName();
                        isRoot = true;
                    }
                }
            }
            if(isRoot)
            {
                Oligosaccharide* oligo = new Oligosaccharide();
                BuildOligosaccharideTreeStructure(key, values, oligo, visited_monos, monos_table, monos_table_linkages, visited_linkages);
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
