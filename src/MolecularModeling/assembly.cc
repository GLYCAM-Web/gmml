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

using namespace std;
using namespace MolecularModeling;
using namespace TopologyFileSpace;
using namespace CoordinateFileSpace;
using namespace PrepFileSpace;
using namespace PdbFileSpace;
using namespace ParameterFileSpace;
using namespace Geometry;
using namespace LibraryFileSpace;
using namespace gmml;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
Assembly::Assembly() : description_(""), center_of_geometry_(0,0,0), model_index_(0) {}

Assembly::Assembly(vector<string> file_paths, gmml::InputFileType type)
{
    source_file_type_ = type;
    description_ = "";
    model_index_ = 0;
    switch(type)
    {
        case gmml::PDB:
            source_file_ = file_paths.at(0);
            BuildAssemblyFromPdbFile(source_file_);
            assemblies_ = AssemblyVector();
            break;
        case gmml::TOP:
            source_file_ = file_paths.at(0);
            BuildAssemblyFromTopologyFile(source_file_);
            assemblies_ = AssemblyVector();
            break;
        case gmml::LIB:
            source_file_ = file_paths.at(0);
            BuildAssemblyFromLibraryFile(source_file_);
            assemblies_ = AssemblyVector();
            break;
        case gmml::PREP:
            source_file_ = file_paths.at(0);
            BuildAssemblyFromPrepFile(source_file_);
            assemblies_ = AssemblyVector();
            break;
        case gmml::TOP_CRD:
            source_file_ = file_paths.at(0)+";"+file_paths.at(1);
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
    for(unsigned int i = 0; i < file_paths.size(); i++)
    {
        vector<string> file = file_paths.at(i);
        gmml::InputFileType input_type = types.at(i);
        Assembly* assembly = new Assembly(file, input_type);
        assembly->SetSequenceNumber(i + 1);
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
double Assembly::GetTotalMass()
{
    return total_mass_;
}
Geometry::Coordinate Assembly::GetCenterOfMass()
{
    return center_of_mass_;
}
Geometry::Coordinate Assembly::GetCenterOfGeometry()
{
    return center_of_geometry_;
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
    if(assemblies_.size() == 0)
    {
        assemblies_.push_back(this);
        assembly->SetSequenceNumber(assemblies_.size() + 1);
        assemblies_.push_back(assembly);
        stringstream ss;
        ss << this->name_ << "-" << assembly->GetName();
        this->name_ = ss.str();
        this->residues_ = ResidueVector();
        this->chemical_type_ = "";
        this->sequence_number_ = 1;
        this->total_mass_ = 0;
        this->center_of_geometry_ = Coordinate();
        this->center_of_mass_ = Coordinate();
        stringstream sss;
        sss << this->source_file_ << "#" << assembly->GetSourceFile();
        this->source_file_ = sss.str();
        source_file_type_ = gmml::MULTIPLE;
        model_index_ = 0;
    }
    else
    {
        stringstream ss;
        ss << this->name_ << "-" << assembly->GetName();
        assembly->SetSequenceNumber(assemblies_.size() + 1);
        assemblies_.push_back(assembly);
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
void Assembly::SetTtoalMass(double total_mass)
{
    total_mass_ = total_mass;
}
void Assembly::SetCenterOfMass(Geometry::Coordinate center_of_mass)
{
    center_of_mass_ = center_of_mass;
}
void Assembly::SetCenterOfGeometry(Geometry::Coordinate center_of_geometry)
{
    center_of_geometry_ = center_of_geometry;
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
void Assembly::BuildAssemblyFromPdbFile(string pdb_file_path)
{
    try
    {
        this->ClearAssembly();
        PdbFile* pdb_file = new PdbFile(pdb_file_path);
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
                ss << residue_name << "_" << chain_id << "_" << sequence_number << "_" << insertion_code << "_" << alternate_location;
                string key = ss.str();
                residue->SetId(key);

                Atom* new_atom = new Atom();
                residue->SetName(residue_name);
                string atom_name = atom->GetAtomName();
                new_atom->SetName(atom_name);
                new_atom->MolecularDynamicAtom::SetCharge(gmml::ConvertString<double>(atom->GetAtomCharge()));
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
                                << matching_insertion_code << "_" << matching_alternate_location;
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
                                 << matching_heterogen_insertion_code << "_" << matching_heterogen_alternate_location;
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

void Assembly::BuildAssemblyFromTopologyFile(string topology_file_path)
{
    this->ClearAssembly();
    TopologyFile* topology_file = new TopologyFile(topology_file_path);
    name_ = topology_file->GetTitle();
    sequence_number_ = 1;
    TopologyAssembly::TopologyResidueMap topology_residues = topology_file->GetAssembly()->GetResidues();
    for(TopologyAssembly::TopologyResidueMap::iterator it = topology_residues.begin(); it != topology_residues.end(); it++)
    {
        Residue* assembly_residue = new Residue();
        assembly_residue->SetAssembly(this);
        string residue_name = (*it).first;
        assembly_residue->SetName(residue_name);
        TopologyResidue* topology_residue = (*it).second;
        stringstream id;
        id << residue_name << "(" << topology_residue->GetIndex() << ")" ;
        assembly_residue->SetId(id.str());

        TopologyResidue::TopologyAtomMap topology_atoms = topology_residue->GetAtoms();
        for(TopologyResidue::TopologyAtomMap::iterator it1 = topology_atoms.begin(); it1 != topology_atoms.end(); it1++)
        {
            Atom* assembly_atom = new Atom();
            string atom_name = (*it1).first;
            assembly_atom->SetName(atom_name);
            stringstream atom_id;
            atom_id << id.str() << ":" << atom_name;
            assembly_atom->SetId(atom_id.str());
            TopologyAtom* topology_atom = (*it1).second;
            assembly_atom->MolecularDynamicAtom::SetCharge(topology_atom->GetAtomCharge());
            assembly_atom->MolecularDynamicAtom::SetMass(topology_atom->GetAtomMass());
            assembly_atom->MolecularDynamicAtom::SetRadius(dNotSet); ///////////////////
            assembly_atom->MolecularDynamicAtom::SetAtomType(topology_atom->GetType());

            assembly_atom->SetResidue(assembly_residue);
            assembly_atom->SetName(topology_atom->GetAtomName());

            assembly_residue->AddAtom(assembly_atom);
        }
        this->AddResidue(assembly_residue);

    }
}

void Assembly::BuildAssemblyFromLibraryFile(string library_file_path)
{
    this->ClearAssembly();
    LibraryFile* library_file = new LibraryFile(library_file_path);
    sequence_number_ = 1;
    LibraryFile::ResidueMap library_residues = library_file->GetResidues();
    stringstream ss;

    for(LibraryFile::ResidueMap::iterator it = library_residues.begin(); it != library_residues.end(); it++)
    {
        Residue* assembly_residue = new Residue();
        assembly_residue->SetAssembly(this);
        string residue_name = (*it).first;
        assembly_residue->SetName(residue_name);
        assembly_residue->SetId(residue_name);
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
            Atom* assembly_atom = new Atom();
            LibraryFileAtom* library_atom = (*it1).second;
            string atom_name = library_atom->GetName();
            assembly_atom->SetName(atom_name);
            stringstream atom_id;
            atom_id << residue_name << ":" << atom_name;
            assembly_atom->SetId(atom_id.str());

            assembly_atom->SetResidue(assembly_residue);
            assembly_atom->SetName(library_atom->GetName());

            assembly_atom->MolecularDynamicAtom::SetCharge(library_atom->GetCharge());
            assembly_atom->MolecularDynamicAtom::SetAtomType(library_atom->GetType());

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

void Assembly::BuildAssemblyFromTopologyCoordinateFile(string topology_file_path, string coordinate_file_path)
{
    this->ClearAssembly();
    TopologyFile* topology_file = new TopologyFile(topology_file_path);
    name_ = topology_file->GetTitle();
    sequence_number_ = 1;
    TopologyAssembly::TopologyResidueMap topology_residues = topology_file->GetAssembly()->GetResidues();
    for(TopologyAssembly::TopologyResidueMap::iterator it = topology_residues.begin(); it != topology_residues.end(); it++)
    {
        Residue* assembly_residue = new Residue();
        assembly_residue->SetAssembly(this);
        string residue_name = (*it).first;
        assembly_residue->SetName(residue_name);
        TopologyResidue* topology_residue = (*it).second;
        stringstream id;
        id << residue_name << "(" << topology_residue->GetIndex() << ")" ;
        assembly_residue->SetId(id.str());

        TopologyResidue::TopologyAtomMap topology_atoms = topology_residue->GetAtoms();
        for(TopologyResidue::TopologyAtomMap::iterator it1 = topology_atoms.begin(); it1 != topology_atoms.end(); it1++)
        {
            Atom* assembly_atom = new Atom();
            string atom_name = (*it1).first;
            assembly_atom->SetName(atom_name);
            stringstream atom_id;
            atom_id << id.str() << ":" << atom_name;
            assembly_atom->SetId(atom_id.str());
            TopologyAtom* topology_atom = (*it1).second;

            assembly_atom->MolecularDynamicAtom::SetCharge(topology_atom->GetAtomCharge());
            assembly_atom->MolecularDynamicAtom::SetAtomType(topology_atom->GetType());
            assembly_atom->MolecularDynamicAtom::SetRadius(dNotSet); ////////////////////////
            assembly_atom->MolecularDynamicAtom::SetMass(topology_atom->GetAtomMass());

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

void Assembly::BuildAssemblyFromPrepFile(string prep_file_path)
{
    this->ClearAssembly();
    PrepFile* prep_file = new PrepFile(prep_file_path);
    sequence_number_ = 1;
    PrepFile::ResidueMap prep_residues = prep_file->GetResidues();
    stringstream ss;

    for(PrepFile::ResidueMap::iterator it = prep_residues.begin(); it != prep_residues.end(); it++)
    {
        CoordinateVector cartesian_coordinate_list = CoordinateVector();
        int head_atom_index = INFINITY;
        int tail_atom_index = -INFINITY;
        Atom* head_atom = new Atom();
        Atom* tail_atom = new Atom();

        Residue* assembly_residue = new Residue();
        assembly_residue->SetAssembly(this);
        string residue_name = (*it).first;
        assembly_residue->SetName(residue_name);
        assembly_residue->SetId(residue_name);
        PrepFileResidue* prep_residue = (*it).second;
        string prep_residue_name = prep_residue->GetName();
        if(distance(prep_residues.begin(), it) == (int)prep_residues.size()-1)
            ss << prep_residue_name;
        else
            ss << prep_residue_name << "-";
        PrepFileResidue::PrepFileAtomVector prep_atoms = prep_residue->GetAtoms();
        for(PrepFileResidue::PrepFileAtomVector::iterator it1 = prep_atoms.begin(); it1 != prep_atoms.end(); it1++)
        {
            Atom* assembly_atom = new Atom();
            PrepFileAtom* prep_atom = (*it1);

            assembly_atom->SetResidue(assembly_residue);
            string atom_name = prep_atom->GetName();
            assembly_atom->SetName(atom_name);
            stringstream atom_id;
            atom_id << residue_name << ":" << atom_name;
            assembly_atom->SetId(atom_id.str());

            assembly_atom->MolecularDynamicAtom::SetAtomType(prep_atom->GetType());
            assembly_atom->MolecularDynamicAtom::SetCharge(prep_atom->GetCharge());

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
                    int great_grabdparent_index = prep_atom->GetDihedralIndex() - 1;
                    int grandparent_index = prep_atom->GetAngleIndex() - 1;
                    int parent_index = prep_atom->GetBondIndex() - 1;
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

void Assembly::ExtractPdbModelCardFromAssembly(PdbModelResidueSet* residue_set, int &serial_number, int &sequence_number, int model_number)
{
    cout << "Creating PDB file" << endl;
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
                stringstream ss;
                ss << fixed << setprecision(1) << atom->GetCharge();
                PdbAtom* pdb_atom = new PdbAtom(serial_number, atom->GetName(), ' ', atom->GetResidue()->GetName(), ' ', sequence_number, ' ',
                                                *((atom->GetCoordinates()).at(model_number)), dNotSet, dNotSet, atom->GetElementSymbol(), ss.str());
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
            stringstream ss;
            ss << fixed << setprecision(1) << atom->GetCharge();
            PdbAtom* pdb_atom = new PdbAtom(serial_number, atom->GetName(), ' ', atom->GetResidue()->GetName(), ' ', sequence_number, ' ',
                                            *((atom->GetCoordinates()).at(model_number)), dNotSet, dNotSet, atom->GetElementSymbol(), ss.str());

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
        }
        sequence_number++;
    }
    atom_card->SetAtoms(atom_map);
    het_atom_card->SetHeterogenAtoms(het_atom_map);
    residue_set->AddAtom(atom_card);
    residue_set->AddHeterogenAtom(het_atom_card);
}

PrepFile* Assembly::BuildPrepFileStructureFromAssembly()
{
    PrepFile* prep_file = new PrepFile();
    ResidueVector assembly_residues = this->GetAllResiduesOfAssembly();
    PrepFile::ResidueMap prep_residues = PrepFile::ResidueMap();
    for(ResidueVector::iterator it = assembly_residues.begin(); it != assembly_residues.end(); it++)
    {
        Residue* assembly_residue = *it;
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
            dummy_atom->SetTopologicalType(PrepFileSpace::kTopTypeM);
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
            prep_atom->SetCharge(assembly_atom->GetCharge());
            prep_atom->SetName(assembly_atom->GetName());
            prep_atom->SetTopologicalType(residue_topological_types.at(atom_index - DEFAULT_DUMMY_ATOMS - 1));
            prep_atom->SetType(assembly_atom->GetAtomType());

            int parent_index = prep_atom->GetBondIndex() - 1;
            int grandparent_index = prep_atom->GetAngleIndex() - 1;
            int great_grandparent_index = prep_atom->GetDihedralIndex() - 1;

            coordinate_list.push_back(cartesian_coordinate_list.at(great_grandparent_index));
            coordinate_list.push_back(cartesian_coordinate_list.at(grandparent_index));
            coordinate_list.push_back(cartesian_coordinate_list.at(parent_index));
            cartesian_coordinate_list.push_back(assembly_atom->GetCoordinates().at(0));

            Coordinate* internal_coordinate = gmml::ConvertCartesianCoordinate2InternalCoordinate(assembly_atom->GetCoordinates().at(0),
                                                                                                  coordinate_list);
            prep_atom->SetBondLength(internal_coordinate->GetX());
            prep_atom->SetAngle(internal_coordinate->GetY());
            prep_atom->SetDihedral(internal_coordinate->GetZ());

            atom_index++;
            prep_atoms.push_back(prep_atom);
        }
//        prep_residue->SetImproperDihedrals();
        prep_residue->SetLoops(loops);
        prep_residue->SetAtoms(prep_atoms);
        prep_residue->SetCharge(prep_residue->CalculatePrepResidueCharge());
        prep_residues[assembly_residue->GetName()] = prep_residue;
    }
    //    prep_file->SetPath();
    prep_file->SetResidues(prep_residues);
    return prep_file;
}

vector<TopologicalType> Assembly::GetAllTopologicalTypesOfAtomsOfResidue(AtomVector assembly_atoms, PrepFileResidue::Loop& loops, vector<int> & bond_index)
{
    vector<TopologicalTypeStackElement> stack = vector<TopologicalTypeStackElement>();
    vector<TopologicalType> topological_types = vector<TopologicalType>();
    TopologicalTypeStackElement top_stack = EMPTY;
    vector<string> visited_atom_name = vector<string>();
    vector<int> atom_index_stack = vector<int>();
    for(AtomVector::iterator it1 = assembly_atoms.begin(); it1 != assembly_atoms.end(); it1++)
    {
        Atom* assembly_atom = (*it1);
        int index = distance(assembly_atoms.begin(), it1) + 1;
        if(stack.empty())
            top_stack = EMPTY;
        else
            top_stack = stack.at(stack.size() - 1);


        int visited_neighbors = 0;
        AtomVector neighbors = assembly_atom->GetNode()->GetNodeNeighbors();
        for(AtomVector::iterator it = neighbors.begin(); it != neighbors.end(); it++)
        {
            Atom* neighbor = *it;
            if(find(visited_atom_name.begin(), visited_atom_name.end(), neighbor->GetName()) != visited_atom_name.end())
                visited_neighbors++;
        }
//        cout << visited_neighbors << " " << assembly_atom->GetName() <<  " " << assembly_atom->GetNode()->GetNodeNeighbors().size() << endl;
        switch ( top_stack )
        {
            case EMPTY:
                if(assembly_atom->GetNode()->GetNodeNeighbors().size() == 1)
                {
                    switch(visited_neighbors)
                    {
                        case 0:
                            stack.push_back(M1);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopTypeM);
                            bond_index.push_back(index + DEFAULT_DUMMY_ATOMS - 1);
                            break;
                        case 1:
                            topological_types.push_back(kTopTypeE);
                            bond_index.push_back(index + DEFAULT_DUMMY_ATOMS - 1);
                            break;
                    }
                }
                else if(assembly_atom->GetNode()->GetNodeNeighbors().size() == 2)
                {
                    switch(visited_neighbors)
                    {
                        case 0:
                            stack.push_back(M2);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopTypeM);
                            bond_index.push_back(index + DEFAULT_DUMMY_ATOMS - 1);
                            break;
                        case 1:
                            stack.push_back(M1);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopTypeM);
                            bond_index.push_back(index + DEFAULT_DUMMY_ATOMS - 1);
                            break;
                        case 2:
                            topological_types.push_back(kTopTypeE);
                            bond_index.push_back(index + DEFAULT_DUMMY_ATOMS - 1);
                            break;
                    }
                }
                else if(assembly_atom->GetNode()->GetNodeNeighbors().size() == 3)
                {
                    switch(visited_neighbors)
                    {
                        case 0:
                            stack.push_back(M3);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopTypeM);
                            bond_index.push_back(index + DEFAULT_DUMMY_ATOMS - 1);
                            break;
                        case 1:
                            stack.push_back(M2);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopTypeM);
                            bond_index.push_back(index + DEFAULT_DUMMY_ATOMS - 1);
                            break;
                        case 2:
                            stack.push_back(M1);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopTypeM);
                            bond_index.push_back(index + DEFAULT_DUMMY_ATOMS - 1);
                            break;
                        case 3:
                            topological_types.push_back(kTopTypeE);
                            bond_index.push_back(index + DEFAULT_DUMMY_ATOMS - 1);
                            break;
                    }
                }
                else if(assembly_atom->GetNode()->GetNodeNeighbors().size() == 4)
                {
                    switch(visited_neighbors)
                    {
                        case 0:
                            stack.push_back(M4);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopTypeM);
                            bond_index.push_back(index + DEFAULT_DUMMY_ATOMS - 1);
                            break;
                        case 1:
                            stack.push_back(M3);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopTypeM);
                            bond_index.push_back(index + DEFAULT_DUMMY_ATOMS - 1);
                            break;
                        case 2:
                            stack.push_back(M2);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopTypeM);
                            bond_index.push_back(index + DEFAULT_DUMMY_ATOMS - 1);
                            break;
                        case 3:
                            stack.push_back(M1);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopTypeM);
                            bond_index.push_back(index + DEFAULT_DUMMY_ATOMS - 1);
                            break;
                        case 4:
                            topological_types.push_back(kTopTypeE);
                            bond_index.push_back(index + DEFAULT_DUMMY_ATOMS - 1);
                            break;
                    }
                }
                break;
            case M1:
                if(assembly_atom->GetNode()->GetNodeNeighbors().size() == 1)
                {
                    top_stack = stack.at(stack.size() - 1);
                    switch(visited_neighbors)
                    {
                        case 0:
                            stack.push_back(M1);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopTypeM);
                            bond_index.push_back(index + DEFAULT_DUMMY_ATOMS - 1);
                            break;
                        case 1:
                            stack.pop_back();
                            int parent_index = atom_index_stack.at(atom_index_stack.size() - 1);
                            atom_index_stack.pop_back();
                            topological_types.push_back(kTopTypeM);
                            bond_index.push_back(parent_index);
                            break;
                    }
                }
                else if(assembly_atom->GetNode()->GetNodeNeighbors().size() == 2)
                {
                    top_stack = stack.at(stack.size() - 1);
                    switch(visited_neighbors)
                    {
                        case 0:
                            stack.push_back(M2);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopTypeM);
                            bond_index.push_back(index + DEFAULT_DUMMY_ATOMS - 1);
                            break;
                        case 1:
                        {
                            stack.pop_back();
                            int parent_index = atom_index_stack.at(atom_index_stack.size() - 1);
                            atom_index_stack.pop_back();
                            stack.push_back(M1);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopTypeM);
                            bond_index.push_back(parent_index);
                            break;
                        }
                        case 2:
                            stack.pop_back();
                            int parent_index = atom_index_stack.at(atom_index_stack.size() - 1);
                            atom_index_stack.pop_back();
                            int top_atom_index_stack = atom_index_stack.at(0);
                            loops[top_atom_index_stack] = index + DEFAULT_DUMMY_ATOMS;
                            atom_index_stack.pop_back();
                            stack.pop_back();
                            topological_types.push_back(kTopTypeE);
                            bond_index.push_back(parent_index);
                            break;
                    }
                }
                else if(assembly_atom->GetNode()->GetNodeNeighbors().size() == 3)
                {
                    top_stack = stack.at(stack.size() - 1);
                    switch(visited_neighbors)
                    {
                        case 0:
                        {
                            stack.push_back(M3);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopTypeM);
                            bond_index.push_back(index + DEFAULT_DUMMY_ATOMS - 1);
                            break;
                        }
                        case 1:
                        {
                            stack.pop_back();
                            int parent_index = atom_index_stack.at(atom_index_stack.size() - 1);
                            atom_index_stack.pop_back();
                            stack.push_back(M2);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopTypeM);
                            bond_index.push_back(parent_index);
                            break;
                        }
                        case 2:
                        {
                            stack.pop_back();
                            int parent_index = atom_index_stack.at(atom_index_stack.size() - 1);
                            atom_index_stack.pop_back();
                            int top_atom_index_stack = atom_index_stack.at(0);
                            loops[top_atom_index_stack] = index + DEFAULT_DUMMY_ATOMS;
                            atom_index_stack.pop_back();
                            stack.pop_back();
                            stack.push_back(M1);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopTypeM);
                            bond_index.push_back(parent_index);
                            break;
                        }
                        case 3:
                        {
                            stack.pop_back();
                            atom_index_stack.pop_back();
                            topological_types.push_back(kTopTypeE);
                            bond_index.push_back(index + DEFAULT_DUMMY_ATOMS - 1);
                            break;
                        }
                    }
                }
                else if(assembly_atom->GetNode()->GetNodeNeighbors().size() == 4)
                {
                    top_stack = stack.at(stack.size() - 1);
                    switch(visited_neighbors)
                    {
                        case 0:
                        {
                            stack.push_back(M4);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopTypeM);
                            bond_index.push_back(index + DEFAULT_DUMMY_ATOMS - 1);
                            break;
                        }
                        case 1:
                        {
                            stack.pop_back();
                            int parent_index = atom_index_stack.at(atom_index_stack.size() - 1);
                            atom_index_stack.pop_back();
                            stack.push_back(M3);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopTypeM);
                            bond_index.push_back(parent_index);
                            break;
                        }
                        case 2:
                        {
                            stack.pop_back();
                            int parent_index = atom_index_stack.at(atom_index_stack.size() - 1);
                            atom_index_stack.pop_back();
                            int top_atom_index_stack = atom_index_stack.at(0);
                            loops[top_atom_index_stack] = index + DEFAULT_DUMMY_ATOMS;
                            atom_index_stack.pop_back();
                            stack.pop_back();
                            stack.push_back(M2);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopTypeM);
                            bond_index.push_back(parent_index);
                            break;
                        }
                        case 3:
                        {
                            stack.pop_back();
                            int parent_index = atom_index_stack.at(atom_index_stack.size() - 1);
                            atom_index_stack.pop_back();
                            stack.push_back(M1);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopTypeM);
                            bond_index.push_back(parent_index);
                            break;
                        }
                        case 4:
                        {
                            stack.pop_back();
                            int parent_index = atom_index_stack.at(atom_index_stack.size() - 1);
                            atom_index_stack.pop_back();
                            topological_types.push_back(kTopTypeE);
                            bond_index.push_back(parent_index);
                            break;
                        }
                    }
                }
                break;
            case M2:
                if(assembly_atom->GetNode()->GetNodeNeighbors().size() == 1)
                {
                    top_stack = stack.at(stack.size() - 1);
                    switch(visited_neighbors)
                    {
                        case 0:
                        {
                            stack.push_back(M1);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopTypeM);
                            bond_index.push_back(index + DEFAULT_DUMMY_ATOMS - 1);
                            break;
                        }
                        case 1:
                        {
                            stack.pop_back();
                            int parent_index = atom_index_stack.at(atom_index_stack.size() - 1);
                            atom_index_stack.pop_back();
                            stack.push_back(M1);
                            topological_types.push_back(kTopTypeE);
                            bond_index.push_back(parent_index);
                            break;
                        }
                    }
                }
                else if(assembly_atom->GetNode()->GetNodeNeighbors().size() == 2)
                {
                    top_stack = stack.at(stack.size() - 1);
                    switch(visited_neighbors)
                    {
                        case 0:
                        {
                            stack.push_back(M2);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopTypeM);
                            bond_index.push_back(index + DEFAULT_DUMMY_ATOMS - 1);
                            break;
                        }
                        case 1:
                        {
                            stack.pop_back();
                            int parent_index = atom_index_stack.at(atom_index_stack.size() - 1);
                            atom_index_stack.pop_back();
                            stack.push_back(M1);
                            stack.push_back(S1);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopTypeS);
                            bond_index.push_back(parent_index);
                            break;
                        }
                        case 2:
                        {
                            stack.pop_back();
                            int parent_index = atom_index_stack.at(atom_index_stack.size() - 1);
                            atom_index_stack.pop_back();
                            int top_atom_index_stack = atom_index_stack.at(0);
                            loops[top_atom_index_stack] = index + DEFAULT_DUMMY_ATOMS;
                            atom_index_stack.pop_back();
                            stack.pop_back();
                            stack.push_back(M1);
                            topological_types.push_back(kTopTypeE);
                            bond_index.push_back(parent_index);
                            break;
                        }
                    }
                }
                else if(assembly_atom->GetNode()->GetNodeNeighbors().size() == 3)
                {
                    top_stack = stack.at(stack.size() - 1);
                    switch(visited_neighbors)
                    {
                        case 0:
                        {
                            stack.push_back(M3);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopTypeM);
                            bond_index.push_back(index + DEFAULT_DUMMY_ATOMS - 1);
                            break;
                        }
                        case 1:
                        {
                            stack.pop_back();
                            int parent_index = atom_index_stack.at(atom_index_stack.size() - 1);
                            atom_index_stack.pop_back();
                            stack.push_back(M1);
                            stack.push_back(B2);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopTypeB);
                            bond_index.push_back(parent_index);
                            break;
                        }
                        case 2:
                        {
                            stack.pop_back();
                            int parent_index = atom_index_stack.at(atom_index_stack.size() - 1);
                            atom_index_stack.pop_back();
                            int top_atom_index_stack = atom_index_stack.at(0);
                            loops[top_atom_index_stack] = index + DEFAULT_DUMMY_ATOMS;
                            atom_index_stack.pop_back();
                            stack.pop_back();
                            stack.push_back(M1);
                            stack.push_back(B1);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopTypeS);
                            bond_index.push_back(parent_index);
                            break;
                        }
                        case 3:
                        {
                            stack.pop_back();
                            int parent_index = atom_index_stack.at(atom_index_stack.size() - 1);
                            atom_index_stack.pop_back();
                            stack.push_back(M1);
                            topological_types.push_back(kTopTypeE);
                            bond_index.push_back(parent_index);
                            break;
                        }
                    }
                }
                else if(assembly_atom->GetNode()->GetNodeNeighbors().size() == 4)
                {
                    top_stack = stack.at(stack.size() - 1);
                    switch(visited_neighbors)
                    {
                        case 0:
                        {
                            stack.push_back(M4);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopTypeM);
                            bond_index.push_back(index + DEFAULT_DUMMY_ATOMS - 1);
                            break;
                        }
                        case 1:
                        {
                            stack.pop_back();
                            int parent_index = atom_index_stack.at(atom_index_stack.size() - 1);
                            atom_index_stack.pop_back();
                            stack.push_back(M1);
                            stack.push_back(T3);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopType3);
                            bond_index.push_back(parent_index);
                            break;
                        }
                        case 2:
                        {
                            stack.pop_back();
                            int parent_index = atom_index_stack.at(atom_index_stack.size() - 1);
                            atom_index_stack.pop_back();
                            int top_atom_index_stack = atom_index_stack.at(0);
                            loops[top_atom_index_stack] = index + DEFAULT_DUMMY_ATOMS;
                            atom_index_stack.pop_back();
                            stack.pop_back();
                            stack.push_back(M1);
                            stack.push_back(T2);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopTypeB);
                            bond_index.push_back(parent_index);
                            break;
                        }
                        case 3:
                        {
                            stack.pop_back();
                            int parent_index = atom_index_stack.at(atom_index_stack.size() - 1);
                            atom_index_stack.pop_back();
                            stack.push_back(M1);
                            stack.push_back(T1);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopTypeS);
                            bond_index.push_back(parent_index);
                            break;
                        }
                        case 4:
                        {
                            stack.pop_back();
                            int parent_index = atom_index_stack.at(atom_index_stack.size() - 1);
                            atom_index_stack.pop_back();
                            stack.push_back(M1);
                            topological_types.push_back(kTopTypeE);
                            bond_index.push_back(parent_index);
                            break;
                        }
                    }
                }
                break;
            case M3:
                if(assembly_atom->GetNode()->GetNodeNeighbors().size() == 1)
                {
                    top_stack = stack.at(stack.size() - 1);
                    switch(visited_neighbors)
                    {
                        case 0:
                        {
                            stack.push_back(M1);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopTypeM);
                            bond_index.push_back(index + DEFAULT_DUMMY_ATOMS - 1);
                            break;
                        }
                        case 1:
                        {
                            stack.pop_back();
                            int parent_index = atom_index_stack.at(atom_index_stack.size() - 1);
                            atom_index_stack.pop_back();
                            stack.push_back(M2);
                            topological_types.push_back(kTopTypeE);
                            bond_index.push_back(parent_index);
                            break;
                        }
                    }
                }
                else if(assembly_atom->GetNode()->GetNodeNeighbors().size() == 2)
                {
                    top_stack = stack.at(stack.size() - 1);
                    switch(visited_neighbors)
                    {
                        case 0:
                        {
                            stack.push_back(M2);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopTypeM);
                            bond_index.push_back(index + DEFAULT_DUMMY_ATOMS - 1);
                            break;
                        }
                        case 1:
                        {
                            stack.pop_back();
                            int parent_index = atom_index_stack.at(atom_index_stack.size() - 1);
                            atom_index_stack.pop_back();
                            stack.push_back(M2);
                            stack.push_back(S1);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopTypeS);
                            bond_index.push_back(parent_index);
                            break;
                        }
                        case 2:
                        {
                            stack.pop_back();
                            int parent_index = atom_index_stack.at(atom_index_stack.size() - 1);
                            atom_index_stack.pop_back();
                            int top_atom_index_stack = atom_index_stack.at(0);
                            loops[top_atom_index_stack] = index + DEFAULT_DUMMY_ATOMS;
                            atom_index_stack.pop_back();
                            stack.pop_back();
                            stack.push_back(M2);
                            topological_types.push_back(kTopTypeE);
                            bond_index.push_back(parent_index);
                            break;
                        }
                    }
                }
                else if(assembly_atom->GetNode()->GetNodeNeighbors().size() == 3)
                {
                    top_stack = stack.at(stack.size() - 1);
                    switch(visited_neighbors)
                    {
                        case 0:
                        {
                            stack.push_back(M3);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopTypeM);
                            bond_index.push_back(index + DEFAULT_DUMMY_ATOMS - 1);
                            break;
                        }
                        case 1:
                        {
                            stack.pop_back();
                            int parent_index = atom_index_stack.at(atom_index_stack.size() - 1);
                            atom_index_stack.pop_back();
                            stack.push_back(M2);
                            stack.push_back(B2);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopTypeB);
                            bond_index.push_back(parent_index);
                            break;
                        }
                        case 2:
                        {
                            stack.pop_back();
                            int parent_index = atom_index_stack.at(atom_index_stack.size() - 1);
                            atom_index_stack.pop_back();
                            int top_atom_index_stack = atom_index_stack.at(0);
                            loops[top_atom_index_stack] = index + DEFAULT_DUMMY_ATOMS;
                            atom_index_stack.pop_back();
                            stack.pop_back();
                            stack.push_back(M2);
                            stack.push_back(B1);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopTypeS);
                            bond_index.push_back(parent_index);
                            break;
                        }
                        case 3:
                        {
                            stack.pop_back();
                            int parent_index = atom_index_stack.at(atom_index_stack.size() - 1);
                            atom_index_stack.pop_back();
                            stack.push_back(M2);
                            topological_types.push_back(kTopTypeE);
                            bond_index.push_back(parent_index);
                            break;
                        }
                    }
                }
                else if(assembly_atom->GetNode()->GetNodeNeighbors().size() == 4)
                {
                    top_stack = stack.at(stack.size() - 1);
                    switch(visited_neighbors)
                    {
                        case 0:
                        {
                            stack.push_back(M4);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopTypeM);
                            bond_index.push_back(index + DEFAULT_DUMMY_ATOMS - 1);
                            break;
                        }
                        case 1:
                        {
                            stack.pop_back();
                            int parent_index = atom_index_stack.at(atom_index_stack.size() - 1);
                            atom_index_stack.pop_back();
                            stack.push_back(M2);
                            stack.push_back(T3);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopType3);
                            bond_index.push_back(parent_index);
                            break;
                        }
                        case 2:
                        {
                            stack.pop_back();
                            int parent_index = atom_index_stack.at(atom_index_stack.size() - 1);
                            atom_index_stack.pop_back();
                            int top_atom_index_stack = atom_index_stack.at(0);
                            loops[top_atom_index_stack] = index + DEFAULT_DUMMY_ATOMS;
                            atom_index_stack.pop_back();
                            stack.pop_back();
                            stack.push_back(M2);
                            stack.push_back(T2);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopTypeB);
                            bond_index.push_back(parent_index);
                            break;
                        }
                        case 3:
                        {
                            stack.pop_back();
                            int parent_index = atom_index_stack.at(atom_index_stack.size() - 1);
                            atom_index_stack.pop_back();
                            stack.push_back(M2);
                            stack.push_back(T1);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopTypeS);
                            bond_index.push_back(parent_index);
                            break;
                        }
                        case 4:
                        {
                            stack.pop_back();
                            int parent_index = atom_index_stack.at(atom_index_stack.size() - 1);
                            atom_index_stack.pop_back();
                            stack.push_back(M2);
                            topological_types.push_back(kTopTypeE);
                            bond_index.push_back(parent_index);
                            break;
                        }
                    }
                }
                break;
            case M4:
                if(assembly_atom->GetNode()->GetNodeNeighbors().size() == 1)
                {
                    top_stack = stack.at(stack.size() - 1);
                    switch(visited_neighbors)
                    {
                        case 0:
                        {
                            stack.push_back(M1);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopTypeM);
                            bond_index.push_back(index + DEFAULT_DUMMY_ATOMS - 1);
                            break;
                        }
                        case 1:
                        {
                            stack.pop_back();
                            int parent_index = atom_index_stack.at(atom_index_stack.size() - 1);
                            atom_index_stack.pop_back();
                            stack.push_back(M3);
                            topological_types.push_back(kTopTypeE);
                            bond_index.push_back(parent_index);
                            break;
                        }
                    }
                }
                else if(assembly_atom->GetNode()->GetNodeNeighbors().size() == 2)
                {
                    top_stack = stack.at(stack.size() - 1);
                    switch(visited_neighbors)
                    {
                        case 0:
                        {
                            stack.push_back(M2);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopTypeM);
                            bond_index.push_back(index + DEFAULT_DUMMY_ATOMS - 1);
                            break;
                        }
                        case 1:
                        {
                            stack.pop_back();
                            int parent_index = atom_index_stack.at(atom_index_stack.size() - 1);
                            atom_index_stack.pop_back();
                            stack.push_back(M3);
                            stack.push_back(S1);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopTypeS);
                            bond_index.push_back(parent_index);
                            break;
                        }
                        case 2:
                        {
                            stack.pop_back();
                            int parent_index = atom_index_stack.at(atom_index_stack.size() - 1);
                            atom_index_stack.pop_back();
                            int top_atom_index_stack = atom_index_stack.at(0);
                            loops[top_atom_index_stack] = index + DEFAULT_DUMMY_ATOMS;
                            atom_index_stack.pop_back();
                            stack.pop_back();
                            stack.push_back(M3);
                            topological_types.push_back(kTopTypeE);
                            bond_index.push_back(parent_index);
                            break;
                        }
                    }
                }
                else if(assembly_atom->GetNode()->GetNodeNeighbors().size() == 3)
                {
                    top_stack = stack.at(stack.size() - 1);
                    switch(visited_neighbors)
                    {
                        case 0:
                        {
                            stack.push_back(M3);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopTypeM);
                            bond_index.push_back(index + DEFAULT_DUMMY_ATOMS - 1);
                            break;
                        }
                        case 1:
                        {
                            stack.pop_back();
                            int parent_index = atom_index_stack.at(atom_index_stack.size() - 1);
                            atom_index_stack.pop_back();
                            stack.push_back(M3);
                            stack.push_back(B2);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopTypeB);
                            bond_index.push_back(parent_index);
                            break;
                        }
                        case 2:
                        {
                            stack.pop_back();
                            int parent_index = atom_index_stack.at(atom_index_stack.size() - 1);
                            atom_index_stack.pop_back();
                            int top_atom_index_stack = atom_index_stack.at(0);
                            loops[top_atom_index_stack] = index + DEFAULT_DUMMY_ATOMS;
                            atom_index_stack.pop_back();
                            stack.pop_back();
                            stack.push_back(M3);
                            stack.push_back(B1);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopTypeS);
                            bond_index.push_back(parent_index);
                            break;
                        }
                        case 3:
                        {
                            stack.pop_back();
                            int parent_index = atom_index_stack.at(atom_index_stack.size() - 1);
                            atom_index_stack.pop_back();
                            stack.push_back(M3);
                            topological_types.push_back(kTopTypeE);
                            bond_index.push_back(parent_index);
                            break;
                        }
                    }
                }
                else if(assembly_atom->GetNode()->GetNodeNeighbors().size() == 4)
                {
                    top_stack = stack.at(stack.size() - 1);
                    switch(visited_neighbors)
                    {
                        case 0:
                        {
                            stack.push_back(M4);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopTypeM);
                            bond_index.push_back(index + DEFAULT_DUMMY_ATOMS - 1);
                            break;
                        }
                        case 1:
                        {
                            stack.pop_back();
                            int parent_index = atom_index_stack.at(atom_index_stack.size() - 1);
                            atom_index_stack.pop_back();
                            stack.push_back(M3);
                            stack.push_back(T3);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopType3);
                            bond_index.push_back(parent_index);
                            break;
                        }
                        case 2:
                        {
                            stack.pop_back();
                            int parent_index = atom_index_stack.at(atom_index_stack.size() - 1);
                            atom_index_stack.pop_back();
                            int top_atom_index_stack = atom_index_stack.at(0);
                            loops[top_atom_index_stack] = index + DEFAULT_DUMMY_ATOMS;
                            atom_index_stack.pop_back();
                            stack.pop_back();
                            stack.push_back(M3);
                            stack.push_back(T2);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopTypeB);
                            bond_index.push_back(parent_index);
                            break;
                        }
                        case 3:
                        {
                            stack.pop_back();
                            int parent_index = atom_index_stack.at(atom_index_stack.size() - 1);
                            atom_index_stack.pop_back();
                            stack.push_back(M3);
                            stack.push_back(T1);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopTypeS);
                            bond_index.push_back(parent_index);
                            break;
                        }
                        case 4:
                        {
                            stack.pop_back();
                            int parent_index = atom_index_stack.at(atom_index_stack.size() - 1);
                            atom_index_stack.pop_back();
                            stack.push_back(M3);
                            topological_types.push_back(kTopTypeE);
                            bond_index.push_back(parent_index);
                            break;
                        }
                    }
                }
                break;
            case S1:
                if(assembly_atom->GetNode()->GetNodeNeighbors().size() == 1)
                {
                    top_stack = stack.at(stack.size() - 1);
                    switch(visited_neighbors)
                    {
                        case 0:
                        {
                            stack.push_back(M1);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopTypeM);
                            bond_index.push_back(index + DEFAULT_DUMMY_ATOMS - 1);
                            break;
                        }
                        case 1:
                        {
                            stack.pop_back();
                            int parent_index = atom_index_stack.at(atom_index_stack.size() - 1);
                            atom_index_stack.pop_back();
                            topological_types.push_back(kTopTypeE);
                            bond_index.push_back(parent_index);
                            break;
                        }
                    }
                }
                else if(assembly_atom->GetNode()->GetNodeNeighbors().size() == 2)
                {
                    top_stack = stack.at(stack.size() - 1);
                    switch(visited_neighbors)
                    {
                        case 0:
                        {
                            stack.push_back(M2);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopTypeM);
                            bond_index.push_back(index + DEFAULT_DUMMY_ATOMS - 1);
                            break;
                        }
                        case 1:
                        {
                            stack.pop_back();
                            int parent_index = atom_index_stack.at(atom_index_stack.size() - 1);
                            atom_index_stack.pop_back();
                            stack.push_back(S1);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopTypeS);
                            bond_index.push_back(parent_index);
                            break;
                        }
                        case 2:
                        {
                            stack.pop_back();
                            int parent_index = atom_index_stack.at(atom_index_stack.size() - 1);
                            atom_index_stack.pop_back();
                            int top_atom_index_stack = atom_index_stack.at(0);
                            loops[top_atom_index_stack] = index + DEFAULT_DUMMY_ATOMS;
                            atom_index_stack.pop_back();
                            stack.pop_back();
                            topological_types.push_back(kTopTypeE);
                            bond_index.push_back(parent_index);
                            break;
                        }
                    }
                }
                else if(assembly_atom->GetNode()->GetNodeNeighbors().size() == 3)
                {
                    top_stack = stack.at(stack.size() - 1);
                    switch(visited_neighbors)
                    {
                        case 0:
                        {
                            stack.push_back(M3);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopTypeM);
                            bond_index.push_back(index + DEFAULT_DUMMY_ATOMS - 1);
                            break;
                        }
                        case 1:
                        {
                            stack.pop_back();
                            int parent_index = atom_index_stack.at(atom_index_stack.size() - 1);
                            atom_index_stack.pop_back();
                            stack.push_back(B2);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopTypeB);
                            bond_index.push_back(parent_index);
                            break;
                        }
                        case 2:
                        {
                            stack.pop_back();
                            int parent_index = atom_index_stack.at(atom_index_stack.size() - 1);
                            atom_index_stack.pop_back();
                            int top_atom_index_stack = atom_index_stack.at(0);
                            loops[top_atom_index_stack] = index + DEFAULT_DUMMY_ATOMS;
                            atom_index_stack.pop_back();
                            stack.pop_back();
                            stack.push_back(B1);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopTypeS);
                            bond_index.push_back(parent_index);
                            break;
                        }
                        case 3:
                        {
                            stack.pop_back();
                            int parent_index = atom_index_stack.at(atom_index_stack.size() - 1);
                            atom_index_stack.pop_back();
                            topological_types.push_back(kTopTypeE);
                            bond_index.push_back(parent_index);
                            break;
                        }
                    }

                }
                else if(assembly_atom->GetNode()->GetNodeNeighbors().size() == 4)
                {
                    top_stack = stack.at(stack.size() - 1);
                    switch(visited_neighbors)
                    {
                        case 0:
                        {
                            stack.push_back(M4);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopTypeM);
                            bond_index.push_back(index + DEFAULT_DUMMY_ATOMS - 1);
                            break;
                        }
                        case 1:
                        {
                            stack.pop_back();
                            int parent_index = atom_index_stack.at(atom_index_stack.size() - 1);
                            atom_index_stack.pop_back();
                            stack.push_back(T3);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopType3);
                            bond_index.push_back(parent_index);
                            break;
                        }
                        case 2:
                        {
                            stack.pop_back();
                            int parent_index = atom_index_stack.at(atom_index_stack.size() - 1);
                            atom_index_stack.pop_back();
                            int top_atom_index_stack = atom_index_stack.at(0);
                            loops[top_atom_index_stack] = index + DEFAULT_DUMMY_ATOMS;
                            atom_index_stack.pop_back();
                            stack.pop_back();
                            stack.push_back(T2);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopTypeB);
                            bond_index.push_back(parent_index);
                            break;
                        }
                        case 3:
                        {
                            stack.pop_back();
                            int parent_index = atom_index_stack.at(atom_index_stack.size() - 1);
                            atom_index_stack.pop_back();
                            stack.push_back(T1);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopTypeS);
                            bond_index.push_back(parent_index);
                            break;
                        }
                        case 4:
                        {
                            stack.pop_back();
                            int parent_index = atom_index_stack.at(atom_index_stack.size() - 1);
                            atom_index_stack.pop_back();
                            topological_types.push_back(kTopTypeE);
                            bond_index.push_back(parent_index);
                            break;
                        }
                    }
                }
                break;
            case B1:
                if(assembly_atom->GetNode()->GetNodeNeighbors().size() == 1)
                {
                    top_stack = stack.at(stack.size() - 1);
                    switch(visited_neighbors)
                    {
                        case 0:
                        {
                            stack.push_back(M1);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopTypeM);
                            bond_index.push_back(index + DEFAULT_DUMMY_ATOMS - 1);
                            break;
                        }
                        case 1:
                        {
                            stack.pop_back();
                            int parent_index = atom_index_stack.at(atom_index_stack.size() - 1);
                            atom_index_stack.pop_back();
                            topological_types.push_back(kTopTypeE);
                            bond_index.push_back(parent_index);
                            break;
                        }
                    }
                }
                else if(assembly_atom->GetNode()->GetNodeNeighbors().size() == 2)
                {
                    top_stack = stack.at(stack.size() - 1);
                    switch(visited_neighbors)
                    {
                        case 0:
                        {
                            stack.push_back(M2);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopTypeM);
                            bond_index.push_back(index + DEFAULT_DUMMY_ATOMS - 1);
                            break;
                        }
                        case 1:
                        {
                            stack.pop_back();
                            int parent_index = atom_index_stack.at(atom_index_stack.size() - 1);
                            atom_index_stack.pop_back();
                            stack.push_back(S1);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopTypeS);
                            bond_index.push_back(parent_index);
                            break;
                        }
                        case 2:
                        {
                            stack.pop_back();
                            int parent_index = atom_index_stack.at(atom_index_stack.size() - 1);
                            atom_index_stack.pop_back();
                            int top_atom_index_stack = atom_index_stack.at(0);
                            loops[top_atom_index_stack] = index + DEFAULT_DUMMY_ATOMS;
                            atom_index_stack.pop_back();
                            stack.pop_back();
                            topological_types.push_back(kTopTypeE);
                            bond_index.push_back(parent_index);
                            break;
                        }
                    }
                }
                else if(assembly_atom->GetNode()->GetNodeNeighbors().size() == 3)
                {
                    top_stack = stack.at(stack.size() - 1);
                    switch(visited_neighbors)
                    {
                        case 0:
                        {
                            stack.push_back(M3);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopTypeM);
                            bond_index.push_back(index + DEFAULT_DUMMY_ATOMS - 1);
                            break;
                        }
                        case 1:
                        {
                            stack.pop_back();
                            int parent_index = atom_index_stack.at(atom_index_stack.size() - 1);
                            atom_index_stack.pop_back();
                            stack.push_back(B2);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopTypeB);
                            bond_index.push_back(parent_index);
                            break;
                        }
                        case 2:
                        {
                            stack.pop_back();
                            int parent_index = atom_index_stack.at(atom_index_stack.size() - 1);
                            atom_index_stack.pop_back();
                            int top_atom_index_stack = atom_index_stack.at(0);
                            loops[top_atom_index_stack] = index + DEFAULT_DUMMY_ATOMS;
                            atom_index_stack.pop_back();
                            stack.pop_back();
                            stack.push_back(B1);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopTypeS);
                            bond_index.push_back(parent_index);
                            break;
                        }
                        case 3:
                        {
                            stack.pop_back();
                            int parent_index = atom_index_stack.at(atom_index_stack.size() - 1);
                            atom_index_stack.pop_back();
                            topological_types.push_back(kTopTypeE);
                            bond_index.push_back(parent_index);
                            break;
                        }
                    }
                }
                else if(assembly_atom->GetNode()->GetNodeNeighbors().size() == 4)
                {
                    top_stack = stack.at(stack.size() - 1);
                    switch(visited_neighbors)
                    {
                        case 0:
                        {
                            stack.push_back(M4);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopTypeM);
                            bond_index.push_back(index + DEFAULT_DUMMY_ATOMS - 1);
                            break;
                        }
                        case 1:
                        {
                            stack.pop_back();
                            int parent_index = atom_index_stack.at(atom_index_stack.size() - 1);
                            atom_index_stack.pop_back();
                            stack.push_back(T3);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopType3);
                            bond_index.push_back(parent_index);
                            break;
                        }
                        case 2:
                        {
                            stack.pop_back();
                            int parent_index = atom_index_stack.at(atom_index_stack.size() - 1);
                            atom_index_stack.pop_back();
                            int top_atom_index_stack = atom_index_stack.at(0);
                            loops[top_atom_index_stack] = index + DEFAULT_DUMMY_ATOMS;
                            atom_index_stack.pop_back();
                            stack.pop_back();
                            stack.push_back(T2);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopTypeB);
                            bond_index.push_back(parent_index);
                            break;
                        }
                        case 3:
                        {
                            stack.pop_back();
                            int parent_index = atom_index_stack.at(atom_index_stack.size() - 1);
                            atom_index_stack.pop_back();
                            stack.push_back(T1);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopTypeS);
                            bond_index.push_back(parent_index);
                            break;
                        }
                        case 4:
                        {
                            stack.pop_back();
                            int parent_index = atom_index_stack.at(atom_index_stack.size() - 1);
                            atom_index_stack.pop_back();
                            topological_types.push_back(kTopTypeE);
                            bond_index.push_back(parent_index);
                            break;
                        }
                    }
                }
                break;
            case B2:
                if(assembly_atom->GetNode()->GetNodeNeighbors().size() == 1)
                {
                    top_stack = stack.at(stack.size() - 1);
                    switch(visited_neighbors)
                    {
                        case 0:
                        {
                            stack.push_back(M1);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopTypeM);
                            bond_index.push_back(index + DEFAULT_DUMMY_ATOMS - 1);
                            break;
                        }
                        case 1:
                        {
                            stack.pop_back();
                            int parent_index = atom_index_stack.at(atom_index_stack.size() - 1);
                            atom_index_stack.pop_back();
                            stack.push_back(B1);
                            topological_types.push_back(kTopTypeE);
                            bond_index.push_back(parent_index);
                            break;
                        }
                    }
                }
                else if(assembly_atom->GetNode()->GetNodeNeighbors().size() == 2)
                {
                    top_stack = stack.at(stack.size() - 1);
                    switch(visited_neighbors)
                    {
                        case 0:
                        {
                            stack.push_back(M2);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopTypeM);
                            bond_index.push_back(index + DEFAULT_DUMMY_ATOMS - 1);
                            break;
                        }
                        case 1:
                        {
                            stack.pop_back();
                            int parent_index = atom_index_stack.at(atom_index_stack.size() - 1);
                            atom_index_stack.pop_back();
                            stack.push_back(B1);
                            stack.push_back(S1);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopTypeS);
                            bond_index.push_back(parent_index);
                            break;
                        }
                        case 2:
                        {
                            stack.pop_back();
                            int parent_index = atom_index_stack.at(atom_index_stack.size() - 1);
                            atom_index_stack.pop_back();
                            int top_atom_index_stack = atom_index_stack.at(0);
                            loops[top_atom_index_stack] = index + DEFAULT_DUMMY_ATOMS;
                            atom_index_stack.pop_back();
                            stack.pop_back();
                            stack.push_back(B1);
                            topological_types.push_back(kTopTypeE);
                            bond_index.push_back(parent_index);
                            break;
                        }
                    }
                }
                else if(assembly_atom->GetNode()->GetNodeNeighbors().size() == 3)
                {
                    top_stack = stack.at(stack.size() - 1);
                    switch(visited_neighbors)
                    {
                        case 0:
                        {
                            stack.push_back(M3);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopTypeM);
                            bond_index.push_back(index + DEFAULT_DUMMY_ATOMS - 1);
                            break;
                        }
                        case 1:
                        {
                            stack.pop_back();
                            int parent_index = atom_index_stack.at(atom_index_stack.size() - 1);
                            atom_index_stack.pop_back();
                            stack.push_back(B1);
                            stack.push_back(B2);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopTypeB);
                            bond_index.push_back(parent_index);
                            break;
                        }
                        case 2:
                        {
                            stack.pop_back();
                            int parent_index = atom_index_stack.at(atom_index_stack.size() - 1);
                            atom_index_stack.pop_back();
                            int top_atom_index_stack = atom_index_stack.at(0);
                            loops[top_atom_index_stack] = index + DEFAULT_DUMMY_ATOMS;
                            atom_index_stack.pop_back();
                            stack.pop_back();
                            stack.push_back(B1);
                            stack.push_back(B1);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopTypeS);
                            bond_index.push_back(parent_index);
                            break;
                        }
                        case 3:
                        {
                            stack.pop_back();
                            int parent_index = atom_index_stack.at(atom_index_stack.size() - 1);
                            atom_index_stack.pop_back();
                            stack.push_back(B1);
                            topological_types.push_back(kTopTypeE);
                            bond_index.push_back(parent_index);
                            break;
                        }
                    }
                }
                else if(assembly_atom->GetNode()->GetNodeNeighbors().size() == 4)
                {
                    top_stack = stack.at(stack.size() - 1);
                    switch(visited_neighbors)
                    {
                        case 0:
                        {
                            stack.push_back(M4);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopTypeM);
                            bond_index.push_back(index + DEFAULT_DUMMY_ATOMS - 1);
                            break;
                        }
                        case 1:
                        {
                            stack.pop_back();
                            int parent_index = atom_index_stack.at(atom_index_stack.size() - 1);
                            atom_index_stack.pop_back();
                            stack.push_back(B1);
                            stack.push_back(T3);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopType3);
                            bond_index.push_back(parent_index);
                            break;
                        }
                        case 2:
                        {
                            stack.pop_back();
                            int parent_index = atom_index_stack.at(atom_index_stack.size() - 1);
                            atom_index_stack.pop_back();
                            int top_atom_index_stack = atom_index_stack.at(0);
                            loops[top_atom_index_stack] = index + DEFAULT_DUMMY_ATOMS;
                            atom_index_stack.pop_back();
                            stack.pop_back();
                            stack.push_back(B1);
                            stack.push_back(T2);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopTypeB);
                            bond_index.push_back(parent_index);
                            break;
                        }
                        case 3:
                        {
                            stack.pop_back();
                            int parent_index = atom_index_stack.at(atom_index_stack.size() - 1);
                            atom_index_stack.pop_back();
                            stack.push_back(B1);
                            stack.push_back(T1);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopTypeS);
                            bond_index.push_back(parent_index);
                            break;
                        }
                        case 4:
                        {
                            stack.pop_back();
                            int parent_index = atom_index_stack.at(atom_index_stack.size() - 1);
                            atom_index_stack.pop_back();
                            stack.push_back(B1);
                            topological_types.push_back(kTopTypeE);
                            bond_index.push_back(parent_index);
                            break;
                        }
                    }
                }
                break;
            case T1:
                if(assembly_atom->GetNode()->GetNodeNeighbors().size() == 1)
                {
                    top_stack = stack.at(stack.size() - 1);
                    switch(visited_neighbors)
                    {
                        case 0:
                        {
                            stack.push_back(M1);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopTypeM);
                            bond_index.push_back(index + DEFAULT_DUMMY_ATOMS - 1);
                            break;
                        }
                        case 1:
                        {
                            stack.pop_back();
                            int parent_index = atom_index_stack.at(atom_index_stack.size() - 1);
                            atom_index_stack.pop_back();
                            topological_types.push_back(kTopTypeE);
                            bond_index.push_back(parent_index);
                            break;
                        }
                    }
                }
                else if(assembly_atom->GetNode()->GetNodeNeighbors().size() == 2)
                {
                    top_stack = stack.at(stack.size() - 1);
                    switch(visited_neighbors)
                    {
                        case 0:
                        {
                            stack.push_back(M2);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopTypeM);
                            bond_index.push_back(index + DEFAULT_DUMMY_ATOMS - 1);
                            break;
                        }
                        case 1:
                        {
                            stack.pop_back();
                            int parent_index = atom_index_stack.at(atom_index_stack.size() - 1);
                            atom_index_stack.pop_back();
                            stack.push_back(S1);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopTypeS);
                            bond_index.push_back(parent_index);
                            break;
                        }
                        case 2:
                        {
                            stack.pop_back();
                            int parent_index = atom_index_stack.at(atom_index_stack.size() - 1);
                            atom_index_stack.pop_back();
                            int top_atom_index_stack = atom_index_stack.at(0);
                            loops[top_atom_index_stack] = index + DEFAULT_DUMMY_ATOMS;
                            atom_index_stack.pop_back();
                            stack.pop_back();
                            topological_types.push_back(kTopTypeE);
                            bond_index.push_back(parent_index);
                            break;
                        }
                    }
                }
                else if(assembly_atom->GetNode()->GetNodeNeighbors().size() == 3)
                {
                    top_stack = stack.at(stack.size() - 1);
                    switch(visited_neighbors)
                    {
                        case 0:
                        {
                            stack.push_back(M3);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopTypeM);
                            bond_index.push_back(index + DEFAULT_DUMMY_ATOMS - 1);
                            break;
                        }
                        case 1:
                        {
                            stack.pop_back();
                            int parent_index = atom_index_stack.at(atom_index_stack.size() - 1);
                            atom_index_stack.pop_back();
                            stack.push_back(B2);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopTypeB);
                            bond_index.push_back(parent_index);
                            break;
                        }
                        case 2:
                        {
                            stack.pop_back();
                            int parent_index = atom_index_stack.at(atom_index_stack.size() - 1);
                            atom_index_stack.pop_back();
                            int top_atom_index_stack = atom_index_stack.at(0);
                            loops[top_atom_index_stack] = index + DEFAULT_DUMMY_ATOMS;
                            atom_index_stack.pop_back();
                            stack.pop_back();
                            stack.push_back(B1);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopTypeS);
                            bond_index.push_back(parent_index);
                            break;
                        }
                        case 3:
                        {
                            stack.pop_back();
                            int parent_index = atom_index_stack.at(atom_index_stack.size() - 1);
                            atom_index_stack.pop_back();
                            topological_types.push_back(kTopTypeE);
                            bond_index.push_back(parent_index);
                            break;
                        }
                    }
                }
                else if(assembly_atom->GetNode()->GetNodeNeighbors().size() == 4)
                {
                    top_stack = stack.at(stack.size() - 1);
                    switch(visited_neighbors)
                    {
                        case 0:
                        {
                            stack.push_back(M4);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopTypeM);
                            bond_index.push_back(index + DEFAULT_DUMMY_ATOMS - 1);
                            break;
                        }
                        case 1:
                        {
                            stack.pop_back();
                            int parent_index = atom_index_stack.at(atom_index_stack.size() - 1);
                            atom_index_stack.pop_back();
                            stack.push_back(T3);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopType3);
                            bond_index.push_back(parent_index);
                            break;
                        }
                        case 2:
                        {
                            stack.pop_back();
                            int parent_index = atom_index_stack.at(atom_index_stack.size() - 1);
                            atom_index_stack.pop_back();
                            int top_atom_index_stack = atom_index_stack.at(0);
                            loops[top_atom_index_stack] = index + DEFAULT_DUMMY_ATOMS;
                            atom_index_stack.pop_back();
                            stack.pop_back();
                            stack.push_back(T2);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopTypeB);
                            bond_index.push_back(parent_index);
                            break;
                        }
                        case 3:
                        {
                            stack.pop_back();
                            int parent_index = atom_index_stack.at(atom_index_stack.size() - 1);
                            atom_index_stack.pop_back();
                            stack.push_back(T1);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopTypeS);
                            bond_index.push_back(parent_index);
                            break;
                        }
                        case 4:
                        {
                            stack.pop_back();
                            int parent_index = atom_index_stack.at(atom_index_stack.size() - 1);
                            atom_index_stack.pop_back();
                            topological_types.push_back(kTopTypeE);
                            bond_index.push_back(parent_index);
                            break;
                        }
                    }
                }
                break;
            case T2:
                if(assembly_atom->GetNode()->GetNodeNeighbors().size() == 1)
                {
                    top_stack = stack.at(stack.size() - 1);
                    switch(visited_neighbors)
                    {
                        case 0:
                        {
                            stack.push_back(M1);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopTypeM);
                            bond_index.push_back(index + DEFAULT_DUMMY_ATOMS - 1);
                            break;
                        }
                        case 1:
                        {
                            stack.pop_back();
                            int parent_index = atom_index_stack.at(atom_index_stack.size() - 1);
                            atom_index_stack.pop_back();
                            stack.push_back(T1);
                            topological_types.push_back(kTopTypeE);
                            bond_index.push_back(parent_index);
                            break;
                        }
                    }
                }
                else if(assembly_atom->GetNode()->GetNodeNeighbors().size() == 2)
                {
                    top_stack = stack.at(stack.size() - 1);
                    switch(visited_neighbors)
                    {
                        case 0:
                        {
                            stack.push_back(M2);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopTypeM);
                            bond_index.push_back(index + DEFAULT_DUMMY_ATOMS - 1);
                            break;
                        }
                        case 1:
                        {
                            stack.pop_back();
                            int parent_index = atom_index_stack.at(atom_index_stack.size() - 1);
                            atom_index_stack.pop_back();
                            stack.push_back(T1);
                            stack.push_back(S1);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopTypeS);
                            bond_index.push_back(parent_index);
                            break;
                        }
                        case 2:
                        {
                            stack.pop_back();
                            int parent_index = atom_index_stack.at(atom_index_stack.size() - 1);
                            atom_index_stack.pop_back();
                            int top_atom_index_stack = atom_index_stack.at(0);
                            loops[top_atom_index_stack] = index + DEFAULT_DUMMY_ATOMS;
                            atom_index_stack.pop_back();
                            stack.pop_back();
                            stack.push_back(T1);
                            topological_types.push_back(kTopTypeE);
                            bond_index.push_back(parent_index);
                            break;
                        }
                    }
                }
                else if(assembly_atom->GetNode()->GetNodeNeighbors().size() == 3)
                {
                    top_stack = stack.at(stack.size() - 1);
                    switch(visited_neighbors)
                    {
                        case 0:
                        {
                            stack.push_back(M3);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopTypeM);
                            bond_index.push_back(index + DEFAULT_DUMMY_ATOMS - 1);
                            break;
                        }
                        case 1:
                        {
                            stack.pop_back();
                            int parent_index = atom_index_stack.at(atom_index_stack.size() - 1);
                            atom_index_stack.pop_back();
                            stack.push_back(T1);
                            stack.push_back(B2);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopTypeB);
                            bond_index.push_back(parent_index);
                            break;
                        }
                        case 2:
                        {
                            stack.pop_back();
                            int parent_index = atom_index_stack.at(atom_index_stack.size() - 1);
                            atom_index_stack.pop_back();
                            int top_atom_index_stack = atom_index_stack.at(0);
                            loops[top_atom_index_stack] = index + DEFAULT_DUMMY_ATOMS;
                            atom_index_stack.pop_back();
                            stack.pop_back();
                            stack.push_back(T1);
                            stack.push_back(B1);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopTypeS);
                            bond_index.push_back(parent_index);
                            break;
                        }
                        case 3:
                        {
                            stack.pop_back();
                            int parent_index = atom_index_stack.at(atom_index_stack.size() - 1);
                            atom_index_stack.pop_back();
                            stack.push_back(T1);
                            topological_types.push_back(kTopTypeE);
                            bond_index.push_back(parent_index);
                            break;
                        }
                    }
                }
                else if(assembly_atom->GetNode()->GetNodeNeighbors().size() == 4)
                {
                    top_stack = stack.at(stack.size() - 1);
                    switch(visited_neighbors)
                    {
                        case 0:
                        {
                            stack.push_back(M4);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopTypeM);
                            bond_index.push_back(index + DEFAULT_DUMMY_ATOMS - 1);
                            break;
                        }
                        case 1:
                        {
                            stack.pop_back();
                            int parent_index = atom_index_stack.at(atom_index_stack.size() - 1);
                            atom_index_stack.pop_back();
                            stack.push_back(T1);
                            stack.push_back(T3);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopType3);
                            bond_index.push_back(parent_index);
                            break;
                        }
                        case 2:
                        {
                            stack.pop_back();
                            int parent_index = atom_index_stack.at(atom_index_stack.size() - 1);
                            atom_index_stack.pop_back();
                            int top_atom_index_stack = atom_index_stack.at(0);
                            loops[top_atom_index_stack] = index + DEFAULT_DUMMY_ATOMS;
                            atom_index_stack.pop_back();
                            stack.pop_back();
                            stack.push_back(T1);
                            stack.push_back(T2);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopTypeB);
                            bond_index.push_back(parent_index);
                            break;
                        }
                        case 3:
                        {
                            stack.pop_back();
                            int parent_index = atom_index_stack.at(atom_index_stack.size() - 1);
                            atom_index_stack.pop_back();
                            stack.push_back(T1);
                            stack.push_back(T1);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopTypeS);
                            bond_index.push_back(parent_index);
                            break;
                        }
                        case 4:
                        {
                            stack.pop_back();
                            int parent_index = atom_index_stack.at(atom_index_stack.size() - 1);
                            atom_index_stack.pop_back();
                            stack.push_back(T1);
                            topological_types.push_back(kTopTypeE);
                            bond_index.push_back(parent_index);
                            break;
                        }
                    }
                }
                break;
            case T3:
                if(assembly_atom->GetNode()->GetNodeNeighbors().size() == 1)
                {
                    top_stack = stack.at(stack.size() - 1);
                    switch(visited_neighbors)
                    {
                        case 0:
                        {
                            stack.push_back(M1);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopTypeM);
                            bond_index.push_back(index + DEFAULT_DUMMY_ATOMS - 1);
                            break;
                        }
                        case 1:
                        {
                            stack.pop_back();
                            int parent_index = atom_index_stack.at(atom_index_stack.size() - 1);
                            atom_index_stack.pop_back();
                            stack.push_back(T2);
                            topological_types.push_back(kTopTypeE);
                            bond_index.push_back(parent_index);
                            break;
                        }
                    }
                }
                else if(assembly_atom->GetNode()->GetNodeNeighbors().size() == 2)
                {
                    top_stack = stack.at(stack.size() - 1);
                    switch(visited_neighbors)
                    {
                        case 0:
                        {
                            stack.push_back(M2);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopTypeM);
                            bond_index.push_back(index + DEFAULT_DUMMY_ATOMS - 1);
                            break;
                        }
                        case 1:
                        {
                            stack.pop_back();
                            int parent_index = atom_index_stack.at(atom_index_stack.size() - 1);
                            atom_index_stack.pop_back();
                            stack.push_back(T2);
                            stack.push_back(S1);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopTypeS);
                            bond_index.push_back(parent_index);
                            break;
                        }
                        case 2:
                        {
                            stack.pop_back();
                            int parent_index = atom_index_stack.at(atom_index_stack.size() - 1);
                            atom_index_stack.pop_back();
                            int top_atom_index_stack = atom_index_stack.at(0);
                            loops[top_atom_index_stack] = index + DEFAULT_DUMMY_ATOMS;
                            atom_index_stack.pop_back();
                            stack.pop_back();
                            stack.push_back(T2);
                            topological_types.push_back(kTopTypeE);
                            bond_index.push_back(parent_index);
                            break;
                        }
                    }
                }
                else if(assembly_atom->GetNode()->GetNodeNeighbors().size() == 3)
                {
                    top_stack = stack.at(stack.size() - 1);
                    switch(visited_neighbors)
                    {
                        case 0:
                        {
                            stack.push_back(M3);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopTypeM);
                            bond_index.push_back(index + DEFAULT_DUMMY_ATOMS - 1);
                            break;
                        }
                        case 1:
                        {
                            stack.pop_back();
                            int parent_index = atom_index_stack.at(atom_index_stack.size() - 1);
                            atom_index_stack.pop_back();
                            stack.push_back(T2);
                            stack.push_back(B2);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopTypeB);
                            bond_index.push_back(parent_index);
                            break;
                        }
                        case 2:
                        {
                            stack.pop_back();
                            int parent_index = atom_index_stack.at(atom_index_stack.size() - 1);
                            atom_index_stack.pop_back();
                            int top_atom_index_stack = atom_index_stack.at(0);
                            loops[top_atom_index_stack] = index + DEFAULT_DUMMY_ATOMS;
                            atom_index_stack.pop_back();
                            stack.pop_back();
                            stack.push_back(T2);
                            stack.push_back(B1);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopTypeS);
                            bond_index.push_back(parent_index);
                            break;
                        }
                        case 3:
                        {
                            stack.pop_back();
                            int parent_index = atom_index_stack.at(atom_index_stack.size() - 1);
                            atom_index_stack.pop_back();
                            stack.push_back(T2);
                            topological_types.push_back(kTopTypeE);
                            bond_index.push_back(parent_index);
                            break;
                        }
                    }
                }
                else if(assembly_atom->GetNode()->GetNodeNeighbors().size() == 4)
                {
                    top_stack = stack.at(stack.size() - 1);
                    switch(visited_neighbors)
                    {
                        case 0:
                        {
                            stack.push_back(M4);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopTypeM);
                            bond_index.push_back(index + DEFAULT_DUMMY_ATOMS - 1);
                            break;
                        }
                        case 1:
                        {
                            stack.pop_back();
                            int parent_index = atom_index_stack.at(atom_index_stack.size() - 1);
                            atom_index_stack.pop_back();
                            stack.push_back(T2);
                            stack.push_back(T3);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopType3);
                            bond_index.push_back(parent_index);
                            break;
                        }
                        case 2:
                        {
                            stack.pop_back();
                            int parent_index = atom_index_stack.at(atom_index_stack.size() - 1);
                            atom_index_stack.pop_back();
                            int top_atom_index_stack = atom_index_stack.at(0);
                            loops[top_atom_index_stack] = index + DEFAULT_DUMMY_ATOMS;
                            atom_index_stack.pop_back();
                            stack.pop_back();
                            stack.push_back(T2);
                            stack.push_back(T2);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopTypeB);
                            bond_index.push_back(parent_index);
                            break;
                        }
                        case 3:
                        {
                            stack.pop_back();
                            int parent_index = atom_index_stack.at(atom_index_stack.size() - 1);
                            atom_index_stack.pop_back();
                            stack.push_back(T2);
                            stack.push_back(T1);
                            atom_index_stack.push_back(index + DEFAULT_DUMMY_ATOMS);
                            topological_types.push_back(kTopTypeS);
                            bond_index.push_back(parent_index);
                            break;
                        }
                        case 4:
                        {
                            stack.pop_back();
                            int parent_index = atom_index_stack.at(atom_index_stack.size() - 1);
                            atom_index_stack.pop_back();
                            stack.push_back(T2);
                            topological_types.push_back(kTopTypeE);
                            bond_index.push_back(parent_index);
                            break;
                        }
                    }
                }
                break;
        }
        visited_atom_name.push_back(assembly_atom->GetName());
    }
    return topological_types;

}

TopologyFile* Assembly::BuildTopologyFileStructureFromAssembly(string parameter_file_path)
{
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
    topology_file->SetNumberOfBondTypes(this->CountNumberOfBondTypes());
    topology_file->SetNumberOfAngleTypes(this->CountNumberOfAngleTypes());
    topology_file->SetNumberOfDihedralTypes(this->CountNumberOfDihedralTypes(parameter_file_path));
    //    topology_file->SetNumberOfAtomTypesInParameterFile();
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
        ParameterFile* parameter_file = new ParameterFile(parameter_file_path);
        ParameterFileSpace::ParameterFile::BondMap bonds = parameter_file->GetBonds();
        ParameterFileSpace::ParameterFile::AngleMap angles = parameter_file->GetAngles();
        ParameterFileSpace::ParameterFile::DihedralMap dihedrals = parameter_file->GetDihedrals();
        ParameterFileSpace::ParameterFile::AtomTypeMap atom_types_map = parameter_file->GetAtomTypes();
        for(AtomVector::iterator it1 = assembly_atoms.begin(); it1 != assembly_atoms.end(); it1++)
        {
            Atom* assembly_atom = (*it1);
            stringstream key1;
            key1 << assembly_atom->GetId();
            TopologyAtom* topology_atom = new TopologyAtom();
            topology_atom->SetAtomName(assembly_atom->GetName());
            topology_atom->SetAtomCharge(assembly_atom->GetCharge());
            topology_atom->SetAtomicNumber(iNotSet);
            topology_atom->SetAtomMass(assembly_atom->GetMass());
            //            topology_atom->SetNumberOfExcludedAtomsForEachAtom();
            topology_atom->SetResidueName(assembly_residue->GetName());
            topology_atom->SetType(assembly_atom->GetAtomType());
            topology_atom->SetTreeChainClasification("0");
            topology_atom->SetRadii(dNotSet);
            topology_atom->SetScreen(dNotSet);
            topology_atom->SetIndex(atom_counter);

            topology_residue->AddAtom(topology_atom);
            atom_counter++;

            ///Pairs
            for(AtomVector::iterator it2 = assembly_atoms.begin(); it2 != assembly_atoms.end(); it2++)
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
                        find(inserted_pairs.begin(), inserted_pairs.end(), reverse_sss.str()) == inserted_pairs.end()  )
                {
                    TopologyAtomPair* topology_atom_pair = new TopologyAtomPair();
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
            }

            ///Bond Types, Bonds
            AtomNode* atom_node = assembly_atom->GetNode();
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
        topology_assembly->AddResidue(topology_residue);
    }



    //    for(ResidueVector::iterator it = assembly_residues.begin(); it != assembly_residues.end(); it++)
    //    {
    //        Residue* assembly_residue = *it;
    //        AtomVector assembly_atoms = assembly_residue->GetAtoms();
    //        for(AtomVector::iterator it1 = assembly_atoms.begin(); it1 != assembly_atoms.end(); it1++)
    //        {
    //            Atom* assembly_atom = (*it1);


    //        }
    //    }
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
            cout << atom_pair_type.at(0) << "-" << atom_pair_type.at(1) << " bond type does not exist in the parameter files" << endl;
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
            cout << atom_pair_type.at(0) << "-" << atom_pair_type.at(1) << " bond type does not exist in the parameter files" << endl;
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
            cout << angle_type.at(0) << "-" << angle_type.at(1) << "-" << angle_type.at(2) << " angle type does not exist in the parameter files" << endl;
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
            cout << angle_type.at(0) << "-" << angle_type.at(1) << "-" << angle_type.at(2) << " angle type does not exist in the parameter files" << endl;
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
                break;
            }
        }
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
                    break;
                }
            }
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
                                || (neighbor->GetName().substr(0,1).compare("H") == 0 ||
                                    (neighbor->GetName().substr(1,1).compare("H") == 0 && isdigit(ConvertString<char>(neighbor->GetName().substr(0,1)))))
                                ||(neighbor_of_neighbor->GetName().substr(0,1).compare("H") == 0 ||
                                   (neighbor_of_neighbor->GetName().substr(1,1).compare("H") == 0 && isdigit(ConvertString<char>(neighbor_of_neighbor->GetName().substr(0,1)))))
                                ||(neighbor_of_neighbor_of_neighbor->GetName().substr(0,1).compare("H") == 0 ||
                                   (neighbor_of_neighbor_of_neighbor->GetName().substr(1,1).compare("H") == 0 && isdigit(ConvertString<char>(neighbor_of_neighbor_of_neighbor->GetName().substr(0,1))))))
                            topology_dihedral->SetIncludingHydrogen(true);
                        else
                            topology_dihedral->SetIncludingHydrogen(false);

                        if(permutation_index == 0 || (permutation_index >= 6 && permutation_index <= 9) || (permutation_index >= 30 && permutation_index <= 35) || (permutation_index >= 66 && permutation_index <= 69))
                        {
                            topology_dihedral->SetResidueNames(residue_names1);
                            topology_dihedral->SetDihedrals(dihedral_atom_names1);
                        }
                        if(permutation_index == 2 || (permutation_index >= 14 && permutation_index <= 17) || (permutation_index >= 42 && permutation_index <= 47) || (permutation_index >= 74 && permutation_index <= 77))
                        {
                            topology_dihedral->SetResidueNames(residue_names2);
                            topology_dihedral->SetDihedrals(dihedral_atom_names2);
                        }
                        if(permutation_index == 4 || (permutation_index >= 22 && permutation_index <= 25) || (permutation_index >= 54 && permutation_index <= 59) || (permutation_index >= 82 && permutation_index <= 85))
                        {
                            topology_dihedral->SetResidueNames(residue_names3);
                            topology_dihedral->SetDihedrals(dihedral_atom_names3);
                        }
                        if(permutation_index == 1 || (permutation_index >= 10 && permutation_index <= 13) || (permutation_index >= 36 && permutation_index <= 41) || (permutation_index >= 70 && permutation_index <= 73))
                        {
                            topology_dihedral->SetResidueNames(reverse_residue_names1);
                            topology_dihedral->SetDihedrals(reverse_dihedral_atom_names1);
                        }
                        if(permutation_index == 3 || (permutation_index >= 18 && permutation_index <= 21) || (permutation_index >= 48 && permutation_index <= 53) || (permutation_index >= 78 && permutation_index <= 81))
                        {
                            topology_dihedral->SetResidueNames(reverse_residue_names2);
                            topology_dihedral->SetDihedrals(reverse_dihedral_atom_names2);
                        }
                        if(permutation_index == 5 || (permutation_index >= 26 && permutation_index <= 29) || (permutation_index >= 60 && permutation_index <= 65) || (permutation_index >= 86 && permutation_index <= 89))
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
                    if(permutation_index == 0 || (permutation_index >= 6 && permutation_index <= 9) || (permutation_index >= 30 && permutation_index <= 35) || (permutation_index >= 66 && permutation_index <= 69))
                    {
                        inserted_dihedrals.push_back(dihedral1);
                    }
                    if(permutation_index == 2 || (permutation_index >= 14 && permutation_index <= 17) || (permutation_index >= 42 && permutation_index <= 47) || (permutation_index >= 74 && permutation_index <= 77))
                    {
                        inserted_dihedrals.push_back(dihedral2);
                    }
                    if(permutation_index == 4 || (permutation_index >= 22 && permutation_index <= 25) || (permutation_index >= 54 && permutation_index <= 59) || (permutation_index >= 82 && permutation_index <= 85))
                    {
                        inserted_dihedrals.push_back(dihedral3);
                    }
                    if(permutation_index == 1 || (permutation_index >= 10 && permutation_index <= 13) || (permutation_index >= 36 && permutation_index <= 41) || (permutation_index >= 70 && permutation_index <= 73))
                    {
                        inserted_dihedrals.push_back(reverse_dihedral1);
                    }
                    if(permutation_index == 3 || (permutation_index >= 18 && permutation_index <= 21) || (permutation_index >= 48 && permutation_index <= 53) || (permutation_index >= 78 && permutation_index <= 81))
                    {
                        inserted_dihedrals.push_back(reverse_dihedral2);
                    }
                    if(permutation_index == 5 || (permutation_index >= 26 && permutation_index <= 29) || (permutation_index >= 60 && permutation_index <= 65) || (permutation_index >= 86 && permutation_index <= 89))
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
    LibraryFile* library_file = new LibraryFile();
    LibraryFile::ResidueMap residue_map = LibraryFile::ResidueMap();
    ResidueVector residues_of_assembly = this->GetAllResiduesOfAssembly();
    for(ResidueVector::iterator it = residues_of_assembly.begin(); it != residues_of_assembly.end(); it++)
    {
        Residue* assembly_residue = *it;
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
            AtomVector atom_neighbours = atom_node->GetNodeNeighbors();
            for(AtomVector::iterator it2 = atom_neighbours.begin(); it2 != atom_neighbours.end(); it2++)
            {
                int bonded_atom_index = distance(assembly_residue_atoms.begin(), it2);
                bonded_atom_indices.push_back(bonded_atom_index);
            }
            LibraryFileAtom* atom = new LibraryFileAtom(residue_atom->GetAtomType(), residue_atom->GetName(), residue_index, atom_index,
                                                        gmml::iNotSet, residue_atom->GetCharge(),
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
    model_index_ = model_index;
    AtomVector all_atoms_of_assembly = this->GetAllAtomsOfAssembly();
    int i = 0;
    for(AtomVector::iterator it = all_atoms_of_assembly.begin(); it != all_atoms_of_assembly.end(); it++)
    {
        Atom* atom = (*it);
        AtomNode* atom_node = new AtomNode();
        atom_node->SetAtom(atom);
        atom_node->SetId(i);
        i++;
        for(AtomVector::iterator it1 = all_atoms_of_assembly.begin(); it1 != all_atoms_of_assembly.end(); it1++)
        {
            if(it != it1)
            {
                Atom* neighbor_atom = (*it1);
                if((atom->GetCoordinates().at(model_index)->Distance(*(neighbor_atom->GetCoordinates().at(model_index)))) < cutoff)
                {
                    atom_node->AddNodeNeighbor(neighbor_atom);
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
        cout << "Building structure ..." << endl;
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
    cout << "Building structure ..." << endl;
    TopologyFile* topology_file = new TopologyFile(this->GetSourceFile());
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
                ss << atom_1->GetId() << "-" << atom_2->GetId();
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
    cout << "Building structure ..." << endl;
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
                        if(assembly_atom_id.compare(library_atom_id) == 0)
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
    cout << "Building structure ..." << endl;
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
                        if(assembly_atom->GetId().compare(ss.str()) == 0)
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
    cout << "Building structure ..." << endl;
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

void Assembly::CalculateCenterOfGeometry()
{
    Coordinate center_of_geometry = Coordinate();
    int counter = 0;
    for(AssemblyVector::iterator it = assemblies_.begin(); it != assemblies_.end(); it++)
    {
        Assembly* assembly = (*it);
        assembly->CalculateCenterOfGeometry();
        counter++;
        center_of_geometry.operator +(assembly->GetCenterOfGeometry());
    }
    for(ResidueVector:: iterator it1 = residues_.begin(); it1 != residues_.end(); it1++)
    {
        Residue* residue = (*it1);
        Residue::AtomVector atoms = residue->GetAtoms();
        for(Residue::AtomVector::iterator it2 = atoms.begin(); it2 != atoms.end(); it2++)
        {
            Atom* atom = (*it2);
            Atom::CoordinateVector coordinates = atom->GetCoordinates();
            Geometry::Coordinate coordinate = *coordinates[this->model_index_];
            center_of_geometry.operator +(coordinate);
            counter++;
        }
    }
    center_of_geometry.operator /(counter);
    center_of_geometry_ = Coordinate(center_of_geometry);
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
        AtomVector node_neighbors = atom_node->GetNodeNeighbors();
        counter += node_neighbors.size();
    }
    return counter/2;
}

int Assembly::CountNumberOfBondTypes()
{
    AtomVector atoms = GetAllAtomsOfAssembly();
    vector<string> type_list = vector<string>();
    for(AtomVector::iterator it = atoms.begin(); it != atoms.end(); it++)
    {
        Atom* atom = (*it);
        string atom_bond_type = atom->GetAtomType();
        AtomNode* atom_node = atom->GetNode();
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
                type_list.push_back(key);
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
    return counter/2;
}

int Assembly::CountNumberOfAngleTypes()
{
    AtomVector atoms = GetAllAtomsOfAssembly();
    vector<string> type_list = vector<string>();
    for(AtomVector::iterator it = atoms.begin(); it != atoms.end(); it++)
    {
        Atom* atom = (*it);
        string type1 = atom->GetAtomType();
        stringstream ss;
        ss << atom->GetId();
        AtomNode* atom_node = atom->GetNode();
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
                        type_list.push_back(ss2.str());
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

void Assembly::CycleDetection()
{
    vector<vector<string> > cycles = vector<vector<string> >();
    vector<string> cycle = vector<string>();
    //    vector<string> visited_atoms = vector<string>();
    vector<vector<string> > cycle_permutations = vector<vector<string> >();
    ResidueVector residues = GetAllResiduesOfAssembly();
    for(ResidueVector:: iterator res = residues.begin(); res != residues.end(); res++)
    {
        Residue* residue = (*res);
        if(residue->GetName().compare("HOH") == 0)
            continue;
        cout << "residue: " << residue->GetId() << endl;
        AtomVector atoms = residue->GetAtoms();
        for(AtomVector:: iterator it1 = atoms.begin(); it1 != atoms.end(); it1++)
        {
            cycle.clear();
            Atom* atom1 = (*it1);
            string id1 = atom1->GetId();
            cycle.push_back(id1);
            AtomNode* node1 = atom1->GetNode();
            AtomVector atom1_neigbors = node1->GetNodeNeighbors();
            for(AtomVector:: iterator it2 = atom1_neigbors.begin(); it2 != atom1_neigbors.end(); it2++)
            {
                Atom* atom2 = (*it2);
                string id2 = atom2->GetId();
                cycle.push_back(id2);
                AtomNode* node2 = atom2->GetNode();
                AtomVector atom2_neigbors = node2->GetNodeNeighbors();
                for(AtomVector:: iterator it3 = atom2_neigbors.begin(); it3 != atom2_neigbors.end(); it3++)
                {
                    Atom* atom3 = (*it3);
                    string id3 = atom3->GetId();
                    if(id1.compare(id3) != 0)
                    {
                        cycle.push_back(id3);
                        AtomNode* node3 = atom3->GetNode();
                        AtomVector atom3_neigbors = node3->GetNodeNeighbors();
                        for(AtomVector:: iterator it4 = atom3_neigbors.begin(); it4 != atom3_neigbors.end(); it4++)
                        {
                            Atom* atom4 = (*it4);
                            string id4 = atom4->GetId();
                            if(id2.compare(id4) != 0)
                            {
                                cycle.push_back(id4);
                                AtomNode* node4 = atom4->GetNode();
                                AtomVector atom4_neigbors = node4->GetNodeNeighbors();
                                for(AtomVector:: iterator it5 = atom4_neigbors.begin(); it5 != atom4_neigbors.end(); it5++)
                                {
                                    Atom* atom5 = (*it5);
                                    string id5 = atom5->GetId();
                                    if(id3.compare(id5) != 0)
                                    {
                                        cycle.push_back(id5);
                                        AtomNode* node5 = atom5->GetNode();
                                        AtomVector atom5_neigbors = node5->GetNodeNeighbors();
                                        for(AtomVector:: iterator it6 = atom5_neigbors.begin(); it6 != atom5_neigbors.end(); it6++)
                                        {
                                            Atom* atom6 = (*it6);
                                            string id6 = atom6->GetId();
                                            if(id4.compare(id6) != 0)
                                            {
                                                AtomNode* node6 = atom6->GetNode();
                                                AtomVector atom6_neigbors = node6->GetNodeNeighbors();
                                                for(AtomVector:: iterator it7 = atom6_neigbors.begin(); it7 != atom6_neigbors.end(); it7++)
                                                {
                                                    Atom* atom7 = (*it7);
                                                    string id7 = atom7->GetId();
                                                    if(id1.compare(id7) == 0)
                                                    {
                                                        cout << id1 << ";" << id7 << endl;
                                                        cycle.push_back(id6);
                                                        cycle_permutations = CreateAllCyclePermutations(id1, id2, id3, id4, id5, id6);
                                                        bool is_cycle_existed = false;
                                                        for(vector<vector<string> >:: iterator perms = cycle_permutations.begin(); perms != cycle_permutations.end(); perms++)
                                                        {
                                                            if(is_cycle_existed)
                                                                break;
                                                            vector<string> perm = (*perms);
                                                            stringstream ss_perm;
                                                            ss_perm << perm.at(0) << "; " << perm.at(1) << "; " << perm.at(2) << "; " <<
                                                                       perm.at(3) << "; " << perm.at(4) << "; " << perm.at(5);
                                                            for(vector<vector<string> >:: iterator c = cycles.begin(); c != cycles.end(); c++)
                                                            {
                                                                vector<string> existing_cycle = (*c);
                                                                stringstream ss_cycle;
                                                                ss_cycle << existing_cycle.at(0) << "; " << existing_cycle.at(1) << "; " << existing_cycle.at(2) << "; " <<
                                                                            existing_cycle.at(3) << "; " << existing_cycle.at(4) << "; " << existing_cycle.at(5);
                                                                if(ss_perm.str().compare(ss_cycle.str()) == 0)
                                                                {
                                                                    is_cycle_existed = true;
                                                                    break;
                                                                }
                                                            }
                                                        }
                                                        if(!is_cycle_existed)
                                                        {
                                                            cycles.push_back(cycle);

                                                            cout << cycle.at(0) << "; " << cycle.at(1) << "; " << cycle.at(2) << "; " <<
                                                                    cycle.at(3) << "; " << cycle.at(4) << "; " << cycle.at(5) << endl;
                                                            break;
                                                        }
                                                    }
                                                    //                                                    cycle.pop_back();
                                                }
                                            }
                                        }
                                    }
                                    if(!cycle.empty())
                                        cycle.pop_back();
                                }
                            }
                            if(!cycle.empty())
                                cycle.pop_back();
                        }
                    }
                    if(!cycle.empty())
                        cycle.pop_back();
                }
                if(!cycle.empty())
                    cycle.pop_back();
            }
        }
    }
}
std::vector<std::vector<std::string> > Assembly::CreateAllCyclePermutations(string id1, string id2, string id3, string id4, string id5, string id6)
{
    vector<vector<string> > all_permutations = vector<vector<string> >();
    vector<string> order = vector<string>();
    vector<string> rev_order = vector<string>();
    order.push_back(id1);
    order.push_back(id2);
    order.push_back(id3);
    order.push_back(id4);
    order.push_back(id5);
    order.push_back(id6);
    rev_order.push_back(id6);
    rev_order.push_back(id5);
    rev_order.push_back(id4);
    rev_order.push_back(id3);
    rev_order.push_back(id2);
    rev_order.push_back(id1);
    all_permutations.push_back(order);
    all_permutations.push_back(rev_order);
    order.clear();
    rev_order.clear();
    order.push_back(id2);
    order.push_back(id3);
    order.push_back(id4);
    order.push_back(id5);
    order.push_back(id6);
    order.push_back(id1);
    rev_order.push_back(id1);
    rev_order.push_back(id6);
    rev_order.push_back(id5);
    rev_order.push_back(id4);
    rev_order.push_back(id3);
    rev_order.push_back(id2);
    all_permutations.push_back(order);
    all_permutations.push_back(rev_order);
    order.clear();
    rev_order.clear();
    order.push_back(id3);
    order.push_back(id4);
    order.push_back(id5);
    order.push_back(id6);
    order.push_back(id1);
    order.push_back(id2);
    rev_order.push_back(id2);
    rev_order.push_back(id1);
    rev_order.push_back(id6);
    rev_order.push_back(id5);
    rev_order.push_back(id4);
    rev_order.push_back(id3);
    all_permutations.push_back(order);
    all_permutations.push_back(rev_order);
    order.clear();
    rev_order.clear();
    order.push_back(id4);
    order.push_back(id5);
    order.push_back(id6);
    order.push_back(id1);
    order.push_back(id2);
    order.push_back(id3);
    rev_order.push_back(id3);
    rev_order.push_back(id2);
    rev_order.push_back(id1);
    rev_order.push_back(id6);
    rev_order.push_back(id5);
    rev_order.push_back(id4);
    all_permutations.push_back(order);
    all_permutations.push_back(rev_order);
    order.clear();
    rev_order.clear();
    order.push_back(id5);
    order.push_back(id6);
    order.push_back(id1);
    order.push_back(id2);
    order.push_back(id3);
    order.push_back(id4);
    rev_order.push_back(id4);
    rev_order.push_back(id3);
    rev_order.push_back(id2);
    rev_order.push_back(id1);
    rev_order.push_back(id6);
    rev_order.push_back(id5);
    all_permutations.push_back(order);
    all_permutations.push_back(rev_order);
    order.clear();
    rev_order.clear();
    order.push_back(id6);
    order.push_back(id1);
    order.push_back(id2);
    order.push_back(id3);
    order.push_back(id4);
    order.push_back(id5);
    rev_order.push_back(id5);
    rev_order.push_back(id4);
    rev_order.push_back(id3);
    rev_order.push_back(id2);
    rev_order.push_back(id1);
    rev_order.push_back(id6);
    all_permutations.push_back(order);
    all_permutations.push_back(rev_order);
    return all_permutations;
}

Assembly::CycleMap Assembly::DetectCyclesByDFS()
{
    int counter = 0;
//    vector<string> size_vector = Split(cycle_size, "|");

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
    cout << "Number of cycles found: " << counter << endl;
    for(AtomIdAtomMap::iterator it = src_dest_map.begin(); it != src_dest_map.end(); it++)
    {
        string source_id = (*it).first;
        Atom* destination = (*it).second;
        cycle.clear();
        stringstream cycle_stream;
        ReturnCycleAtoms(source_id, destination, atom_parent_map, cycle, cycle_stream);
//        if(find(size_vector.begin(), size_vector.end(), ConvertT(cycle.size())) != size_vector.end())
        cycles[cycle_stream.str()] = cycle;
    }
    for(CycleMap::iterator it = cycles.begin(); it != cycles.end(); it++)
    {
        string cycle_atoms_str = (*it).first;
        AtomVector cycle_atoms = (*it).second;
        cout << cycle_atoms_str << endl;
//        for(AtomVector::iterator it1 = cycle_atoms.begin(); it1 != cycle_atoms.end(); it1++)
//        {
//            Atom* cycle_atom = (*it1);
//            cout << cycle_atom->GetId() << ";";
//        }
//        cout << endl;
        Atom* anomeric = FindAnomericCarbon(cycle_atoms);
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
                    src_dest_map[neighbor->GetId()] = atom;
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

void Assembly::RemoveFusedCycles(CycleMap &cycles)
{
//    for(CycleMap::iterator it = cycles.begin(); it != cycles.end(); it++)
//    {
//        AtomVector cycle_i = (*it).second;
//        for(AtomVector::iterator it1 = cycle_i.begin(); it1 != cycle_i.end(); it1++)
//        {
//            Atom* atom_i = (*it1);
//        }

//        for(CycleMap::iterator it1 = cycles.begin(); it1 != cycles.end(); it1++)
//        {
//            if(it != it1)///comparing cycle i with all the other cycles except itself
//            {
//                AtomVector cycle_j = (*it1).second;
//                for(AtomVector::iterator it2 = cycle_j.begin(); it2 != cycle_j.end(); it2++)
//                {
//                    Atom* atom_j = (*it2);
//                }

//                stringstream atom_pair;
//                stringstream atom_pair_rev;
//                for(int i = 0; i < cycle_i.size(); i++)
//                {
//                    if(i == cycle_i.size() - 1)
//                    {
//                        Atom* a1 = cycle_i.at(i);
//                        Atom* a2 = cycle_i.at(0);
//                    }
//                    else
//                    {
//                        Atom* a1 = cycle_i.at(i);
//                        Atom* a2 = cycle_i.at(i + 1);
//                    }
//                    atom_pair << a1->GetId() << "-" << a2->GetId();
//                    atom_pair_rev << a2->GetId() << "-" << a1->GetId();

//                    if(cycle_i_str.find(atom_pair.str()) != string::npos)
//                    {
//                        to_be_deleted_cycles.push_back(cycle_i);
//                        to_be_deleted_cycles.push_back(cycle_j);
//                        break;
//                    }
//                    else if(cycle_i_str.find(atom_pair_rev.str()) != string::npos)
//                    {
//                        to_be_deleted_cycles.push_back(cycle_i);
//                        to_be_deleted_cycles.push_back(cycle_j);
//                        break;
//                    }

//                }
//            }
//            }
//        }
//    }
}

Atom* Assembly::FindAnomericCarbon(AtomVector cycle)
{
    Atom* anomeric_carbon = new Atom();
    for(AtomVector::iterator it = cycle.begin(); it != cycle.end(); it++)
    {
        Atom* cycle_atom = (*it);
        if((cycle_atom->GetName().substr(0,1).compare("O") == 0 && isdigit(ConvertString<char>(cycle_atom->GetName().substr(1,1)))))///find oxygen in ring
        {
            AtomNode* node = cycle_atom->GetNode();
            AtomVector neighbors = node->GetNodeNeighbors();

            Atom* o_neighbor1 = neighbors.at(0);
            AtomNode* o_neighbor1_node = o_neighbor1->GetNode();
            AtomVector o_neighbor1_neighbors = o_neighbor1_node->GetNodeNeighbors();
            for(AtomVector::iterator it1 = o_neighbor1_neighbors.begin(); it1 != o_neighbor1_neighbors.end(); it1++)///check if neighbor1 of oxygen has another oxygen neighbor
            {
                Atom* neighbor1_neighbor = (*it1);
                if((neighbor1_neighbor->GetName().substr(0,1).compare("O") == 0 && isdigit(ConvertString<char>(neighbor1_neighbor->GetName().substr(1,1)))))
                {
                    anomeric_carbon = o_neighbor1;
                    cout << "anomeric carbon is: " << anomeric_carbon->GetName() << endl;
                    return anomeric_carbon;
                }
            }

            Atom* o_neighbor2 = neighbors.at(1);
            AtomNode* o_neighbor2_node = o_neighbor2->GetNode();
            AtomVector o_neighbor2_neighbors = o_neighbor2_node->GetNodeNeighbors();
            for(AtomVector::iterator it2 = o_neighbor2_neighbors.begin(); it2 != o_neighbor2_neighbors.end(); it2++)///check if neighbor2 of oxygen has another oxygen neighbor
            {
                Atom* neighbor2_neighbor = (*it2);
                if((neighbor2_neighbor->GetName().substr(0,1).compare("O") == 0 && isdigit(ConvertString<char>(neighbor2_neighbor->GetName().substr(1,1)))))
                {
                    anomeric_carbon = o_neighbor2;
                    cout << "anomeric carbon is: " << anomeric_carbon->GetName() << endl;
                    return anomeric_carbon;
                }
            }
        }
    }
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
