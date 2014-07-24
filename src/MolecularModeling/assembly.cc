#include "../../includes/MolecularModeling/assembly.hpp"
#include "../../includes/MolecularModeling/residue.hpp"
#include "../../includes/MolecularModeling/atom.hpp"
#include "../../includes/MolecularModeling/atomnode.hpp"

#include "../../includes/FileSet/TopologyFileSpace/topologyfile.hpp"
#include "../../includes/FileSet/TopologyFileSpace/topologyassembly.hpp"
#include "../../includes/FileSet/TopologyFileSpace/topologyresidue.hpp"
#include "../../includes/FileSet/TopologyFileSpace/topologyatom.hpp"
#include "../../includes/FileSet/CoordinateFileSpace/coordinatefile.hpp"
#include "../../includes/ParameterSet/PrepFileSpace/prepfile.hpp"
#include "../../includes/ParameterSet/PrepFileSpace/prepfileresidue.hpp"
#include "../../includes/ParameterSet/PrepFileSpace/prepfileatom.hpp"
#include "../../includes/FileSet/PdbFileSpace/pdbfile.hpp"
#include "../../includes/FileSet/PdbFileSpace/pdbmodelcard.hpp"
#include "../../includes/FileSet/PdbFileSpace/pdbmodel.hpp"
#include "../../includes/FileSet/PdbFileSpace/pdbmodelresidueset.hpp"
#include "../../includes/FileSet/PdbFileSpace/pdbatomcard.hpp"
#include "../../includes/FileSet/PdbFileSpace/pdbheterogenatomcard.hpp"
#include "../../includes/FileSet/PdbFileSpace/pdbatom.hpp"
#include "../../includes/ParameterSet/LibraryFileSpace/libraryfile.hpp"
#include "../../includes/ParameterSet/LibraryFileSpace/libraryfileatom.hpp"
#include "../../includes/ParameterSet/LibraryFileSpace/libraryfileresidue.hpp"
#include "../../includes/utils.hpp"

using namespace std;
using namespace MolecularModeling;
using namespace TopologyFileSpace;
using namespace CoordinateFileSpace;
using namespace PrepFileSpace;
using namespace PdbFileSpace;
using namespace Geometry;
using namespace LibraryFileSpace;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
Assembly::Assembly() : description_("") {}

Assembly::Assembly(vector<string> file_paths, gmml::InputFileType type)
{
    source_file_type_ = type;
    description_ = "";
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
    }
}

Assembly::Assembly(vector<vector<string> > file_paths, vector<gmml::InputFileType> types)
{
    stringstream name;
    stringstream source_file;
    for(int i = 0; i < file_paths.size(); i++)
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
    assemblies_.push_back(assembly);
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
//////////////////////////////////////////////////////////
//                       FUNCTIONS                      //
//////////////////////////////////////////////////////////
void Assembly::BuildAssemblyFromPdbFile(string pdb_file_path)
{
    PdbFile* pdb_file = new PdbFile(pdb_file_path);
    PdbFile::PdbResidueAtomsMap residue_atoms_map = pdb_file->GetAllAtomsOfResidues();

    for(PdbFile::PdbResidueAtomsMap::iterator it = residue_atoms_map.begin(); it != residue_atoms_map.end(); it++)
    {
        PdbFile::PdbAtomVector* atoms = (*it).second;
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

            Atom* new_atom = new Atom();
            residue->SetName(residue_name);
            new_atom->SetName(atom->GetAtomName());
            new_atom->SetResidue(residue);
            PdbModelCard* models = pdb_file->GetModels();
            PdbModelCard::PdbModelMap model_maps = models->GetModels();
            if(model_maps.size() == 1)
            {
                new_atom->AddCoordinate(new Geometry::Coordinate(atom->GetAtomOrthogonalCoordinate()));
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
                        }
                    }
                }
            }
            residue->AddAtom(new_atom);
        }
        this->AddResidue(residue);
    }

}

void Assembly::BuildAssemblyFromTopologyFile(string topology_file_path)
{
    TopologyFile* topology_file = new TopologyFile(topology_file_path);
    name_ = topology_file->GetTitle();
    sequence_number_ = 1;
    TopologyAssembly::TopologyResidueMap topology_residues = topology_file->GetAssembly()->GetResidues();
    for(TopologyAssembly::TopologyResidueMap::iterator it = topology_residues.begin(); it != topology_residues.end(); it++)
    {
        Residue* assembly_residue = new Residue();
        assembly_residue->SetAssembly(this);
        assembly_residue->SetName((*it).first);
        TopologyResidue* topology_residue = (*it).second;

        TopologyResidue::TopologyAtomMap topology_atoms = topology_residue->GetAtoms();
        for(TopologyResidue::TopologyAtomMap::iterator it1 = topology_atoms.begin(); it1 != topology_atoms.end(); it1++)
        {
            Atom* assembly_atom = new Atom();
            assembly_atom->SetName((*it1).first);
            TopologyAtom* topology_atom = (*it1).second;

            assembly_atom->SetResidue(assembly_residue);
            assembly_atom->SetName(topology_atom->GetAtomName());

            assembly_residue->AddAtom(assembly_atom);
        }
        residues_.push_back(assembly_residue);

    }
}

void Assembly::BuildAssemblyFromLibraryFile(string library_file_path)
{
    LibraryFile* library_file = new LibraryFile(library_file_path);
    sequence_number_ = 1;
    LibraryFile::ResidueMap library_residues = library_file->GetResidues();
    stringstream ss;

    for(LibraryFile::ResidueMap::iterator it = library_residues.begin(); it != library_residues.end(); it++)
    {
        Residue* assembly_residue = new Residue();
        assembly_residue->SetAssembly(this);
        assembly_residue->SetName((*it).first);
        LibraryFileResidue* library_residue = (*it).second;
        string library_residue_name = library_residue->GetName();
        if(distance(library_residues.begin(), it) == library_residues.size()-1)
            ss << library_residue_name;
        else
            ss << library_residue_name << "-";

        LibraryFileResidue::AtomMap library_atoms = library_residue->GetAtoms();
        for(LibraryFileResidue::AtomMap::iterator it1 = library_atoms.begin(); it1 != library_atoms.end(); it1++)
        {
            Atom* assembly_atom = new Atom();
            LibraryFileAtom* library_atom = (*it1).second;
            assembly_atom->SetName(library_atom->GetName());

            assembly_atom->SetResidue(assembly_residue);
            assembly_atom->SetName(library_atom->GetName());

            Coordinate* coordinate = new Coordinate(library_atom->GetCoordinate());
            assembly_atom->AddCoordinate(coordinate);
            assembly_residue->AddAtom(assembly_atom);
        }
        residues_.push_back(assembly_residue);
    }
    name_ = ss.str();
}

void Assembly::BuildAssemblyFromTopologyCoordinateFile(string topology_file_path, string coordinate_file_path)
{
    TopologyFile* topology_file = new TopologyFile(topology_file_path);
    name_ = topology_file->GetTitle();
    sequence_number_ = 1;
    TopologyAssembly::TopologyResidueMap topology_residues = topology_file->GetAssembly()->GetResidues();
    for(TopologyAssembly::TopologyResidueMap::iterator it = topology_residues.begin(); it != topology_residues.end(); it++)
    {
        Residue* assembly_residue = new Residue();
        assembly_residue->SetAssembly(this);
        assembly_residue->SetName((*it).first);
        TopologyResidue* topology_residue = (*it).second;

        TopologyResidue::TopologyAtomMap topology_atoms = topology_residue->GetAtoms();
        for(TopologyResidue::TopologyAtomMap::iterator it1 = topology_atoms.begin(); it1 != topology_atoms.end(); it1++)
        {
            Atom* assembly_atom = new Atom();
            assembly_atom->SetName((*it1).first);
            TopologyAtom* topology_atom = (*it1).second;

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
    PrepFile* prep_file = new PrepFile(prep_file_path);
    sequence_number_ = 1;
    PrepFile::ResidueMap prep_residues = prep_file->GetResidues();
    stringstream ss;

    for(PrepFile::ResidueMap::iterator it = prep_residues.begin(); it != prep_residues.end(); it++)
    {
        Residue* assembly_residue = new Residue();
        assembly_residue->SetAssembly(this);
        assembly_residue->SetName((*it).first);
        PrepFileResidue* prep_residue = (*it).second;
        string prep_residue_name = prep_residue->GetName();
        if(distance(prep_residues.begin(), it) == prep_residues.size()-1)
            ss << prep_residue_name;
        else
            ss << prep_residue_name << "-";

        PrepFileResidue::PrepFileAtomVector prep_atoms = prep_residue->GetAtoms();
        vector<Coordinate*> coordinates = vector<Coordinate*>();
        for(PrepFileResidue::PrepFileAtomVector::iterator it1 = prep_atoms.begin(); it1 != prep_atoms.end(); it1++)
        {
            Atom* assembly_atom = new Atom();
            PrepFileAtom* prep_atom = (*it1);

            assembly_atom->SetResidue(assembly_residue);
            assembly_atom->SetName(prep_atom->GetName());

            if(prep_residue->GetCoordinateType() == PrepFileSpace::kINT)
            {
                vector<Coordinate*> coordinate_list = vector<Coordinate*>();
                int index = distance(prep_atoms.begin(), it1);
                if(index == 0)
                {

                }
                if(index == 1)
                {
                    coordinate_list.push_back(coordinates.at(index-1));
                }
                if(index == 2)
                {
                    coordinate_list.push_back(coordinates.at(index-2));
                    coordinate_list.push_back(coordinates.at(index-1));
                }
                if(index > 2)
                {
                    coordinate_list.push_back(coordinates.at(index-3));
                    coordinate_list.push_back(coordinates.at(index-2));
                    coordinate_list.push_back(coordinates.at(index-1));
                }
                Coordinate* coordinate = gmml::ConvertInternalCoordinate2CartesianCoordinate(coordinate_list, prep_atom->GetBondLength(),
                                                                                             prep_atom->GetAngle(), prep_atom->GetDihedral());
                coordinates.push_back(coordinate);
                assembly_atom->AddCoordinate(coordinate);
            }
            else if(prep_residue->GetCoordinateType() == PrepFileSpace::kXYZ)
            {
                assembly_atom->AddCoordinate(new Coordinate(prep_atom->GetBondLength(), prep_atom->GetAngle(), prep_atom->GetDihedral()));
            }
            assembly_residue->AddAtom(assembly_atom);
        }
        residues_.push_back(assembly_residue);
    }
    name_ = ss.str();
}

void Assembly::BuildStructure(gmml::BuildingStructureOption building_option, vector<string> options, vector<string> file_paths)
{
    switch(building_option)
    {
        case gmml::DISTANCE:
            if(options.size() == 1)
            {
                vector<string> tokens = gmml::Split(options.at(0), ":");
                if(tokens.at(0).compare("cutoff") == 0)
                {
                    double cutoff = gmml::ConvertString<double>(tokens.at(1));
                    stringstream ss(description_);
                    ss << "Building option: Distance;";
                    description_ = ss.str();
                    BuildStructureByDistance(cutoff);
                }
            }
            else
            {
                BuildStructureByDistance();
            }
            break;
        case gmml::ORIGINAL:
            if(options.size() == 1)
            {
                vector<string> tokens = gmml::Split(options.at(0), ":");
                if(tokens.at(0).compare("type") == 0)
                {
                    gmml::InputFileType type = gmml::ConvertString2AssemblyInputFileType(tokens.at(1));
                    if(file_paths.size() == 1)
                    {
                        stringstream ss(description_);
                        ss << "Building option: Original;";
                        ss << "File type: " << type << ";" << "File path: " << file_paths.at(0) << ";";
                        description_ = ss.str();
                        BuildStructureByOriginalFileBondingInformation(type, file_paths.at(0));
                    }
                }
            }
            break;
        case gmml::DATABASE:
            vector<gmml::InputFileType> types = vector<gmml::InputFileType>();
            for(int i = 0; i < options.size(); i++)
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
                for(int i = 0; i < types.size(); i++)
                {
                    ss << "File type: " << types.at(i) << ";" << "File path: " << file_paths.at(i) << ";";
                }
                description_ = ss.str();
                BuildStructureByDatabaseFilesBondingInformation(types, file_paths);
            }
            break;
    }
}

void Assembly::BuildStructureByOriginalFileBondingInformation(gmml::InputFileType type, string file_path)
{
    switch(type)
    {
    }
}

void Assembly::BuildStructureByDatabaseFilesBondingInformation(vector<gmml::InputFileType> types, vector<string> file_paths)
{

}

Assembly::AtomVector Assembly::GetAllAtomsOfAssembly()
{
    AtomVector all_atoms_of_assembly;
    ResidueVector residues = GetResidues();
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

void Assembly::BuildStructureByDistance(double cutoff, int model_index)
{
    AtomVector all_atoms_of_assembly = GetAllAtomsOfAssembly();
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
            Atom* neighbor_atom = (*it1);
            if(it != it1)
            {
                if((atom->GetCoordinates().at(model_index)->Distance(*(neighbor_atom->GetCoordinates().at(model_index)))) < cutoff)
                    atom_node->AddNodeNeighbor(neighbor_atom);
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
