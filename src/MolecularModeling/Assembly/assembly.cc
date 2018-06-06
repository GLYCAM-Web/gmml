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
#include "../../../includes/MolecularModeling/residuenode.hpp"         //Added by ayush on 11/16/17 for identifying residuenodes in assembly
#include "../../../includes/MolecularModeling/Molecule.hpp"            //Added by ayush on 11/12/17 for molecules in assembly

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

Assembly::Assembly(string file_path, gmml::InputFileType type)
{
    source_file_type_ = type;
    description_ = "";
    model_index_ = 0;
    sequence_number_ = 1;
    source_file_ = file_path;
    residues_ = ResidueVector();
    assemblies_ = AssemblyVector();
    id_ = "1";
    switch(type)
    {
    case gmml::PDB:
        BuildAssemblyFromPdbFile(source_file_);
        break;
    case gmml::PDBQT:
        BuildAssemblyFromPdbqtFile(source_file_);
        break;
    case gmml::TOP:
        BuildAssemblyFromTopologyFile(source_file_);
        break;
    case gmml::LIB:
        BuildAssemblyFromLibraryFile(source_file_);
        break;
    case gmml::PREP:
        BuildAssemblyFromPrepFile(source_file_);
        break;
    case gmml::MULTIPLE:
        break;
    case gmml::UNKNOWN:
        break;
    default:
        std::cout << "Error, problem with input file type in Assembly Constructor" << std::endl;
    }
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
    default:
        std::cout << "Error, input type not recognized in Assembly Constructor" << std::endl;
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

Assembly::AtomVector Assembly::GetAllAtomsOfAssemblyWithinProteinBackbone()
{ // written by OG, so syntax a bit different from other functions.

    AtomVector selection_from_assembly = AtomVector();
    AtomVector all_protein_atoms = this->GetAllAtomsOfAssemblyWithinProteinResidues();
    for(AtomVector::iterator it = all_protein_atoms.begin(); it != all_protein_atoms.end(); it++)
    {
        Atom *atom = *it;
        if ( (atom->GetName().compare("N")==0) || (atom->GetName().compare("C")==0) || (atom->GetName().compare("O")==0) || (atom->GetName().compare("CA")==0) )
        {
            selection_from_assembly.push_back(atom);
        }
    }
    return selection_from_assembly;
}

Assembly::AtomVector Assembly::GetAllAtomsOfAssemblyWithinProteinSidechain()
{ // written by OG, so syntax a bit different from other functions.

    AtomVector selection_from_assembly = AtomVector();
    AtomVector all_protein_atoms = this->GetAllAtomsOfAssemblyWithinProteinResidues();
    for(AtomVector::iterator it = all_protein_atoms.begin(); it != all_protein_atoms.end(); it++)
    {
        Atom *atom = *it;
        if ( (atom->GetName().compare("N")!=0) && (atom->GetName().compare("C")!=0) && (atom->GetName().compare("O")!=0) && (atom->GetName().compare("CA")!=0) )
        {
            selection_from_assembly.push_back(atom);
        }
    }
    return selection_from_assembly;
}

Assembly::AtomVector Assembly::GetAllAtomsOfAssemblyWithinProteinResidues()
{ // written by OG, so syntax a bit different from other functions.
    AtomVector selection_from_assembly = AtomVector();
    AssemblyVector assemblies = this->GetAssemblies();
    //std::cout << "size is " << assemblies.size() << std::endl;
    for(AssemblyVector::iterator it = assemblies.begin(); it != assemblies.end(); it++)
    {
        Assembly* assembly = (*it);
        ResidueVector residues = assembly->GetResidues();
        for(ResidueVector::iterator it1= residues.begin(); it1 != residues.end(); it1++)
        {
            Residue* residue = (*it1);
            if (residue->CheckIfProtein())
            {
                //std::cout << "Found a protein" << std::endl;
                AtomVector atoms = residue->GetAtoms();
                for(AtomVector::iterator it2 = atoms.begin(); it2 != atoms.end(); it2++)
                {
                    Atom *atom = *it2;
                    selection_from_assembly.push_back(atom);
                }

            }

        }
    }
    // This is unintuitive, but GetAssemblies does not return "this" assembly, just additonal "sub-assemblies" contained within this assembly. Horrific.
    ResidueVector residues =  this->GetResidues();
    for(ResidueVector::iterator it1= residues.begin(); it1 != residues.end(); it1++)
    {
        Residue* residue = (*it1);
        //std::cout << "Checking for a protein" << std::endl;
        if (residue->CheckIfProtein())
        {
            //std::cout << "Found a protein" << std::endl;
            AtomVector atoms = residue->GetAtoms();
            for(AtomVector::iterator it2 = atoms.begin(); it2 != atoms.end(); it2++)
            {
                Atom *atom = *it2;
                selection_from_assembly.push_back(atom);
            }

        }

    }
    return selection_from_assembly;
}

Assembly::AtomVector Assembly::GetAllAtomsOfAssemblyNotWithinProteinResidues()
{ // written by OG, so syntax a bit different from other functions.
    AtomVector selection_from_assembly = AtomVector();
    AssemblyVector assemblies = this->GetAssemblies();
    //std::cout << "size is " << assemblies.size() << std::endl;
    for(AssemblyVector::iterator it = assemblies.begin(); it != assemblies.end(); it++)
    {
        Assembly* assembly = (*it);
        ResidueVector residues = assembly->GetResidues();
        for(ResidueVector::iterator it1= residues.begin(); it1 != residues.end(); it1++)
        {
            Residue* residue = (*it1);
            if (residue->CheckIfProtein()) {}
            else
            {
                //std::cout << "Found a non-protein" << std::endl;
                AtomVector atoms = residue->GetAtoms();
                for(AtomVector::iterator it2 = atoms.begin(); it2 != atoms.end(); it2++)
                {
                    Atom *atom = *it2;
                    selection_from_assembly.push_back(atom);
                }

            }

        }
    }
    // This is unintuitive, but GetAssemblies does not return "this" assembly, just additonal "sub-assemblies" contained within this assembly. Horrific.
    ResidueVector residues =  this->GetResidues();
    for(ResidueVector::iterator it1= residues.begin(); it1 != residues.end(); it1++)
    {
        Residue* residue = (*it1);
        //std::cout << "Checking for a protein" << std::endl;
        if (residue->CheckIfProtein()) {}
        else
        {
            //std::cout << "Found a non-protein" << std::endl;
            AtomVector atoms = residue->GetAtoms();
            for(AtomVector::iterator it2 = atoms.begin(); it2 != atoms.end(); it2++)
            {
                Atom *atom = *it2;
                selection_from_assembly.push_back(atom);
            }

        }

    }
    return selection_from_assembly;
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

Assembly::CoordinateVector Assembly::GetCycleAtomCoordinates( Monosaccharide* mono ) {
  CoordinateVector coordinates;
  for( AtomVector::iterator it1 = mono->cycle_atoms_.begin(); it1 != mono->cycle_atoms_.end(); it1++ ) {
    Atom* atom = ( *it1 );
    CoordinateVector atom_coordinates = atom->GetCoordinates();
    for( CoordinateVector::iterator it2 = atom_coordinates.begin(); it2 != atom_coordinates.end(); it2++ ) {
      coordinates.push_back( ( *it2 ) );
    }
  }
  return coordinates;
}

Assembly::NoteVector Assembly::GetNotes()
{
    return notes_;
}

  //Added by ayush on 11/16/17 for residuenodes in assembly
Assembly::ResidueNodeVector Assembly::GetAllResidueNodesOfAssembly()
{
        return residuenodes_;
}

 //Added by ayush on 11/12/17 for molecules in assembly
Assembly::MoleculeVector Assembly::GetMolecules()
{
      return molecules_;
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


//Added by ayush on 11/16/17 for residuenodes in assembly
void Assembly::SetResidueNodes(ResidueNodeVector residuenodes)
{
    residuenodes_.clear();
    for(ResidueNodeVector::iterator it = residuenodes.begin(); it != residuenodes.end(); it++)
    {
        residuenodes_.push_back(*it);
    }
}


//Added by ayush on 11/12/17 for molecules in assembly
void Assembly::SetMolecules(MoleculeVector molecules)
{
    molecules_.clear();
    for(MoleculeVector::iterator it = molecules.begin(); it != molecules.end(); it++)
    {
        molecules_.push_back(*it);
	}
}
void Assembly::MergeAssembly(Assembly *other) // Added by Oliver. He is unsure and this may well cause problems.
{
    ResidueVector residues = other->GetResidues();
    for (ResidueVector::iterator it = residues.begin(); it != residues.end(); ++it)
    {
        Residue *residue = *it;
        this->AddResidue(residue);
    }
}

//////////////////////////////////////////////////////////
//                       FUNCTIONS                      //
//////////////////////////////////////////////////////////
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


//Added by ayush on 11/16/17 for residuenodes in assembly
void Assembly::AddResidueNode(ResidueNode* residuenode)
{
    residuenodes_.push_back(residuenode);
}

Assembly::ResidueNodeVector Assembly::GenerateResidueNodesInAssembly()
{

    int residue_counter=1;
    ResidueVector assembly_residues=this->GetAllResiduesOfAssembly();
    for(ResidueVector::iterator it = assembly_residues.begin(); it != assembly_residues.end(); it++)
    {
        ResidueNode* residuenode = new ResidueNode();

       Residue* residue = (*it);

        residuenode->SetResidue(residue);

        residuenode->SetId(residue_counter);

        string residue_id=residue->GetId();

        AtomVector residue_atoms = residue->GetAtoms();

        for(AtomVector::iterator it1 = residue_atoms.begin(); it1 != residue_atoms.end(); it1++)
        {
                  Atom* residue_atom = *it1;
                  AtomVector atom_neighbors = residue_atom->GetNode()->GetNodeNeighbors();

                    for(AtomVector::iterator it2 =atom_neighbors.begin(); it2 != atom_neighbors.end(); it2++)
                    {
                             Atom* neighbor_atom = *it2;

                              if((neighbor_atom->GetResidue()->GetId()).compare(residue_id)!=0)
                              {
                                  residuenode->AddResidueNodeConnectingAtom(residue_atom);

                              }
                    }

        }
            residue->SetNode(residuenode);
            this->AddResidueNode(residuenode);
            residue_counter++;

    }

    //For adding residuenode neighbours
    for(ResidueNodeVector::iterator it3 = residuenodes_.begin(); it3 != residuenodes_.end(); it3++)
    {
            ResidueNode* residuenode = (*it3);

             AtomVector residue_connecting_atoms = residuenode->GetResidueNodeConnectingAtoms();

            for(AtomVector::iterator it4 = residue_connecting_atoms.begin(); it4 != residue_connecting_atoms.end(); it4++)
            {
                Atom* connecting_atom = (*it4);


                AtomVector connecting_atom_neighbors= connecting_atom->GetNode()->GetNodeNeighbors();
                for(AtomVector::iterator it5 = connecting_atom_neighbors.begin(); it5 != connecting_atom_neighbors.end(); it5++)
                {

                        Atom* connecting_atom_neighbor = (*it5);


                        if((connecting_atom_neighbor->GetResidue()->GetId()).compare(residuenode->GetResidue()->GetId())!=0)
                        {
                                ResidueNode* neighbor_residue=connecting_atom_neighbor->GetResidue()->GetNode();
                                 residuenode->AddResidueNodeNeighbor(neighbor_residue);
                        }

                }

            }
    }
        return residuenodes_;
  }

//Added by ayush on 12/8/17 for molecules in assembly
void Assembly::GenerateMoleculesInAssembly()
{
    int moleculecount=0;

    for(ResidueNodeVector::iterator it = residuenodes_.begin(); it != residuenodes_.end(); it++)
    {
         ResidueNode* residuenode = (*it);
         if(residuenode->GetIsVisited()==false)
         {
             GenerateMoleculesDFSUtil(residuenode);
             moleculecount++;
         }
    }

    cout<<"Total number of molecules in Assembly is:"<<moleculecount<<endl;

}

void Assembly::GenerateMoleculesDFSUtil(ResidueNode* DFSresiduenode)
{
        DFSresiduenode->SetIsVisited(true);
        ResidueNodeVector residuenode_neighbors=DFSresiduenode->GetResidueNodeNeighbors();

        for(ResidueNodeVector::iterator it2 = residuenode_neighbors.begin(); it2 != residuenode_neighbors.end(); it2++)
        {
               ResidueNode* current_residuenode_neighbor = (*it2);
               if(current_residuenode_neighbor->GetIsVisited()==false)
               {
                    GenerateMoleculesDFSUtil(current_residuenode_neighbor);
               }

        }


}

//Added by ayush on 11/12/17 for molecules in assembly
void Assembly::AddMolecule(Molecule *molecule)
{
    int max_index_=molecules_.size()+1;
    molecule->SetMoleculeIndex(max_index_);
    molecules_.push_back(molecule);
}


 //Added by ayush on 04/11/18 for TopologyFix in assembly
Assembly::AtomVector Assembly::GetAllBondedAtomsByStartDirection(Atom* start_atom, Atom* direction_atom , AtomVector ignore_list)
{

        if(start_atom == NULL || direction_atom == NULL){
            cout<<"Start Atom or Direction Atom is null"<<endl;}
        else if(CheckIfAtomExistInAssembly(start_atom)==false||CheckIfAtomExistInAssembly(direction_atom)==false){
            cout<<"Start Atom or Direction Atom does not exist in Assembly"<<endl;
        }else{
                bool isNeighbor=false;
                 AtomVector start_atom_neighbors = start_atom->GetNode()->GetNodeNeighbors();

                 for(AtomVector::iterator it = start_atom_neighbors.begin(); it != start_atom_neighbors.end(); it++)
                 {
                     Atom* atom = *it;
                    if(atom->GetId().compare(direction_atom->GetId())==0)
                    { isNeighbor=true;}
                 }
                 if(isNeighbor==true)
                 {
                   
                    start_atom->GetNode()->SetIsVisited(true);                   
                    direction_atom->GetNode()->SetIsVisited(true);
                    bonded_atoms_bystartdirection_.push_back(direction_atom);
                    AtomVector direction_atom_neighbors=direction_atom->GetNode()->GetNodeNeighbors();
                    for(AtomVector::iterator it = direction_atom_neighbors.begin(); it != direction_atom_neighbors.end(); it++)
                    {
                          Atom* neighbor_atom = *it;
                          if(neighbor_atom->GetNode()->GetIsVisited()==false)
                          {
                                BondedAtomsByStartDirectionDFSUtil(neighbor_atom, direction_atom_neighbors, ignore_list);
                          }
                    }

                 }else{
                     cout<<"Direction Atom is invalid. Not a neighbor of Start Atom."<<endl;
                 }
        }

        return bonded_atoms_bystartdirection_;
}

void Assembly::BondedAtomsByStartDirectionDFSUtil(Atom* DFSatom, AtomVector start_atom_neighbors, AtomVector ignore_list)
{
    bool toConsider=false; //Check if the atom is part of the start atom neighor or ignore list and set toConsider accordingly. If toConsider is true, skip the atom.

    DFSatom->GetNode()->SetIsVisited(true);

    for(AtomVector::iterator it1 = start_atom_neighbors.begin(); it1!= start_atom_neighbors.end(); it1++)
    {
        Atom* atom = *it1;
       if(atom->GetId().compare(DFSatom->GetId())==0)
       {
           toConsider = true;
       }
    }

     for(AtomVector::iterator it2 = ignore_list.begin(); it2!= ignore_list.end(); it2++)
    {
        Atom* atom = *it2;
       if(atom->GetId().compare(DFSatom->GetId())==0)
       {
           toConsider = false;
       }
    }

    if(toConsider==false)
    {
        AtomVector DFSatom_neighbors = DFSatom->GetNode()->GetNodeNeighbors();

        for(AtomVector::iterator it3 = DFSatom_neighbors.begin(); it3!= DFSatom_neighbors.end(); it3++)
        {
            Atom* atom = *it3;
                if(atom->GetNode()->GetIsVisited()==false)
                {
                    BondedAtomsByStartDirectionDFSUtil(DFSatom, start_atom_neighbors, ignore_list);
                }
        }

        bonded_atoms_bystartdirection_.push_back(DFSatom);
    }
}


bool Assembly::CheckIfAtomExistInAssembly(Atom* toCheckAtom)
{
    bool status = false;
    AtomVector all_atoms_of_assembly = this->GetAllAtomsOfAssembly();
    for(AtomVector::iterator it = all_atoms_of_assembly.begin(); it != all_atoms_of_assembly.end(); it++)
    {
        Atom* atom = *it;
       if(atom->GetId().compare(toCheckAtom->GetId())==0)
       {
            status=true;
       }
    }
    return status;
}

 //Added by ayush on 04/16/18 for TopologyFix in assembly

Assembly::CoordinateVector Assembly::GetCoordinatesFromAtomVector(AtomVector atomList, int CoordinateIndex)
{
    CoordinateVector coordinatesByIndex= CoordinateVector();
    for( AtomVector::iterator it1 = atomList.begin(); it1 != atomList.end(); it1++ ) {
      Atom* atom = ( *it1 );
      CoordinateVector atom_coordinates = atom->GetCoordinates();
      if(CoordinateIndex < atom_coordinates.size())
      {
        Coordinate* coordinate = atom_coordinates[CoordinateIndex];
        coordinatesByIndex.push_back(coordinate);
      }else{
          cout<<"Atom ID "<<atom->GetId()<<" does not have coordinate at index: "<<CoordinateIndex<<endl;
      }
    }
    return coordinatesByIndex;
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
