#include <math.h>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <set>
#include <queue>
#include <stack>
#include <algorithm>

#include "../../../includes/MolecularModeling/assembly.hpp"
#include "../../../includes/MolecularModeling/residue.hpp"
#include "../../../includes/MolecularModeling/residuenode.hpp"
#include "../../../includes/MolecularModeling/atom.hpp"
#include "../../../includes/MolecularModeling/atomnode.hpp"
#include "../../../includes/MolecularModeling/molecule.hpp"
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
#include "../../../includes/ParameterSet/OffFileSpace/offfile.hpp"         //Added by ayush on 06/16/18 for OffFile
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

int local_debug = 0;

using MolecularModeling::Assembly;
using MolecularModeling::ResidueNode;
//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////

Assembly::Assembly() : sequence_number_(1), id_("1"), description_(""), model_index_(0)
{
    residues_ = ResidueVector();
    assemblies_ = AssemblyVector();
}

Assembly::Assembly(std::string file_path, gmml::InputFileType type)
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
//        std::cout << "Error, problem with input file type in Assembly Constructor" << std::endl;
        break;
    }
}

Assembly::Assembly(std::vector<std::string> file_paths, gmml::InputFileType type)
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
//        std::cout << "Error, input type not recognized in Assembly Constructor" << std::endl;
        break;
    }
}

Assembly::Assembly(Assembly *assembly) : sequence_number_(1), id_("1"), description_(""), model_index_(0)
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

Assembly::Assembly(std::vector<std::vector<std::string> > file_paths, std::vector<gmml::InputFileType> types)
{
    std::stringstream name;
    std::stringstream source_file;
    sequence_number_ = 1;
    id_ = "1";
    for(unsigned int i = 0; i < file_paths.size(); i++)
    {
        std::vector<std::string> file = file_paths.at(i);
        gmml::InputFileType input_type = types.at(i);
        Assembly* assembly = new Assembly(file, input_type);
        assembly->SetSequenceNumber(i + 1);
        std::stringstream ss;
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

Assembly::Assembly(std::stringstream& atomStream)
{
  BuildAssemblyFromAtomStream(atomStream);
}

Assembly::Assembly(std::vector<MolecularModeling::Residue*> residueVector)
{
  residues_ = residueVector;
}
//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
std::string Assembly::GetName()
{
    return name_;
}

MolecularModeling::AssemblyVector Assembly::GetAssemblies()
{
    return assemblies_;
}

MolecularModeling::ResidueVector Assembly::GetResidues()
{
    return residues_;
}

std::string Assembly::GetChemicalType()
{
    return chemical_type_;
}

int Assembly::GetSequenceNumber()
{
    return sequence_number_;
}

std::string Assembly::GetId()
{
    return id_;
}

std::string Assembly::GetDescription()
{
    return description_;
}

std::string Assembly::GetSourceFile()
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

MolecularModeling::AtomVector Assembly::GetAllAtomsOfAssembly()
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

MolecularModeling::AtomVector Assembly::GetAllAtomsOfAssemblyWithinProteinBackbone()
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

MolecularModeling::AtomVector Assembly::GetAllAtomsOfAssemblyWithinProteinSidechain()
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

MolecularModeling::AtomVector Assembly::GetAllAtomsOfAssemblyWithinProteinResidues()
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

MolecularModeling::AtomVector Assembly::GetAllAtomsOfAssemblyNotWithinProteinResidues()
{ // written by OG, so syntax a bit different from other functions.
    AtomVector selection_from_assembly = AtomVector();
    MolecularModeling::AssemblyVector assemblies = this->GetAssemblies();
    //std::cout << "size is " << assemblies.size() << std::endl;
    for(MolecularModeling::AssemblyVector::iterator it = assemblies.begin(); it != assemblies.end(); it++)
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


MolecularModeling::AtomVector Assembly::GetAllAtomsOfAssemblyExceptProteinWaterResiduesAtoms()
{
    AtomVector all_atoms_of_assembly = AtomVector();
    AssemblyVector assemblies = this->GetAssemblies();
    if ( local_debug > 0 )
    {
      gmml::log(__LINE__, __FILE__, gmml::INF, "Assemblies size: ");
      gmml::log(__LINE__, __FILE__, gmml::INF, std::to_string(assemblies.size()));
    }
    
    for(AssemblyVector::iterator it = assemblies.begin(); it != assemblies.end(); it++)
    {
        Assembly* assembly = (*it);
        AtomVector atoms_of_assembly = assembly->GetAllAtomsOfAssembly();
        if ( local_debug > 0 )
        {
          gmml::log(__LINE__, __FILE__, gmml::INF, "Atomvector size: ");
          gmml::log(__LINE__, __FILE__, gmml::INF, std::to_string(atoms_of_assembly.size()));
        }
        for(AtomVector::iterator it1 = atoms_of_assembly.begin(); it1 != atoms_of_assembly.end(); it1++)
        {
            Atom* atom = (*it1);
            all_atoms_of_assembly.push_back(atom);
        }
    }
    ResidueVector residues = this->GetResidues();
    if ( local_debug > 0 )
    {
      gmml::log(__LINE__, __FILE__, gmml::INF, "Residue vector size: ");
      gmml::log(__LINE__, __FILE__, gmml::INF, std::to_string(residues.size()));
    }
    for(ResidueVector::iterator it = residues.begin(); it != residues.end(); it++)
    {
        Residue* residue = (*it);
        if((residue->GetName().compare("HOH") != 0) && (residue->CheckIfProtein() != true) &&
           (residue->GetName().compare("DOD") != 0))
        {
            AtomVector atoms = residue->GetAtoms();
            for(AtomVector::iterator it1 = atoms.begin(); it1 != atoms.end(); it1++)
            {
                Atom* atom = (*it1);
                // if(atom->GetDescription().find("Het;") != std::string::npos)
                // There is a bug in this code somewhere, so I commented it out
                // Won't find monosaccharides in monosaccharide only pdbs as they are in the ATOM section.
                all_atoms_of_assembly.push_back(atom);
            }
        }
    }
    return all_atoms_of_assembly;
}

MolecularModeling::ResidueVector Assembly::GetAllResiduesOfAssembly()
{
    ResidueVector all_residues_of_assembly = ResidueVector();
    AssemblyVector sub_assemblies = this->GetAssemblies();
    for(AssemblyVector::iterator it = sub_assemblies.begin(); it != sub_assemblies.end(); it++)
    {
        Assembly* sub_assembly = (*it);
        ResidueVector residues_of_assembly = sub_assembly->GetAllResiduesOfAssembly();
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

MolecularModeling::ResidueVector Assembly::GetAllProteinResiduesOfAssembly()
{
    ResidueVector protein_residues;
    ResidueVector all_residues = this->GetAllResiduesOfAssembly();
    for (ResidueVector::iterator it1 = all_residues.begin(); it1 != all_residues.end(); ++it1)
    {
        Residue *current_residue = *it1;
        if (current_residue->CheckIfProtein()==1) // the current residue is an amino acid
        {
            protein_residues.push_back(current_residue);
        }
    }
    return protein_residues;
}
GeometryTopology::CoordinateVector Assembly::GetAllCoordinates()
{
    GeometryTopology::CoordinateVector coordinates;
    for(AssemblyVector::iterator it = this->assemblies_.begin(); it != this->assemblies_.end(); it++)
    {
        Assembly* assembly = (*it);
        GeometryTopology::CoordinateVector assembly_coordinate = assembly->GetAllCoordinates();
        if(assembly_coordinate.size() == 0)
        {
//            std::cout << "Central data structure is not complete in order for generating this type of file: Missing coordinate(s)" << std::endl;
            gmml::log(__LINE__, __FILE__, gmml::ERR, "Central data structure is not complete in order for generating this type of file: Missing coordinate(s)");
            return GeometryTopology::CoordinateVector();
        }
        for(GeometryTopology::CoordinateVector::iterator it1 = assembly_coordinate.begin(); it1 != assembly_coordinate.end(); it1++)
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
//                std::cout << "Central data structure is not complete in order for generating this type of file: Missing coordinate(s)" << std::endl;
                gmml::log(__LINE__, __FILE__, gmml::ERR, "Central data structure is not complete in order for generating this type of file: Missing coordinate(s)");
                return GeometryTopology::CoordinateVector();
            }
            else
            {
                coordinates.push_back(atom->GetCoordinates()[model_index_]);
            }
        }
    }
    return coordinates;
}

GeometryTopology::CoordinateVector Assembly::GetCycleAtomCoordinates( Glycan::Monosaccharide* mono ) {
  GeometryTopology::CoordinateVector coordinates;
  for( AtomVector::iterator it1 = mono->cycle_atoms_.begin(); it1 != mono->cycle_atoms_.end(); it1++ ) {
    MolecularModeling::Atom* atom = ( *it1 );
    GeometryTopology::CoordinateVector atom_coordinates = atom->GetCoordinates();
    for( GeometryTopology::CoordinateVector::iterator it2 = atom_coordinates.begin(); it2 != atom_coordinates.end(); it2++ ) {
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
MolecularModeling::ResidueNodeVector Assembly::GetAllResidueNodesOfAssembly()
{
        return residuenodes_;
}

 //Added by ayush on 11/12/17 for molecules in assembly
MolecularModeling::MoleculeVector Assembly::GetMolecules()
{
      return molecules_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void Assembly::SetName(std::string name)
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
    std::stringstream ss;
    ss << this->name_ << "-" << assembly->GetName();
    this->name_ = ss.str();
    std::stringstream sss;
    sss << this->source_file_ << "#" << assembly->GetSourceFile();
    this->source_file_ = sss.str();
    source_file_type_ = gmml::MULTIPLE;
    assembly->SetSequenceNumber(assemblies_.size() + 1);
    std::stringstream ssss;
    ssss << this->id_ << "." << assemblies_.size() + 1;
    assembly->UpdateIds(ssss.str());
    assembly->SetId(ssss.str());
    this->assemblies_.push_back(assembly);
}

void Assembly::RemoveAssembly(Assembly *assembly)
{
    assemblies_.erase(std::remove(assemblies_.begin(), assemblies_.end(), assembly), assemblies_.end());
    //Following code added by Yao 08/20/2019. In addtion to removing assembly, also remove its residues& atoms from current assembly.
    ResidueVector all_residues = assembly->GetResidues();
    for (ResidueVector::iterator resit = all_residues.begin(); resit != all_residues.end(); resit++){
	this->RemoveResidue(*resit);
    }
}

void Assembly::UpdateIds(std::string new_id)
{
    for(AssemblyVector::iterator  it = assemblies_.begin(); it != assemblies_.end(); it++)
    {
        Assembly* assembly = *it;
        std::stringstream ss;
        ss << new_id << "." << assembly->GetId().substr(assembly->GetId().find_first_of(this->id_) + this->id_.size() + 1);
        (*it)->UpdateIds(ss.str());
        (*it)->SetId(ss.str());
    }
    for(ResidueVector::iterator it = residues_.begin(); it != residues_.end(); it++)
    {
        std::vector<std::string> id_tokens = gmml::Split((*it)->GetId(), "_");
        std::stringstream ss;
        for(unsigned int i = 0; i < id_tokens.size() - 1; i++)
        {
            ss << id_tokens.at(i) << "_";
        }
        ss << new_id;
        (*it)->SetId(ss.str());

        AtomVector atoms = (*it)->GetAtoms();
        for(AtomVector::iterator it1 = atoms.begin(); it1 != atoms.end(); it1++)
        {
            std::vector<std::string> id_tokens = gmml::Split((*it1)->GetId(), "_");
            std::stringstream ss;
            for(unsigned int i = 0; i < id_tokens.size() - 1; i++)
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

void Assembly::InsertResidue(Residue* point_of_insertion,  Residue *to_be_inserted)
{
    int distance = std::distance(residues_.begin(), std::find(residues_.begin(), residues_.end(), point_of_insertion) );
    residues_.insert(residues_.begin() + distance, to_be_inserted);
}

void Assembly::RemoveResidue(Residue *residue) // Added back in by Oliver so that Glycoprotein builder will compile.
{
//    ResidueVector newResidues = ResidueVector();
//    newResidues.resize(residues_.size() - 1); // Resizing each push_back is inefficient, set size to be current - 1.
//    for(ResidueVector::iterator it = residues_.begin(); it != residues_.end(); ++it)
//    {
//        Residue* r = *it;
//        if(r->GetId().compare(residue->GetId()) != 0)
//        {
//            //if(r->GetNode() != NULL)
//            //{
//            //    r->GetNode()->RemoveNodeNeighbor(residue); // OG to Ayush: make it so.
//            //}
//            newResidues.push_back(r);
//        }
//    }
//    this->SetResidues(newResidues);

    // Next part is cleaner way to do it, but when removing residue node too you'd have to iterate anyway
    residues_.erase(std::remove(residues_.begin(), residues_.end(), residue), residues_.end()); // Note need #include <algorithm>
}

void Assembly::SetChemicalType(std::string chemical_type)
{
    chemical_type_ = chemical_type;
}

void Assembly::SetSequenceNumber(int sequence_number)
{
    sequence_number_ = sequence_number;
}

void Assembly::SetId(std::string id)
{
    id_ = id;
}

void Assembly::SetDescription(std::string description)
{
    description_ = description;
}

void Assembly::SetSourceFile(std::string source_file)
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

void Assembly::AddNote(Glycan::Note *note)
{
    notes_.push_back(note);
}


//Added by ayush on 11/16/17 for residuenodes in assembly
void Assembly::SetResidueNodes(MolecularModeling::ResidueNodeVector residuenodes)
{
    residuenodes_.clear();
    for(MolecularModeling::ResidueNodeVector::iterator it = residuenodes.begin(); it != residuenodes.end(); it++)
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
void Assembly::MergeAssembly(Assembly *other) // Added by Oliver. He is unsure and this may well cause problems. Edit: it did.
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
    //    this->total_mass_ = gmml::dNotSet;
    //    this->center_of_geometry_ = GeometryTopology::Coordinate();
    //    this->center_of_mass_ = GeometryTopology::Coordinate();
    //    this->model_index_ = 0;
}

double Assembly::GetTotalCharge()
{
    double charge = 0;
    AtomVector all_atoms_of_assembly = this->GetAllAtomsOfAssembly();
    for(AtomVector::iterator it = all_atoms_of_assembly.begin(); it != all_atoms_of_assembly.end(); it++)
    {
        Atom* atom = *it;
        if(atom->MolecularDynamicAtom::GetCharge() != gmml::dNotSet)
            charge += atom->MolecularDynamicAtom::GetCharge();

    }
    return charge;
}

double Assembly::GetRadius()
{
    double radius = -INFINITY;
    GeometryTopology::Coordinate* geometric_center = new GeometryTopology::Coordinate();
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
        if(atom_radius == gmml::dNotSet)
            atom_radius = gmml::MINIMUM_RADIUS;
        double dist_to_edge = 0;
        if(atom_radius != gmml::dNotSet)
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

void Assembly::GetHierarchicalMapOfAssembly(HierarchicalContainmentMap &hierarchical_map, std::stringstream &index)
{
    hierarchical_map[index.str()] = this->residues_;
    if(this->assemblies_.size() == 0)
        return;
    for(AssemblyVector::iterator it = this->assemblies_.begin(); it != this->assemblies_.end(); it++)
    {
        std::stringstream i;
        i << index.str() << "." << (*it)->GetSequenceNumber();
        Assembly* assembly = (*it);
        assembly->GetHierarchicalMapOfAssembly(hierarchical_map, i);
    }
}

LibraryFileSpace::LibraryFile::ResidueMap Assembly::GetAllResiduesFromMultipleLibFilesMap(std::vector<std::string> lib_files)
{
    LibraryFileSpace::LibraryFile::ResidueMap all_residues;
    LibraryFileSpace::LibraryFile::ResidueMap residues;
    for(std::vector<std::string>::iterator it = lib_files.begin(); it != lib_files.end(); it++)
    {
        LibraryFileSpace::LibraryFile* lib_file = new LibraryFileSpace::LibraryFile(*it);
        residues = lib_file->GetResidues();
        for(LibraryFileSpace::LibraryFile::ResidueMap::iterator it1 = residues.begin(); it1 != residues.end(); it1++)
        {
            std::string lib_residue_name = (*it1).first;
            all_residues[lib_residue_name] = (*it1).second;
        }
    }
    return all_residues;
}

PrepFileSpace::PrepFile::ResidueMap Assembly::GetAllResiduesFromMultiplePrepFilesMap(std::vector<std::string> prep_files)
{
    PrepFileSpace::PrepFile::ResidueMap all_residues;
    PrepFileSpace::PrepFile::ResidueMap residues;
    for(std::vector<std::string>::iterator it = prep_files.begin(); it != prep_files.end(); it++)
    {
        PrepFileSpace::PrepFile* prep_file = new PrepFileSpace::PrepFile(*it);
        residues = prep_file->GetResidues();
        for(PrepFileSpace::PrepFile::ResidueMap::iterator it1 = residues.begin(); it1 != residues.end(); it1++)
        {
            std::string prep_residue_name = (*it1).first;
            all_residues[prep_residue_name] = (*it1).second;
        }
    }
    return all_residues;
}

gmml::ResidueNameMap Assembly::GetAllResidueNamesFromMultipleLibFilesMap(std::vector<std::string> lib_files)
{
    gmml::ResidueNameMap all_residue_names;
    std::vector<std::string> residue_names;
    for(std::vector<std::string>::iterator it = lib_files.begin(); it != lib_files.end(); it++)
    {
        LibraryFileSpace::LibraryFile* lib_file = new LibraryFileSpace::LibraryFile(*it);
        residue_names = lib_file->GetAllResidueNames();
        for(std::vector<std::string>::iterator it1 = residue_names.begin(); it1 != residue_names.end(); it1++)
        {
            std::string residue_name = (*it1);
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

MolecularModeling::ResidueNodeVector Assembly::GenerateResidueNodesInAssembly()
{

    int residue_counter=1;
    ResidueVector assembly_residues=this->GetAllResiduesOfAssembly();
    for(ResidueVector::iterator it = assembly_residues.begin(); it != assembly_residues.end(); it++)
    {
        ResidueNode* residuenode = new ResidueNode();

       Residue* residue = (*it);

        residuenode->SetResidue(residue);

        residuenode->SetId(residue_counter);

        std::string residue_id=residue->GetId();

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
    for(MolecularModeling::ResidueNodeVector::iterator it3 = residuenodes_.begin(); it3 != residuenodes_.end(); it3++)
    {
            MolecularModeling::ResidueNode* residuenode = (*it3);

            MolecularModeling::AtomVector residue_connecting_atoms = residuenode->GetResidueNodeConnectingAtoms();

            for(MolecularModeling::AtomVector::iterator it4 = residue_connecting_atoms.begin(); it4 != residue_connecting_atoms.end(); it4++)
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

    for(MolecularModeling::ResidueNodeVector::iterator it = residuenodes_.begin(); it != residuenodes_.end(); it++)
    {
         MolecularModeling::ResidueNode* residuenode = (*it);
         if(residuenode->GetIsVisited()==false)
         {
             GenerateMoleculesDFSUtil(residuenode);
             moleculecount++;
         }
    }

//    std::cout<<"Total number of molecules in Assembly is:"<<moleculecount<<std::endl;

}

void Assembly::GenerateMoleculesDFSUtil(MolecularModeling::ResidueNode* DFSresiduenode)
{
        DFSresiduenode->SetIsVisited(true);
        MolecularModeling::ResidueNodeVector residuenode_neighbors=DFSresiduenode->GetResidueNodeNeighbors();

        for(MolecularModeling::ResidueNodeVector::iterator it2 = residuenode_neighbors.begin(); it2 != residuenode_neighbors.end(); it2++)
        {
               MolecularModeling::ResidueNode* current_residuenode_neighbor = (*it2);
               if(current_residuenode_neighbor->GetIsVisited()==false)
               {
                    GenerateMoleculesDFSUtil(current_residuenode_neighbor);
               }

        }


}

//Added by ayush on 11/12/17 for molecules in assembly
void Assembly::AddMolecule(MolecularModeling::Molecule *molecule)
{
    int max_index_=molecules_.size()+1;
    molecule->SetMoleculeIndex(max_index_);
    molecules_.push_back(molecule);
}


 //Added by ayush on 04/11/18 for TopologyFix in assembly
MolecularModeling::AtomVector Assembly::GetAllBondedAtomsByStartDirection(MolecularModeling::Atom* start_atom, MolecularModeling::Atom* direction_atom , AtomVector ignore_list)
{

        if(start_atom == NULL || direction_atom == NULL){
//            std::cout<<"Start Atom or Direction Atom is null"<<std::endl;
            }
        else if(CheckIfAtomExistInAssembly(start_atom)==false||CheckIfAtomExistInAssembly(direction_atom)==false){
//            std::cout<<"Start Atom or Direction Atom does not exist in Assembly"<<std::endl;
        }else{
                bool isNeighbor=false;
                 AtomVector start_atom_neighbors = start_atom->GetNode()->GetNodeNeighbors();

                 for(AtomVector::iterator it = start_atom_neighbors.begin(); it != start_atom_neighbors.end(); it++)
                 {
                     MolecularModeling::Atom* atom = *it;
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
                          MolecularModeling::Atom* neighbor_atom = *it;
                          if(neighbor_atom->GetNode()->GetIsVisited()==false)
                          {
                                BondedAtomsByStartDirectionDFSUtil(neighbor_atom, direction_atom_neighbors, ignore_list);
                          }
                    }

                 }else{
//                     std::cout<<"Direction Atom is invalid. Not a neighbor of Start Atom."<<std::endl;
                 }
        }

        return bonded_atoms_bystartdirection_;
}

void Assembly::BondedAtomsByStartDirectionDFSUtil(MolecularModeling::Atom* DFSatom, AtomVector start_atom_neighbors, AtomVector ignore_list)
{
    bool toConsider=false; //Check if the atom is part of the start atom neighor or ignore list and set toConsider accordingly. If toConsider is true, skip the atom.

    DFSatom->GetNode()->SetIsVisited(true);

    for(AtomVector::iterator it1 = start_atom_neighbors.begin(); it1!= start_atom_neighbors.end(); it1++)
    {
        MolecularModeling::Atom* atom = *it1;
       if(atom->GetId().compare(DFSatom->GetId())==0)
       {
           toConsider = true;
       }
    }

     for(AtomVector::iterator it2 = ignore_list.begin(); it2!= ignore_list.end(); it2++)
    {
        MolecularModeling::Atom* atom = *it2;
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
            MolecularModeling::Atom* atom = *it3;
                if(atom->GetNode()->GetIsVisited()==false)
                {
                    BondedAtomsByStartDirectionDFSUtil(DFSatom, start_atom_neighbors, ignore_list);
                }
        }

        bonded_atoms_bystartdirection_.push_back(DFSatom);
    }
}


bool Assembly::CheckIfAtomExistInAssembly(MolecularModeling::Atom* toCheckAtom)
{
    bool status = false;
    AtomVector all_atoms_of_assembly = this->GetAllAtomsOfAssembly();
    for(AtomVector::iterator it = all_atoms_of_assembly.begin(); it != all_atoms_of_assembly.end(); it++)
    {
        MolecularModeling::Atom* atom = *it;
       if(atom->GetId().compare(toCheckAtom->GetId())==0)
       {
            status=true;
       }
    }
    return status;
}

 //Added by ayush on 04/16/18 for TopologyFix in assembly

GeometryTopology::CoordinateVector Assembly::GetCoordinatesFromAtomVector(AtomVector atomList, int CoordinateIndex)
{
    GeometryTopology::CoordinateVector coordinatesByIndex;
    for( AtomVector::iterator it1 = atomList.begin(); it1 != atomList.end(); it1++ ) {
      MolecularModeling::Atom* atom = ( *it1 );
      GeometryTopology::CoordinateVector atom_coordinates = atom->GetCoordinates();
      if(CoordinateIndex < (int) atom_coordinates.size())
      {
        GeometryTopology::Coordinate* coordinate = atom_coordinates[CoordinateIndex];
        coordinatesByIndex.push_back(coordinate);
      }else{
//          std::cout<<"Atom ID "<<atom->GetId()<<" does not have coordinate at index: "<<CoordinateIndex<<std::endl;
      }
    }
    return coordinatesByIndex;
}

// //Added by ayush on 06/16/18 for OffFile
void Assembly::CreateOffFileFromAssembly(std::string file_name, int CoordinateIndex)
{ 
    OffFileSpace::OffFile* off_file = new OffFileSpace::OffFile();
    off_file->Write(file_name, CoordinateIndex, this);
//    std::cout<<"end of assembly"<<std::endl;

}

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void Assembly::Print(std::ostream &out)
{
    out << "===================== " << name_ << " ============================" << std::endl;
    out << "Source file: " << source_file_ << std::endl;
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
            MolecularModeling::Residue* residue = (*it);
            residue->Print(out);
        }
    }
}

void Assembly::PrettyPrintHet(std::ostream &out)
{
    out << "===================== " << "PDB" << " ============================" << std::endl;
    out << "PDB file name: " << source_file_ << std::endl;
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
            std::string name = residue->GetName();
            if(name.compare("HOH") != 0)
                residue->PrettyPrintHet(out);
        }
    }
}

void Assembly::PrintHetResidues(std::ostream &out)
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
            std::string name = residue->GetName();
            if(name.compare("HOH") != 0)
                residue->PrintHetResidues(out);
        }
    }
}

void Assembly::PrintHetAtoms(std::ostream &out)
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
            std::string name = residue->GetName();
            if(name.compare("HOH") != 0)
                residue->PrintHetAtoms(out);
        }
    }
}

void Assembly::WriteHetResidues(std::string file_name)
{
    std::ofstream out_file;
    out_file.open(file_name.c_str());

    for(ResidueVector::iterator it = residues_.begin(); it != residues_.end(); it++)
    {
        Residue* residue = (*it);
        std::string name = residue->GetName();
        if(name.compare("HOH") != 0)
            residue->WriteHetResidues(out_file);
    }
    out_file.close();
}

void Assembly::WriteHetAtoms(std::string file_name)
{
    std::ofstream out_file;
    out_file.open(file_name.c_str());

    for(ResidueVector::iterator it = residues_.begin(); it != residues_.end(); it++)
    {
        Residue* residue = (*it);
        std::string name = residue->GetName();
        if(name.compare("HOH") != 0)
            residue->WriteHetAtoms(out_file);
    }
}
