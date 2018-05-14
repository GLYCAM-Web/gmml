#include <math.h>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <set>
#include <queue>
#include <stack>
#include <iostream>
#include <algorithm>
#include <iterator>
#include <cstdlib>

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

#include <unistd.h>
#include <errno.h>
#include <string.h>

using MolecularModeling::Assembly;

//////////////////////////////////////////////////////////
//                       FUNCTIONS                      //
//////////////////////////////////////////////////////////
gmml::GlycamResidueNamingMap Assembly::ExtractResidueGlycamNamingMap(std::vector<Glycan::Oligosaccharide*> oligosaccharides)
{
    //TODO: Done
    //Update the map
    //Key: std::string -> residue_id or atom_id
    //Value: std::vector<std::string> -> list of possible glycam naming
    gmml::GlycamResidueNamingMap pdb_glycam_residue_map = gmml::GlycamResidueNamingMap();
    //Iterates on all oligosaccharides and update the naming map
    for(std::vector<Glycan::Oligosaccharide*>::iterator it = oligosaccharides.begin(); it != oligosaccharides.end(); it++)
    {
        int index = 0;
        Glycan::Oligosaccharide* oligo = *it;
        std::string oligo_name = oligo->oligosaccharide_name_;
        //In case that there is no terminal attached to the reducing end adds a temporary terminal residue to make the sequence parser able to parse the sequence
        if(oligo->terminal_.compare("") == 0)
            oligo_name = oligo_name + "1-OH";
        CondensedSequenceSpace::CondensedSequence* condensed_sequence = new CondensedSequenceSpace::CondensedSequence(oligo_name);
        //Gets the three letter code of all carbohydrates involved in current oligosaccharide
        CondensedSequenceSpace::CondensedSequence::CondensedSequenceGlycam06ResidueTree condensed_sequence_glycam06_residue_tree = condensed_sequence->GetCondensedSequenceGlycam06ResidueTree();
        if(oligo->terminal_.compare("") != 0)
        {
            MolecularModeling::Atom* anomeric_o = NULL;
            MolecularModeling::Atom* anomeric_c = NULL;
            if(oligo->root_->cycle_atoms_.at(0) != NULL)
                anomeric_c = oligo->root_->cycle_atoms_.at(0);
            /*if(oligo->root_->side_atoms_.at(0).at(1) != NULL)
                anomeric_o = oligo->root_->side_atoms_.at(0).at(1);*/
	    //Anomeric group has now been removed from sidegroup if it shouldn't be part of oligosaccharide (e.g. NLN nitrogen). So now, find anomeric oxygen from node neighbors of anomeric carbon.
	    AtomVector anomeric_c_neighbors = anomeric_c->GetNode()->GetNodeNeighbors();
	    for (int i=0; i< anomeric_c_neighbors.size(); i++){
		if (!anomeric_c_neighbors[i]-> GetIsCycle() && (anomeric_c_neighbors[i]->GetName().substr(0,1) == "O" || anomeric_c_neighbors[i]->GetName().substr(0,1) == "S" ||
		      anomeric_c_neighbors[i]->GetName().substr(0,1) == "N") ){

		    anomeric_o = anomeric_c_neighbors[i];
		}
	    }

            if(anomeric_o != NULL && anomeric_c != NULL)
            {
                AtomVector terminal_atoms = AtomVector();
                std::string terminal = CheckTerminals(anomeric_o, terminal_atoms);
                //Separating terminal residue in glycam naming
                if(terminal.compare("") != 0)
                {
                    //TODO:
                    //Add the residue mismatch into a structure for the Ontology usage

                    for(AtomVector::iterator it1 = terminal_atoms.begin(); it1 != terminal_atoms.end(); it1++)
                    {
                        MolecularModeling::Atom* terminal_atom = *it1;
                        std::string terminal_atom_id = terminal_atom->GetId();
                        std::string terminal_residue_id = terminal_atom->GetResidue()->GetId();
                        if(pdb_glycam_residue_map.find(terminal_atom_id) == pdb_glycam_residue_map.end())
                            pdb_glycam_residue_map[terminal_atom_id] == std::vector<std::string>();
                        pdb_glycam_residue_map[terminal_atom_id].push_back(condensed_sequence_glycam06_residue_tree.at(index)->GetName());
			
                        if(pdb_glycam_residue_map.find(terminal_residue_id) == pdb_glycam_residue_map.end())
                            pdb_glycam_residue_map[terminal_residue_id] = std::vector<std::string>();
                        pdb_glycam_residue_map[terminal_residue_id].push_back(condensed_sequence_glycam06_residue_tree.at(index)->GetName());

			if ( !(terminal_atom->GetResidue()->CheckIfProtein()) )  //If this terminal atom is not part of protein,for example,NLN, then it should be in a new glycam residue, for example, ROH.
			{
			    cout << "Not protein: " << terminal_atom->GetResidue()->GetName() << ": "<< terminal_atom->GetName() << endl;
			    ResidueVector AllResiduesInAssembly = this->GetResidues();
			    Residue* OldResidueForThisAtom = terminal_atom->GetResidue();
			    string new_terminal_residue_name = condensed_sequence_glycam06_residue_tree.at(index)->GetName();
			    Residue* new_terminal_residue = new Residue(this,new_terminal_residue_name);
			    int DistanceFromStartOfResidueVector = std::distance( AllResiduesInAssembly.begin(), std::find(AllResiduesInAssembly.begin(),AllResiduesInAssembly.end(),
										OldResidueForThisAtom) );

                            this ->InsertResidue(DistanceFromStartOfResidueVector +1 ,new_terminal_residue);
			    string new_terminal_residue_id = terminal_atom->GetId(); //This new residue takes the atom id as residue id.
			    new_terminal_residue->SetId(new_terminal_residue_id);
                            new_terminal_residue->AddAtom(terminal_atom);
                            terminal_atom->SetResidue(new_terminal_residue);
                            OldResidueForThisAtom->RemoveAtom(terminal_atom);
			}
                    }

                }
            }
	}
        index++;
        this->ExtractOligosaccharideNamingMap(pdb_glycam_residue_map, oligo, condensed_sequence_glycam06_residue_tree, index);
    }

    return pdb_glycam_residue_map;
}


void Assembly::ExtractOligosaccharideNamingMap(gmml::GlycamResidueNamingMap& pdb_glycam_map, Glycan::Oligosaccharide *oligosaccharide,
                                               CondensedSequenceSpace::CondensedSequence::CondensedSequenceGlycam06ResidueTree condensed_sequence_glycam06_residue_tree, int &index)
{
    std::string name = condensed_sequence_glycam06_residue_tree.at(index)->GetName();
    //TODO: Done
    //Update to hold all possible three letter names for the specific residue_id
    std::string residue_id = oligosaccharide->root_->cycle_atoms_.at(0)->GetResidue()->GetId();
    if(pdb_glycam_map.find(residue_id) == pdb_glycam_map.end())
        pdb_glycam_map[residue_id] = std::vector<std::string>();
    pdb_glycam_map[residue_id].push_back(name);
    index++;
    //Separating SO3 and PO3 residues in glycam naming
    while(index < (int)condensed_sequence_glycam06_residue_tree.size() && condensed_sequence_glycam06_residue_tree.at(index)->GetIsDerivative())
    {
        int parent_index = condensed_sequence_glycam06_residue_tree.at(index)->GetParentId();
        int carbon_index = gmml::ConvertString<int>(condensed_sequence_glycam06_residue_tree.at(index)->GetAnomericCarbon().substr(1));
        MolecularModeling::Atom* carbon_atom = NULL;
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
            std::string n_linkage_derivative_string = CheckxC_NxO_SO3(carbon_atom, oligosaccharide->root_->cycle_atoms_str_, 'N', n_linkage_derivative_atoms);
            if(n_linkage_derivative_string.compare("xC-N-SO3") == 0)
            {
                for(AtomVector::iterator it = n_linkage_derivative_atoms.begin(); it != n_linkage_derivative_atoms.end(); it++)
                {
                    MolecularModeling::Atom* atom = *it;
                    std::string derivative_atom_id = atom->GetId();
                    if(pdb_glycam_map.find(derivative_atom_id) == pdb_glycam_map.end())
                        pdb_glycam_map[derivative_atom_id] = std::vector<std::string>();
                    pdb_glycam_map[derivative_atom_id].push_back("SO3");
                }
                std::string new_name = condensed_sequence_glycam06_residue_tree.at(parent_index)->GetName();
                std::string derivative_residue_id = oligosaccharide->root_->cycle_atoms_.at(0)->GetResidue()->GetId();
                if(pdb_glycam_map.find(derivative_residue_id) == pdb_glycam_map.end())
                    pdb_glycam_map[derivative_residue_id] = std::vector<std::string>();
                pdb_glycam_map[derivative_residue_id].push_back(gmml::ConvertT<int>(carbon_index) + new_name.substr(1));
            }
            AtomVector o_linkage_derivative_atoms = AtomVector();
            std::string o_linkage_derivative_string = CheckxC_NxO_SO3(carbon_atom, oligosaccharide->root_->cycle_atoms_str_, 'O', o_linkage_derivative_atoms);
            if(o_linkage_derivative_string.compare("xC-O-SO3") == 0)
            {
                for(AtomVector::iterator it = o_linkage_derivative_atoms.begin(); it != o_linkage_derivative_atoms.end(); it++)
                {
                    MolecularModeling::Atom* atom = *it;
                    std::string derivative_atom_id = atom->GetId();
                    if(pdb_glycam_map.find(derivative_atom_id) == pdb_glycam_map.end())
                        pdb_glycam_map[derivative_atom_id] = std::vector<std::string>();
                    pdb_glycam_map[derivative_atom_id].push_back("SO3");
                }
                std::string new_name = condensed_sequence_glycam06_residue_tree.at(parent_index)->GetName();
                std::string derivative_residue_id = oligosaccharide->root_->cycle_atoms_.at(0)->GetResidue()->GetId();
                if(pdb_glycam_map.find(derivative_residue_id) == pdb_glycam_map.end())
                    pdb_glycam_map[derivative_residue_id] = std::vector<std::string>();
                pdb_glycam_map[derivative_residue_id].push_back(gmml::ConvertT<int>(carbon_index) + new_name.substr(1));
            }
            AtomVector n_linkage_derivative_atoms_1 = AtomVector();
            std::string n_linkage_derivative_string_1 = CheckxC_NxO_PO3(carbon_atom, oligosaccharide->root_->cycle_atoms_str_, 'N', n_linkage_derivative_atoms_1);
            if(n_linkage_derivative_string_1.compare("xC-N-PO3") == 0)
            {
                for(AtomVector::iterator it = n_linkage_derivative_atoms_1.begin(); it != n_linkage_derivative_atoms_1.end(); it++)
                {
                    MolecularModeling::Atom* atom = *it;
                    std::string derivative_atom_id = atom->GetId();
                    if(pdb_glycam_map.find(derivative_atom_id) == pdb_glycam_map.end())
                        pdb_glycam_map[derivative_atom_id] = std::vector<std::string>();
                    pdb_glycam_map[derivative_atom_id].push_back("PO3");
                }
                std::string new_name = condensed_sequence_glycam06_residue_tree.at(parent_index)->GetName();
                std::string derivative_residue_id = oligosaccharide->root_->cycle_atoms_.at(0)->GetResidue()->GetId();
                if(pdb_glycam_map.find(derivative_residue_id) == pdb_glycam_map.end())
                    pdb_glycam_map[derivative_residue_id] = std::vector<std::string>();
                pdb_glycam_map[derivative_residue_id].push_back(gmml::ConvertT<int>(carbon_index) + new_name.substr(1));
            }
            AtomVector o_linkage_derivative_atoms_1 = AtomVector();
            std::string o_linkage_derivative_string_1 = CheckxC_NxO_PO3(carbon_atom, oligosaccharide->root_->cycle_atoms_str_, 'O', o_linkage_derivative_atoms_1);
            if(o_linkage_derivative_string_1.compare("xC-O-PO3") == 0)
            {
                for(AtomVector::iterator it = o_linkage_derivative_atoms_1.begin(); it != o_linkage_derivative_atoms_1.end(); it++)
                {
                    MolecularModeling::Atom* atom = *it;
                    std::string derivative_atom_id = atom->GetId();
                    if(pdb_glycam_map.find(derivative_atom_id) == pdb_glycam_map.end())
                        pdb_glycam_map[derivative_atom_id] = std::vector<std::string>();
                    pdb_glycam_map[derivative_atom_id].push_back("PO3");
                }
                std::string new_name = condensed_sequence_glycam06_residue_tree.at(parent_index)->GetName();
                std::string derivative_residue_id = oligosaccharide->root_->cycle_atoms_.at(0)->GetResidue()->GetId();
                if(pdb_glycam_map.find(derivative_residue_id) == pdb_glycam_map.end())
                    pdb_glycam_map[derivative_residue_id] = std::vector<std::string>();
                pdb_glycam_map[derivative_residue_id].push_back(gmml::ConvertT<int>(carbon_index) + new_name.substr(1));
            }
        }
        index++;
    }

    //Recursively assign glycam naming to the monosaccharides of an oligosaccharide
    for(unsigned int i = 0; i < oligosaccharide->child_oligos_.size(); i++)
    {
        this->ExtractOligosaccharideNamingMap(pdb_glycam_map, oligosaccharide->child_oligos_.at(i), condensed_sequence_glycam06_residue_tree, index);
    }
}

void Assembly::UpdateResidueName2GlycamName(gmml::GlycamResidueNamingMap residue_glycam_map, std::string prep_file)
{
    for(AssemblyVector::iterator it = this->GetAssemblies().begin(); it != this->GetAssemblies().end(); it++)
        (*it)->UpdateResidueName2GlycamName(residue_glycam_map, prep_file);

    PrepFileSpace::PrepFile* prep = new PrepFileSpace::PrepFile(prep_file);
    PrepFileSpace::PrepFile::ResidueMap prep_residues = prep->GetResidues();
    ResidueVector residues = this->GetResidues();

    ResidueVector updated_residues = ResidueVector();
    int residue_sequence_number = 0;
    ResidueVector terminal_residues = ResidueVector();
    for(ResidueVector::iterator it2 = residues.begin(); it2 != residues.end(); it2++)
    {
        Residue* residue = *it2;
        std::string residue_name = residue->GetName();
        std::string residue_id = residue->GetId();
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
            std::vector<std::string> glycam_names = residue_glycam_map[residue_id];
            gmml::GlycamAtomNameMap pdb_glycam_map = gmml::GlycamAtomNameMap();
            gmml::GlycamAtomNameMap glycam_pdb_map = gmml::GlycamAtomNameMap();
            std::string glycam_name = *glycam_names.begin();
            ResidueVector query_residues = ResidueVector();
            for(std::vector<std::string>::iterator name_it = glycam_names.begin(); name_it != glycam_names.end(); name_it++)
            {
                std::string glycam_residue_name = *name_it;
                PrepFileSpace::PrepFile::ResidueMap customized_prep_residues = PrepFileSpace::PrepFile::ResidueMap();
                customized_prep_residues[glycam_residue_name] = prep_residues[glycam_residue_name];
                PrepFileSpace::PrepFile* temp_prep_file = new PrepFileSpace::PrepFile();
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
            std::map<std::string, Residue*> residue_set = std::map<std::string, Residue*>();
            for(AtomVector::iterator it1 = atoms.begin(); it1 != atoms.end(); it1++)
            {
                MolecularModeling::Atom* atom = *it1;
                std::string atom_id = atom->GetId();
                //Terminal residue naming

                if(residue_glycam_map.find(residue_id) != residue_glycam_map.end())
                {
                    terminal_residues.at(terminal_residues.size() - 1)->SetAssembly(this);
		    if (residue_glycam_map[residue_id].size() == 1) //if residue_id and glycam names is one to one relationship
		    {
		        glycam_name = *residue_glycam_map[residue_id].begin();
		    }
		    else if (residue_glycam_map[residue_id].size() > 1) // if not one to one
		    {
			//prevent residue from being named to its attached aglycones,whose names are also in the name vector
			std::vector<std::string> known_aglycon_names = std::vector<std::string>();
			known_aglycon_names.push_back("ROH");
			known_aglycon_names.push_back("OME");
			for (std::vector<std::string>::iterator it= residue_glycam_map[residue_id].begin(); it!= residue_glycam_map[residue_id].end(); it++){
			    string possible_glycam_name = *it;
			    //if this name is not a known aglycone
			    if (std::find(known_aglycon_names.begin(),known_aglycon_names.end(),possible_glycam_name) == known_aglycon_names.end()){
		                glycam_name = possible_glycam_name ;

			    }
			}
		    }
                    terminal_residues.at(terminal_residues.size() - 1)->SetName(glycam_name);
                    terminal_residues.at(terminal_residues.size() - 1)->AddAtom(atom);
                    std::vector<std::string> residue_id_tokens = gmml::Split(residue_id, "_");
                    std::string terminal_residue_id = residue_id_tokens.at(0) + "_" + residue_id_tokens.at(1) + "_" + gmml::ConvertT<int>(residue_sequence_number) + "_"
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
                    std::vector<std::string> atom_id_tokens = gmml::Split(atom->GetId(), "_");
                    std::string atom_id_new = atom_id_tokens.at(0) + "_" + atom_id_tokens.at(1) + "_" + atom_id_tokens.at(2) + "_" + atom_id_tokens.at(3) + "_" +
                            gmml::ConvertT<int>(residue_sequence_number) + "_" + atom_id_tokens.at(5) + "_" + atom_id_tokens.at(6);
                    atom->SetId(atom_id_new);
                }
                //Non-terminal residues glycam naming
                else
                {
                    //TODO:
                    //Update to match the atoms with the corresponding prep residue to change the atom names with respect to prep file
                    //Use the map to update the atom naming
                    //Add the atom name mismatch into a structure for the Ontology usage
                    std::string prep_atom_id;
                    if(pdb_glycam_map.find(atom_id) != pdb_glycam_map.end())
                        prep_atom_id = pdb_glycam_map[atom_id];
                    else
                        prep_atom_id = atom_id;
                    glycam_name = gmml::Split(prep_atom_id,"_")[2];
                    if(residue_set.find(glycam_name) == residue_set.end())
                    {
                        residue_set[glycam_name] = (new Residue(this, glycam_name));
                        residue_sequence_number--;
                        std::vector<std::string> res_id_tokens = gmml::Split(residue_id, "_");
                        std::string res_id = glycam_name + "_" + res_id_tokens.at(1) + "_" + gmml::ConvertT<int>(residue_sequence_number) + "_"
                                + res_id_tokens.at(3) + "_" + res_id_tokens.at(4);
                        residue_set[glycam_name]->SetId(res_id);
                    }
                    std::string atom_name = atom->GetName();
                    std::string new_atom_name = gmml::Split(prep_atom_id,"_")[0];
                    std::string new_atom_id = atom_id;
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
            for(std::map<std::string, Residue*>::iterator it1 = residue_set.begin(); it1 != residue_set.end(); it1++)
            {
                updated_residues.push_back((*it1).second);
		
            }
        }
        else
            updated_residues.push_back(residue);

    }
    this->SetResidues(updated_residues);
}


bool Assembly::PatternMatching(Residue *residue, ResidueVector query_residues, gmml::GlycamAtomNameMap &pdb_glycam_map, gmml::GlycamAtomNameMap& glycam_atom_map)
{
    CreatePrunedMatchingGraph(residue, query_residues);
    for(ResidueVector::iterator it = query_residues.begin(); it != query_residues.end(); it++)
    {
        Residue* query_residue = *it;
        AtomVector lowest_degree_atoms = query_residue->GetAtomsWithLowestIntraDegree();
        std::stack<BacktrackingElements*> backtracking_stack_1 = std::stack<BacktrackingElements*>();
        for(int i = lowest_degree_atoms.size() - 1; i >= 0; i--)
            backtracking_stack_1.push(new BacktrackingElements(pdb_glycam_map, glycam_atom_map, lowest_degree_atoms, i));
        if(!backtracking_stack_1.empty())
        {
            // backtracking point
            BacktrackingElements* backtracking_point_1 = backtracking_stack_1.top();
            backtracking_stack_1.pop();
            pdb_glycam_map = backtracking_point_1->pdb_glycam_map_;
            glycam_atom_map = backtracking_point_1->glycam_atom_map_;
            lowest_degree_atoms = backtracking_point_1->atoms_;
            MolecularModeling::Atom* source_atom = lowest_degree_atoms.at(backtracking_point_1->index_);
            AtomVector source_atom_intra_neighbors = source_atom->GetNode()->GetIntraNodeNeighbors();
            std::stack<BacktrackingElements*> backtracking_stack_2 = std::stack<BacktrackingElements*>();
            for(int i = source_atom_intra_neighbors.size() - 1; i >= 0; i--)
                backtracking_stack_2.push(new BacktrackingElements(pdb_glycam_map, glycam_atom_map, source_atom_intra_neighbors, i));
            if(!backtracking_stack_2.empty())
            {
                // backtracking point
                BacktrackingElements* backtracking_point_2 = backtracking_stack_2.top();
                backtracking_stack_2.pop();
                pdb_glycam_map = backtracking_point_2->pdb_glycam_map_;
                glycam_atom_map = backtracking_point_2->glycam_atom_map_;
                source_atom_intra_neighbors = backtracking_point_2->atoms_;
                MolecularModeling::Atom* mapped_atom = source_atom_intra_neighbors.at(backtracking_point_2->index_);
                pdb_glycam_map[mapped_atom->GetId()] = source_atom->GetId();
                glycam_atom_map[source_atom->GetId()] = mapped_atom->GetId();
                AtomVector source_atom_neighbors = source_atom->GetNode()->GetNodeNeighbors();
                std::queue<MolecularModeling::Atom*> to_visit = std::queue<MolecularModeling::Atom*>();
                std::queue<MolecularModeling::Atom*> parents = std::queue<MolecularModeling::Atom*>();
                for(unsigned int i = 0; i < source_atom_neighbors.size(); i++)
                {
                    if(glycam_atom_map.find(source_atom_neighbors.at(i)->GetId()) == glycam_atom_map.end())
                    {
                        to_visit.push(source_atom_neighbors.at(i));
                        parents.push(mapped_atom);
                    }
                }
                std::stack<BacktrackingElements*> backtracking_stack_3 = std::stack<BacktrackingElements*>();
                while(!to_visit.empty())
                {
                    MolecularModeling::Atom* atom = to_visit.front();
                    MolecularModeling::Atom* parent = parents.front();
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
                    for(int i = atom_intra_neighbors_filtered_with_parent.size() - 1; i >= 0; i--)
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
                        MolecularModeling::Atom* mapped = atom_intra_neighbors_filtered_with_parent.at(backtracking_point_3->index_);
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
        MolecularModeling::Atom* query_atom = *it;
        AtomNode* query_atom_node = query_atom->GetNode();
        AtomVector query_atom_intra_neighbors = AtomVector();
        if(query_atom_node == NULL)
            query_atom_node = new AtomNode();
        for(AtomVector::iterator it1 = atoms.begin(); it1 != atoms.end(); it1++)
        {
            MolecularModeling::Atom* atom = *it1;
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
        MolecularModeling::Atom* query_atom = *it;
        AtomNode* query_atom_node = query_atom->GetNode();
        if(query_atom_node != NULL)
        {
            std::string query_atom_neighborhood_label = query_atom_node->CreateNeighboringLabel();
            AtomVector query_atom_intra_neighbors = query_atom_node->GetIntraNodeNeighbors();
            AtomVector updated_query_atom_intra_neighbors = AtomVector();
            for(AtomVector::iterator it1 = query_atom_intra_neighbors.begin(); it1 != query_atom_intra_neighbors.end(); it1++)
            {
                MolecularModeling::Atom* query_atom_intra_neighbor = *it1;
                AtomNode* query_atom_intra_neighbor_node = query_atom_intra_neighbor->GetNode();
                if(query_atom_intra_neighbor_node != NULL)
                {
                    std::string query_atom_intra_neighbor_neighborhood_label = query_atom_intra_neighbor_node->CreateNeighboringLabel();
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
