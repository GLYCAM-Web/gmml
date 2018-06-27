#include <math.h>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <set>
#include <queue>
#include <stack>
#include <sstream>
#include <string>
#include <iostream>

#include "../../../includes/InputSet/PdbFileSpace/inputfile.hpp"
#include "../../../includes/MolecularModeling/assembly.hpp"
#include "../../../includes/MolecularModeling/residue.hpp"
#include "../../../includes/MolecularModeling/residuenode.hpp"
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
#include "../../../includes/MolecularMetadata/GLYCAM/glycam06residueinfo.hpp"
#include "../../../includes/GeometryTopology/grid.hpp"
#include "../../../includes/GeometryTopology/cell.hpp"

#include <unistd.h>
#include <errno.h>
#include <string.h>
#include <iostream>
#include <cstdlib> //Yao: to use the exit() function. Right now, rather than throwing exceptions, I use exit().

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

Assembly::TemplateAssembly* Assembly::BuildTemplateAssemblyFromPrepFile (CondensedSequenceSpace::CondensedSequence::CondensedSequenceGlycam06ResidueTree& glycam06_residue_tree, 
										PrepFileSpace::PrepFile* prep_file)
{
    //Sorting query_residue_names to remove duplicate residue names.
    std::vector<std::string> query_residue_names = std::vector<std::string>();
    for (CondensedSequenceSpace::CondensedSequence::CondensedSequenceGlycam06ResidueTree::iterator it = glycam06_residue_tree.begin(); it != glycam06_residue_tree.end(); it++){

	CondensedSequenceSpace::CondensedSequenceGlycam06Residue* glycam06_residue = *it;
	std::string residue_name = glycam06_residue->GetName();
	query_residue_names.push_back(residue_name);
    }	

    std::set<std::string> query_residue_names_unique = std::set<std::string>();
    query_residue_names_unique.insert(query_residue_names.begin(), query_residue_names.end());
    TemplateAssembly* template_assembly = new TemplateAssembly();

    ResidueVector template_assembly_residues = ResidueVector();
    PrepFileSpace::PrepFile::ResidueMap prep_residue_map = prep_file->GetResidues();
    std::vector<std::string> all_prep_residue_names = std::vector<std::string>();
    //Getting all residue names in prep file
    for (PrepFileSpace::PrepFile::ResidueMap::iterator it = prep_residue_map.begin(); it != prep_residue_map.end(); it++)
    {
	std::string prep_residue_name = it-> first;
	all_prep_residue_names.push_back(prep_residue_name);
    }
    //Check each prep file residue, if their name match query residue name, build assembly residue from this prep residue, and add assembly residue to a template assembly.
    for (std::set<std::string>::iterator it = query_residue_names_unique.begin(); it != query_residue_names_unique.end(); it++)
    {
	if (std::find(all_prep_residue_names.begin(), all_prep_residue_names.end(), *it) != all_prep_residue_names.end() )
 	{
	    PrepFileSpace::PrepFileResidue* prep_residue = prep_residue_map[*it];
	    Residue* assembly_residue = new Residue();
	    assembly_residue->BuildResidueFromPrepFileResidue(prep_residue); 
	    template_assembly_residues.push_back(assembly_residue);
	}
    }

    template_assembly->SetResidues(template_assembly_residues);

    //Tag cycle atoms and sidechain atoms based on MolecularMetadata lookup map
    for( ResidueVector::iterator it = template_assembly_residues.begin(); it != template_assembly_residues.end(); it++ ) {
	MolecularModeling::Residue* residue = (*it);
	std::string residue_name = residue->GetName();
	gmml::AtomVector all_atoms = residue->GetAtoms();
	std::pair<std::multimap<std::string, std::string>::const_iterator, std::multimap<std::string, std::string>::const_iterator> key_range = gmml::MolecularMetadata::GLYCAM::Glycam06NamesToTypesLookupMap.			equal_range(residue_name);

	if (key_range.first == key_range.second){
	    std::cout << "Warning: no match exists in metadata map for template residue " << residue_name << std::endl;
	    std::cout << "Cannot tag ring/sidechain atoms for this residue. This might affect accuracy of setting geometry." << std::endl;
	}
	else{
	    std::vector<std::string> all_types = std::vector<std::string>();
	    for (std::multimap<std::string, std::string>::const_iterator it = key_range.first; it != key_range.second; it++){
		std::string type =  it->second;
		all_types.push_back(type);
	    }
	    bool is_sugar = false;
	    std::string ring_atom_str = "";
	    //If template residue is a aldofuranose, ring atom is: C1,C2,C3,C4,O4
	    if (std::find(all_types.begin(), all_types.end(), "aldose") != all_types.end() && std::find(all_types.begin(), all_types.end(), "furanose") != all_types.end()){
		ring_atom_str = "C1_C2_C3_C4_O4";
		is_sugar = true;
	    }
	    //aldopyranose: C1,C2,C3,C4,C5,O5
	    if (std::find(all_types.begin(), all_types.end(), "aldose") != all_types.end() && std::find(all_types.begin(), all_types.end(), "pyranose") != all_types.end()){
		ring_atom_str = "C1_C2_C3_C4_C5_O5";
		is_sugar = true;
	    }
	    //ketofuranose: C2,C3,C4,C5,O5
	    if (std::find(all_types.begin(), all_types.end(), "ketose") != all_types.end() && std::find(all_types.begin(), all_types.end(), "furanose") != all_types.end()){
		ring_atom_str = "C2_C3_C4_C5_O5";
		is_sugar = true;
	    }
	    //ketopyranose: C2,C3,C4,C5,C6,O5 (shouldnt' it be O6?)
	    if (std::find(all_types.begin(), all_types.end(), "ketose") != all_types.end() && std::find(all_types.begin(), all_types.end(), "pyranose") != all_types.end()){
		ring_atom_str = "C2_C3_C4_C5_C6_O6";
		is_sugar = true;
	    }
	    //ulosonate: C2,C3,C4,C5,C6,O5 , just like ketopyranose (shouldn't it be O6?)
	    if (std::find(all_types.begin(), all_types.end(), "ulosonate") != all_types.end()){
		ring_atom_str = "C2_C3_C4_C5_C6_O6";
		is_sugar = true;
	    }
		for (gmml::AtomVector::iterator it2 = all_atoms.begin(); it2 != all_atoms.end(); it2++){
		    MolecularModeling::Atom* atom = *it2;
		    std::string atom_name = atom->GetName();
		    if (ring_atom_str.find(atom_name) != std::string::npos){
			atom -> SetIsCycle(true);
		    }
		    else if (is_sugar){
			atom -> SetIsSideChain(true);
		    }
		}
	}
    }
//testing
/*
for (ResidueVector::iterator it3 = template_assembly_residues.begin(); it3 != template_assembly_residues.end(); it3++){
    std:: cout << "Residues: " << (*it3)->GetName() <<std::endl;
    gmml::AtomVector atoms = (*it3)->GetAtoms();
    for (gmml::AtomVector::iterator it4 = atoms.begin(); it4 != atoms.end();it4++){
	if ((*it4) ->GetIsCycle()){
	    std::cout << "is cycle: " << (*it4)->GetName() << std::endl;
	} 
	if ((*it4)->GetIsSideChain()){
	    //std::cout << "is sidechain: " << (*it4)->GetName()<< std::endl;
	}
    }
}
*/
//testing

    return template_assembly;
    
}

std::map<int, std::pair<CondensedSequenceSpace::CondensedSequenceGlycam06Residue*, MolecularModeling::Residue*> > 
Assembly::ConvertCondensedSequence2AssemblyResidues(CondensedSequenceSpace::CondensedSequence::CondensedSequenceGlycam06ResidueTree& glycam06_residue_tree, TemplateAssembly* template_assembly)
{
    
    ResidueVector all_template_residues = template_assembly->GetResidues();
    ResidueVector newly_added_residues = ResidueVector();
    std::map<CondensedSequenceSpace::CondensedSequenceGlycam06Residue*, CondensedSequenceSpace::CondensedSequenceGlycam06Residue*> condensed_sequence_child_parent_map = 
	std::map<CondensedSequenceSpace::CondensedSequenceGlycam06Residue*, CondensedSequenceSpace::CondensedSequenceGlycam06Residue*>();

    std::map<int, std::pair<CondensedSequenceSpace::CondensedSequenceGlycam06Residue*, MolecularModeling::Residue*> > index_condensed_sequence_assembly_residue_map = 
	std::map<int, std::pair<CondensedSequenceSpace::CondensedSequenceGlycam06Residue*, MolecularModeling::Residue*> >();
    std::map<Atom*, AtomNode*> new_atom_template_node_map = std::map<Atom*, AtomNode*>();

    for (unsigned int i = 0; i< glycam06_residue_tree.size(); i++) 
    {
	CondensedSequenceSpace::CondensedSequenceGlycam06Residue* condensed_sequence_residue = glycam06_residue_tree[i];
        std::string condensed_sequence_residue_name = condensed_sequence_residue->GetName();
	int residue_serial_number = i;
	for (unsigned int j = 0; j < all_template_residues.size(); j++){
	    //Copy residues
	    if (condensed_sequence_residue_name == all_template_residues[j]->GetName())
	    {
		Residue* template_residue = all_template_residues[j];
		Residue* assembly_residue = new Residue();
		assembly_residue->SetIsSugarDerivative (condensed_sequence_residue->GetIsDerivative());
		assembly_residue->SetNode(NULL);	//When residue object is constructed, its node attribute is not set to NULL, causing it to contain garbage node address.So I need to do this.
		this->AddResidue(assembly_residue);
		assembly_residue->SetAssembly(this);
		assembly_residue->SetName(template_residue->GetName());
		std::stringstream id_stream;
		id_stream << template_residue->GetName() << "_" << residue_serial_number << "_" << template_residue->GetId() << std::endl;
		//std::string residue_id = template_residue->GetName() + "_" + serial_number_stream.str() + "_" + template_residue->GetId();
		assembly_residue->SetId (id_stream.str());
		index_condensed_sequence_assembly_residue_map[residue_serial_number].first = condensed_sequence_residue;
		index_condensed_sequence_assembly_residue_map[residue_serial_number].second = assembly_residue;
		newly_added_residues.push_back(assembly_residue);
	
		//Copy atoms in residue
		int atom_serial_number = 0;
		AtomVector all_template_atoms = template_residue->GetAtoms();
		for (unsigned int k = 0; k < all_template_atoms.size(); k++)
		{
		    Atom* template_atom = all_template_atoms[k];
		    Atom* template_atom_copy = new Atom();
		    assembly_residue->AddAtom(template_atom_copy);
		    template_atom_copy->SetResidue(assembly_residue);
		    template_atom_copy->SetName(template_atom->GetName());
		    template_atom_copy->SetNaming(template_atom->GetNaming());
	    //Attention: SetAtomType()function is overloaded as MolecularModeling::Atom::SetAtomType() and MolecularModeling::MolecularDynamicAtom::SetAtomType(). You don't really know which one to use.
	    //Likiwise, GetAtomType() is also overloaded.
            //In my situation, I called MolecularModeling::MolecularModelingAtom::SetAtomType(), but later called MolecularModeling::Atom::GetAtomType(). The result is empty.
            //We need to talk about this later
		    template_atom_copy->MolecularDynamicAtom::SetAtomType(template_atom->MolecularDynamicAtom::GetAtomType());
		    template_atom_copy->SetCharge(template_atom->GetCharge());
	    //Create coordinate object for new atom
		    GeometryTopology::Coordinate* copy_coordinate = new GeometryTopology::Coordinate(template_atom->GetCoordinates().at(0));
		    template_atom_copy->AddCoordinate(copy_coordinate);
		    std::stringstream atom_serial_number_stream;
		    atom_serial_number_stream << atom_serial_number;
		    std::string atom_id = template_atom->GetName() + "_" + atom_serial_number_stream.str() + "_" + template_atom->GetId();
		    template_atom_copy->SetId(atom_id);
	    //Tag cycle/sidechain atoms
		    template_atom_copy->MolecularModeling::OligoSaccharideDetectionAtom::SetIsCycle(template_atom->GetIsCycle());
		    template_atom_copy->MolecularModeling::OligoSaccharideDetectionAtom::SetIsSideChain(template_atom->GetIsSideChain());

		    AtomVector template_head_atoms = template_residue->GetHeadAtoms();
		    if (std::find(template_head_atoms.begin(), template_head_atoms.end(), template_atom) != template_head_atoms.end() ){
			assembly_residue->AddHeadAtom(template_atom);
		    }
		    AtomVector template_tail_atoms = template_residue->GetTailAtoms();
		    if (std::find(template_tail_atoms.begin(), template_tail_atoms.end(), template_atom) != template_tail_atoms.end() ){
			assembly_residue->AddTailAtom(template_atom);
		    }
		    AtomNode* template_atom_node = template_atom ->GetNode(); 
		    new_atom_template_node_map [template_atom_copy] = template_atom_node;
		    atom_serial_number++;
		}
		//Copy atom nodes
		for (std::map<Atom*, AtomNode*>::iterator it = new_atom_template_node_map.begin(); it != new_atom_template_node_map.end(); it++){
		    Atom* new_atom = it->first;
		    AtomNode* template_node = it ->second;
		    AtomNode* new_atom_node = new AtomNode();
		    new_atom_node->SetAtom(new_atom);
		    new_atom->SetNode(new_atom_node);
		    std::vector<std::string> all_node_neighbor_names = std::vector<std::string>();
		    AtomVector template_node_neighbors = template_node->GetNodeNeighbors();
		    for (AtomVector::iterator it2 = template_node_neighbors.begin(); it2 != template_node_neighbors.end(); it2++){
			all_node_neighbor_names.push_back((*it2)->GetName());
		    }
		    AtomVector atoms_in_new_residue = new_atom->GetResidue()->GetAtoms();
		    for (AtomVector::iterator it2 = atoms_in_new_residue.begin(); it2 != atoms_in_new_residue.end(); it2++){
			Atom* new_atom_in_residue = *it2;
			std::string atom_name = new_atom_in_residue->GetName(); 
			if (std::find(all_node_neighbor_names.begin(), all_node_neighbor_names.end(), atom_name) != all_node_neighbor_names.end()){
			    new_atom_node->AddNodeNeighbor(new_atom_in_residue);
			}

		    }

		}

		
	    }
	}

    }
//testing

    return index_condensed_sequence_assembly_residue_map;

}//ConvertCondensedSequence2AssemblyResidues

void Assembly::SetGlycam06ResidueBonding (std::map<int, std::pair<CondensedSequenceSpace::CondensedSequenceGlycam06Residue*, MolecularModeling::Residue*> >& index_condensed_sequence_assembly_residue_map)
{
    //index_condensed_sequence_assembly_residue_map : map <int ,pair <06 residue, assembly residue> >. The int key is to maintain the order of residue in 06 residue tree.
    std::map<CondensedSequenceSpace::CondensedSequenceGlycam06Residue*, MolecularModeling::Residue*> condensed_assembly_residue_map = 
	std::map<CondensedSequenceSpace::CondensedSequenceGlycam06Residue*, MolecularModeling::Residue*>();

    std::map<CondensedSequenceSpace::CondensedSequenceGlycam06Residue*, CondensedSequenceSpace::CondensedSequence::CondensedSequenceGlycam06ResidueTree> condensed_residue_chilren_map = 
	std::map<CondensedSequenceSpace::CondensedSequenceGlycam06Residue*, CondensedSequenceSpace::CondensedSequence::CondensedSequenceGlycam06ResidueTree> ();

    std::map<CondensedSequenceSpace::CondensedSequenceGlycam06Residue*, CondensedSequenceSpace::CondensedSequence::CondensedSequenceGlycam06ResidueTree> condensed_residue_parent_map = 
	std::map<CondensedSequenceSpace::CondensedSequenceGlycam06Residue*, CondensedSequenceSpace::CondensedSequence::CondensedSequenceGlycam06ResidueTree> ();

    //Setting condensed_assembly_residue_map, and initiate the values of the other two maps above to empty.
     for (std::map<int, std::pair<CondensedSequenceSpace::CondensedSequenceGlycam06Residue*, MolecularModeling::Residue*> >::iterator it = index_condensed_sequence_assembly_residue_map.begin() ;
         it != index_condensed_sequence_assembly_residue_map.end(); it++){
	condensed_assembly_residue_map[it->second.first] = it->second.second;
	condensed_residue_chilren_map[it->second.first] = CondensedSequenceSpace::CondensedSequence::CondensedSequenceGlycam06ResidueTree();
	condensed_residue_parent_map[it->second.first] = CondensedSequenceSpace::CondensedSequence::CondensedSequenceGlycam06ResidueTree();
    }
    //Associating a condensed sequence residue with its children and/or parent, using two nested iterations through index_condensed_sequence_assembly_residue_map;
    for (std::map<int, std::pair<CondensedSequenceSpace::CondensedSequenceGlycam06Residue*, MolecularModeling::Residue*> >::iterator it = index_condensed_sequence_assembly_residue_map.begin() ;
	 it != index_condensed_sequence_assembly_residue_map.end(); it++){
	int index1 = it -> first;
	CondensedSequenceSpace::CondensedSequenceGlycam06Residue* glycam_06_residue = it->second.first;
	for (std::map<int, std::pair<CondensedSequenceSpace::CondensedSequenceGlycam06Residue*, MolecularModeling::Residue*> >::iterator it2 = index_condensed_sequence_assembly_residue_map.begin() ;
	     it2 != index_condensed_sequence_assembly_residue_map.end(); it2++){
	    int index2 = it2 -> first;
	    CondensedSequenceSpace::CondensedSequenceGlycam06Residue* glycam_06_residue_II = it2->second.first;
	    if (glycam_06_residue_II -> GetParentId() == index1){	//If "glycam_06_residue_II" 's parent is "glycam_06_residue", or the latter is the child of the former
		condensed_residue_chilren_map[glycam_06_residue].push_back(glycam_06_residue_II);
	    }
	    if (glycam_06_residue -> GetParentId() == index2){		//If "glycam_06_residue" 's parent is ""glycam_06_residue_II""
		condensed_residue_parent_map[glycam_06_residue].push_back(glycam_06_residue_II);
	    }
	    
	}
    }

    //Will set head and tail atoms based on child-parent relationship in the big for loop below.So now empty whatever the default setting is.
    for (std::map<CondensedSequenceSpace::CondensedSequenceGlycam06Residue*, MolecularModeling::Residue*>::iterator map_it = condensed_assembly_residue_map.begin(); 
	 map_it != condensed_assembly_residue_map.end(); map_it++){
	MolecularModeling::Residue* residue = map_it->second;
	residue->SetHeadAtoms(gmml::AtomVector());
	residue->SetTailAtoms(gmml::AtomVector());
    }

    //Set Assembly Residue Nodes based on information extracted from condensed sequence.
    for (std::map<CondensedSequenceSpace::CondensedSequenceGlycam06Residue*, MolecularModeling::Residue*> ::iterator it = condensed_assembly_residue_map.begin();
		it != condensed_assembly_residue_map.end(); it++){
	CondensedSequenceSpace::CondensedSequenceGlycam06Residue* condensed_residue = it->first;
	MolecularModeling::Residue* corresponding_assembly_residue = it->second;	
	ResidueVector corresponding_assembly_residue_neighbors = ResidueVector();
	//find all assembly residue children for "corresponding_assembly_residue"
	CondensedSequenceSpace::CondensedSequence::CondensedSequenceGlycam06ResidueTree condensed_residue_children = condensed_residue_chilren_map[condensed_residue];
	for (unsigned int i = 0; i < condensed_residue_children.size(); i++){

	    MolecularModeling::Residue* corresponding_assembly_residue_child = condensed_assembly_residue_map[condensed_residue_children[i] ];
	    corresponding_assembly_residue_neighbors.push_back(corresponding_assembly_residue_child); 
	}

	//find all assembly residue parents for "corresponing_assembly_residue"
	CondensedSequenceSpace::CondensedSequence::CondensedSequenceGlycam06ResidueTree condensed_residue_parents = condensed_residue_parent_map[condensed_residue];
	for (unsigned int i = 0; i < condensed_residue_parents.size(); i++){
	    MolecularModeling::Residue* corresponding_assembly_residue_parent = condensed_assembly_residue_map[condensed_residue_parents[i] ];
	    corresponding_assembly_residue_neighbors.push_back(corresponding_assembly_residue_parent); 
	}

	//Create Residue Nodes and add Node neighbors
	ResidueNode* corresponding_assembly_residue_node = NULL;
	if (corresponding_assembly_residue->GetNode() == NULL){
	    corresponding_assembly_residue_node = new ResidueNode();
	    corresponding_assembly_residue_node->SetResidue(corresponding_assembly_residue);
	    corresponding_assembly_residue->SetNode(corresponding_assembly_residue_node);
	}
	else{
	    corresponding_assembly_residue_node = corresponding_assembly_residue->GetNode();
	}
	
	ResidueNodeVector existing_assembly_residue_neighbor_nodes = corresponding_assembly_residue_node->GetResidueNodeNeighbors();
	//Iterate through all neighboring residues, and bond this neighbor to curent assembly residue.
	for (unsigned int i = 0; i < corresponding_assembly_residue_neighbors.size(); i++){
	    MolecularModeling::Residue* neighbor = corresponding_assembly_residue_neighbors[i];
	    ResidueNode* neighbor_node = NULL;
	    if (neighbor->GetNode() == NULL){
		neighbor_node = new ResidueNode();
		neighbor_node->SetResidue(neighbor);
		neighbor->SetNode(neighbor_node);
	    }
	    else{
		neighbor_node = neighbor->GetNode();
	    }

	    if (std::find(existing_assembly_residue_neighbor_nodes.begin(), existing_assembly_residue_neighbor_nodes.end(), neighbor_node) == existing_assembly_residue_neighbor_nodes.end() ){
	        corresponding_assembly_residue_node->AddResidueNodeNeighbor(neighbor_node);
	    }

	    ResidueNodeVector existing_neighbor_neighbor_nodes = neighbor_node->GetResidueNodeNeighbors();
	    if (std::find(existing_neighbor_neighbor_nodes.begin(),existing_neighbor_neighbor_nodes.end(), corresponding_assembly_residue_node) == existing_neighbor_neighbor_nodes.end() ){
	    	neighbor_node->AddResidueNodeNeighbor(corresponding_assembly_residue_node);
	    }
	}//for

    }//for
  
    //Set Residue Node Connecting Atoms, and set a bond between connecting connecting atoms in adjacent residue nodes.
    for (std::map<CondensedSequenceSpace::CondensedSequenceGlycam06Residue*, MolecularModeling::Residue*>::iterator it =condensed_assembly_residue_map.begin();
                it != condensed_assembly_residue_map.end(); it++){
	CondensedSequenceSpace::CondensedSequenceGlycam06Residue* condensed_residue = it->first;
        MolecularModeling::Residue* corresponding_assembly_residue = it->second;
	CondensedSequenceSpace::CondensedSequence::CondensedSequenceGlycam06ResidueTree condensed_residue_children = condensed_residue_chilren_map[condensed_residue];

	std::map<MolecularModeling::Atom*, MolecularModeling::Atom*> parent_tail_child_head_map = std::map<MolecularModeling::Atom*, MolecularModeling::Atom*>();
	//find out the name of parent oxygen
	for (unsigned int i = 0; i < condensed_residue_children.size(); i++){
	    CondensedSequenceSpace::CondensedSequenceGlycam06Residue* condensed_residue_child = condensed_residue_children[i];
	    std::string parent_oxygen_name = condensed_residue_child -> GetParentOxygen();
	    std::string anomeric_carbon_name = condensed_residue_child-> GetAnomericCarbon();
	    MolecularModeling::Residue* assembly_residue_child = condensed_assembly_residue_map[condensed_residue_child];

	    gmml::AtomVector all_atoms_in_assembly_residue = corresponding_assembly_residue->GetAtoms();
	    gmml::AtomVector all_atoms_in_child_assembly_residue = assembly_residue_child->GetAtoms();
	    for (unsigned int j = 0; j < all_atoms_in_assembly_residue.size(); j++){
		if (all_atoms_in_assembly_residue[j]->GetName() == parent_oxygen_name){
		    for (unsigned int k = 0; k < all_atoms_in_child_assembly_residue.size(); k++){
			if (all_atoms_in_child_assembly_residue[k]->GetName() == anomeric_carbon_name){
			    MolecularModeling::Atom* parent_tail_atom = all_atoms_in_assembly_residue[j];
			    MolecularModeling::Atom* child_head_atom = all_atoms_in_child_assembly_residue[k];
			    parent_tail_child_head_map[parent_tail_atom] = child_head_atom;
			}
		    }
		}
	    }
	}

	//Also set a bond between child anomeric carbon and linking parent tail atom.
	for (std::map<MolecularModeling::Atom*, MolecularModeling::Atom*>::iterator tail_head_mapit = parent_tail_child_head_map.begin(); 
		tail_head_mapit != parent_tail_child_head_map.end(); tail_head_mapit++){
	    MolecularModeling::Atom* parent_tail_atom = tail_head_mapit->first;
	    MolecularModeling::Atom* child_head_atom = tail_head_mapit->second;
	    //Add parent tail to child head's node neighbors.
	    AtomNode* head_atom_node = NULL;
	    if (child_head_atom->GetNode() == NULL){
	        head_atom_node = new AtomNode();
	        head_atom_node->SetAtom(child_head_atom);
	        child_head_atom->SetNode(head_atom_node);
	    }
	    else
		head_atom_node = child_head_atom->GetNode();

	    gmml::AtomVector existing_head_atom_neighbors = head_atom_node->GetNodeNeighbors();
	    if (std::find (existing_head_atom_neighbors.begin(), existing_head_atom_neighbors.end(), parent_tail_atom) == 
	         existing_head_atom_neighbors.end() ){
		head_atom_node->AddNodeNeighbor(parent_tail_atom);
	    }
			    
	    //Add child head to parent tail's node neighbors.
	    AtomNode*tail_atom_node = NULL;
	    if (parent_tail_atom->GetNode() == NULL){
                tail_atom_node = new AtomNode();
                tail_atom_node->SetAtom(parent_tail_atom);
                parent_tail_atom->SetNode(tail_atom_node);
            }
            else
                tail_atom_node = parent_tail_atom->GetNode();

            gmml::AtomVector existing_tail_atom_neighbors = tail_atom_node->GetNodeNeighbors();
            if (std::find (existing_tail_atom_neighbors.begin(), existing_tail_atom_neighbors.end(), child_head_atom) == 
                existing_tail_atom_neighbors.end() ){
                tail_atom_node->AddNodeNeighbor(child_head_atom);
            }
	    //If child residue is a derivative, must remove a hydrogen neighbor of the parent tail atom.
	    if (child_head_atom->GetResidue()->GetIsSugarDerivative()){
	        gmml::AtomVector parent_tail_neighbors = tail_atom_node->GetNodeNeighbors();
	        for (unsigned int k = 0; k < parent_tail_neighbors.size(); k++){
	            if (parent_tail_neighbors[k]->GetName().substr(0,1) == "H"){
	                parent_tail_atom -> GetResidue() -> RemoveAtom(parent_tail_neighbors[k]);
	            }
	        }
	    }
	
	}

	//Add all tail/head atoms to corresponding residue's tail/head atoms, as well as corresponding residueNode connecting atoms
	for (std::map<MolecularModeling::Atom*, MolecularModeling::Atom*>::iterator it2 = parent_tail_child_head_map.begin(); it2 != parent_tail_child_head_map.end(); it2++){
	    MolecularModeling::Atom* parent_tail_atom = it2->first;
	    MolecularModeling::Atom* child_head_atom = it2->second;
	    parent_tail_atom->GetResidue()->AddTailAtom(parent_tail_atom);
	    parent_tail_atom->GetResidue()->GetNode()->AddResidueNodeConnectingAtom(parent_tail_atom);
	    child_head_atom->GetResidue()->AddHeadAtom(child_head_atom);
	    child_head_atom->GetResidue()->GetNode()->AddResidueNodeConnectingAtom(child_head_atom);
	}
    }//for

}//SetGlycam06ResidueBonding

void Assembly::RecursivelySetGeometry (MolecularModeling::Residue* parent_residue, int&a)
{
    gmml::AtomVector all_tail_atoms = parent_residue->GetTailAtoms();
    std::cout << "parent: " << parent_residue->GetName() <<std::endl;
    for (unsigned int i = 0; i < all_tail_atoms.size(); i++){
	MolecularModeling::Atom* tail_atom = all_tail_atoms[i];
	gmml::AtomVector tail_atom_neighbors = tail_atom->GetNode()->GetNodeNeighbors();
	gmml::AtomVector all_atoms_in_residue = parent_residue->GetAtoms();
	for (unsigned int j = 0; j < tail_atom_neighbors.size(); j++){
	    MolecularModeling::Atom* neighbor_atom = tail_atom_neighbors[j];
	    //if a neighbor is outside of a parent residue, it must be the head atom of a child residue.There should only exist one such atom, otherwise, something is amiss.
	    if (std::find(all_atoms_in_residue.begin(), all_atoms_in_residue.end(), neighbor_atom) == all_atoms_in_residue.end()){
	        a++;
		MolecularModeling::Atom* head_atom_of_child_residue = neighbor_atom;
		MolecularModeling::Residue* child_residue = head_atom_of_child_residue->GetResidue();
		gmml::AtomVector all_atoms_in_child_residue = child_residue->GetAtoms();
		//Right now, all residues are at the position of the template residue. That is, they are all around the orgin and stacked upon each other.
		//SetResidueResidueBondDistance function: takes a pair of parent tail/child head atoms as argument. This function keeps the parent residue intact,but
		//finds out the new position of child head atom, and move atoms of child residue accordingly.(i.e. grafting)
		this->SetResidueResidueBondDistance(tail_atom, head_atom_of_child_residue);
		if (a >= 999)
		break;
		
		std::cout << "Tail:"<< tail_atom->GetName() <<"--" << tail_atom->GetResidue()->GetName() << std::endl;
		std::cout << "Head:" << head_atom_of_child_residue->GetName() << "--" << head_atom_of_child_residue->GetResidue()->GetName() << std::endl;
		//Set C-O(tail atom)-C(head atom ) angle to 120 deg. The first C is the neighbor of tail atom that is not the head atom && not a hydrogen
		//The only exception is ROH, where you have to use that hydrogen
		//Right now,rely on the first letter of atom name to determine element type (if hydrogen or not). Better solution is the rule class.
		MolecularModeling::Atom* non_hydrogen_tail_atom_neighbor = NULL;
		for (unsigned int k = 0; k < tail_atom_neighbors.size(); k++){
		    if (parent_residue->GetName() == "ROH" && tail_atom_neighbors[k] != head_atom_of_child_residue)
		   	non_hydrogen_tail_atom_neighbor = tail_atom_neighbors[k];

		    else if (tail_atom_neighbors[k] != head_atom_of_child_residue && tail_atom_neighbors[k]->GetName().substr(0,1) != "H"){
		        non_hydrogen_tail_atom_neighbor = tail_atom_neighbors[k];
		    }
		}

		//Exit if atom 1 cannot be be found
		if (non_hydrogen_tail_atom_neighbor == NULL){ 
		    std::cout << "SetAngle: Cannot find atom 1.Skipping" << std::endl;
		}
		//If atom 1 is identified, proceed with setting angle
		else{
		    const double angle_to_set = 109.4;	//assuming sp3 tetrahedral
		    this->SetAngle(non_hydrogen_tail_atom_neighbor, tail_atom, head_atom_of_child_residue, angle_to_set);
		}
		
		//If child is a derivative,or parent is an  aglycon, skip setting dihedrals
		if (child_residue->GetIsSugarDerivative()  || parent_residue->GetName() == "OME" || parent_residue->GetName() == "TBT"){
		    std::cout << "Warning: skip setting dihedrals for derivative/aglycon residue " << child_residue->GetName() << std::endl;
		}
		//Otherwise,attempt setting dihedrals
		else{
		    //Set Dihedral phi:
		    MolecularModeling::Atom* phi_atom_1 = non_hydrogen_tail_atom_neighbor;	
		    MolecularModeling::Atom* phi_atom_2 = tail_atom;	
		    MolecularModeling::Atom* phi_atom_3 = head_atom_of_child_residue;	
		    MolecularModeling::Atom* phi_atom_4 = NULL;

		    gmml::AtomVector head_atom_neighbors = head_atom_of_child_residue->GetNode()->GetNodeNeighbors();
		    std::string anomeric_carbon_index_str = head_atom_of_child_residue->GetName().substr(1,1); //The "1" in C1, the "2" in C2
		    std::stringstream s1;
		    s1 << anomeric_carbon_index_str;
		    int anomeric_carbon_index;
		    s1 >> anomeric_carbon_index;
		    for (unsigned int k = 0; k< head_atom_neighbors.size(); k++){
			MolecularModeling::Atom* neighbor = head_atom_neighbors[k];
			
			//Phi_atom_4 should be C anomeric plus 1. For example, if anomeric is C1, then C2. If anomeirc is C2, then C3. Set this dihedral to 180 deg.
			//For now: If child residue is furanose,keto pyranose i.e. sialic acid etc, phi_atom_4 is C1. If is other pyronaoses,phi_atom_4 should be a hydrogen
			if ( neighbor->GetIsCycle() && neighbor != tail_atom){
			    std::string neighbor_index_str = neighbor->GetName().substr(1,1);
			    std::stringstream s2;
			    s2 << neighbor_index_str;
			    int neighbor_index;
			    s2 >> neighbor_index;
			    if(neighbor_index == anomeric_carbon_index + 1){
			        phi_atom_4 = head_atom_neighbors[k];
			    }
			}
		    }//for
		
			//If phi_atom_4 cannot be found, exit.
		    if (phi_atom_4 == NULL || phi_atom_3 == NULL || phi_atom_2 == NULL || phi_atom_1 == NULL ){
			std::cout << "SetPhiDihedral: cannot find all four phi atoms. Skipping." << std::endl;
		    }
			//If phi_atom_4 can be found, set phi.
		    else{
			std::cout << "phi: " << phi_atom_1->GetName() << "-" << phi_atom_2->GetName() << "-" << phi_atom_3->GetName() << "-" << phi_atom_4->GetName() << std::endl;
			const double dihedral_phi = 180.0;
			this->SetDihedral(phi_atom_1, phi_atom_2, phi_atom_3, phi_atom_4, dihedral_phi);
		    }
		    //Set psi, for example: psi H4-C4-O4-C1 to 0 deg
		    MolecularModeling::Atom* psi_atom_4 = head_atom_of_child_residue;
		    MolecularModeling::Atom* psi_atom_3 = tail_atom;
		    MolecularModeling::Atom* psi_atom_2 = non_hydrogen_tail_atom_neighbor;
		    MolecularModeling::Atom* psi_atom_1 = NULL;
		    gmml::AtomVector psi_atom_2_neighbors = psi_atom_2->GetNode()->GetNodeNeighbors();
		    //Psi atom 1 should be a exocyclic (normally hydrogen) neighbor of the neighbor of tail atom (neighbor of neighbor of tail oxygen), for example: psi H4-C4-O4-C1
		    for (unsigned int k = 0; k < psi_atom_2_neighbors.size(); k++){
			if (!psi_atom_2_neighbors[k]->GetIsCycle() && psi_atom_2_neighbors[k] != psi_atom_3){
			    psi_atom_1 = psi_atom_2_neighbors[k];
			}
		    }

		    if (psi_atom_1 == NULL || psi_atom_2 == NULL || psi_atom_3 == NULL || psi_atom_4 == NULL ){
		        std::cout << "SetPsiDihedral: cannot find all four psi atoms. Skipping." << std::endl;
		    }
		    else {
			std::cout << "psi: " << psi_atom_1->GetName() << "-" << psi_atom_2->GetName() << "-" << psi_atom_3->GetName() << "-" << psi_atom_4->GetName() << std::endl;
		        const double dihedral_psi = 0.0;
		        this->SetDihedral(psi_atom_1, psi_atom_2, psi_atom_3, psi_atom_4, dihedral_psi);
		    }
		
		    //Set Omega (if exists) , for example, C4-C5-C6-O6. Set this to 180 deg
			
		    MolecularModeling::Atom* omega_atom_4 = psi_atom_3;
		    MolecularModeling::Atom* omega_atom_3 = psi_atom_2;
		    MolecularModeling::Atom* omega_atom_2 = NULL;
		    MolecularModeling::Atom* omega_atom_1 = NULL;
		    //omega atom 3 should be a exocyclic non-hydrogen(probably carbon) atom that 's connects to the atoms that connects to tail atom (neighbor of neighbor of tail oxygen)
		    //If such an atom is already on the ring, then there is no omega angle. If it is exocyclic, then omega exists.
		    if (omega_atom_3->GetIsCycle()){
			std::cout << "Tail atom is directly attached to ring atom. So omega does not exist. Skip setting omega." << std::endl;
		    }
		    else{
			//Choose omega atom 2 from the neighbors of omega atom 3. It can't be omega_atom_4, and it shouldn't be a hydrogen
		        gmml::AtomVector omega_atom_3_neighbors = omega_atom_3->GetNode()->GetNodeNeighbors();
		        for (gmml::AtomVector::iterator atom_it4 = omega_atom_3_neighbors.begin(); atom_it4 != omega_atom_3_neighbors.end(); atom_it4++){
			    MolecularModeling::Atom* neighbor = *atom_it4;
			    if (neighbor->GetName().substr(0,1) != "H" && neighbor != omega_atom_4){
			        omega_atom_2 = neighbor;
			    }
		        }
		    }
		    //Once omega atom 2 is identified, get its non-hydrogen node neighbor, this should be omega atom 1.
		    if (omega_atom_2 !=NULL){	//if there is such an exocyclic atom, then omega atom 2 exists.
			gmml::AtomVector omega_atom_2_neighbors = omega_atom_2->GetNode()->GetNodeNeighbors();
			for (gmml::AtomVector::iterator atom_it5 = omega_atom_2_neighbors.begin(); atom_it5 != omega_atom_2_neighbors.end(); atom_it5++){
			    MolecularModeling::Atom* neighbor = *atom_it5;
			    if (neighbor->GetName().substr(0,1) == "H" && neighbor != omega_atom_3){
				omega_atom_1 = neighbor;
			    }
			}
		    }
		    //If all four omega atoms exist, set omega to -60 deg, assuming gt.
		    if (omega_atom_4 == NULL || omega_atom_3 == NULL || omega_atom_2 == NULL || omega_atom_1 == NULL){
			std::cout << "SetOmegaDihedral: cannot find all four omega atoms. Skipping." << std::endl;
		    }
		    else {
			std::cout << "omega: " << omega_atom_1->GetName() << "-" << omega_atom_2->GetName() << "-" << omega_atom_3->GetName() << "-" << omega_atom_4->GetName() << std::endl;
			const double dihedral_omega = -60.0;
			this ->SetDihedral(omega_atom_1, omega_atom_2, omega_atom_3, omega_atom_4, dihedral_omega);
		    }
		   	
		}//else Done setting phi,phi, omega(if exists)
		//Start new recursion
		MolecularModeling::Residue* new_parent_residue = child_residue;
		//Leave a blank line for easy visualization
		std::cout << std::endl;
		this-> RecursivelySetGeometry(new_parent_residue, a);
	    }//if
	}//for

    }//for
}

//New BuildAssemblyFromCondensedSequence() created by Yao on 06/25/2018. This will replace the old version below.
void Assembly::BuildAssemblyFromCondensedSequence(std::string condensed_sequence, PrepFileSpace::PrepFile* prep_file)
{
    CondensedSequenceSpace::CondensedSequence sequence (condensed_sequence);
    CondensedSequenceSpace::CondensedSequence::CondensedSequenceGlycam06ResidueTree glycam06_residues = sequence.GetCondensedSequenceGlycam06ResidueTree();

    MolecularModeling::Assembly::TemplateAssembly* template_assembly = this-> BuildTemplateAssemblyFromPrepFile (glycam06_residues, prep_file);

    std::map<int, std::pair<CondensedSequenceSpace::CondensedSequenceGlycam06Residue*, MolecularModeling::Residue*> > glycam06_assembly_residue_map = 
	this -> ConvertCondensedSequence2AssemblyResidues (glycam06_residues, template_assembly);

    this -> SetGlycam06ResidueBonding (glycam06_assembly_residue_map);
    
    int a=0; //temporary, for testing
    for (std::map<int,std::pair<CondensedSequenceSpace::CondensedSequenceGlycam06Residue*, MolecularModeling::Residue*> >::iterator it = 
	glycam06_assembly_residue_map.begin(); it != glycam06_assembly_residue_map.end(); it++){

	CondensedSequenceSpace::CondensedSequenceGlycam06Residue* glycam_06_res = it->second.first;
	MolecularModeling::Residue* corresponding_assembly_residue = it->second.second;
        //The atom and the only atom without a parent is the absolute parent(terminal).
        if (glycam_06_res->GetParentId() == -1 ){
            MolecularModeling::Residue* root = corresponding_assembly_residue;
            this->RecursivelySetGeometry(root, a);
	    break;
        }
    }
//test
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
    PdbFileSpace::PdbFile* pdb_file;
    try
    {
        gmml::log(__LINE__, __FILE__, gmml::INF, "Reading PDB file into PdbFileSpace::PdbFile structure ...");
        pdb_file = new PdbFileSpace::PdbFile(pdb_file_path);
    }
    catch(PdbFileSpace::PdbFileProcessingException &ex)
    {
        std::cout << "Generating PdbFileSpace::PdbFile structure from " << pdb_file_path << "failed." << std::endl;
    }
    this->BuildAssemblyFromPdbFile(pdb_file, amino_lib_files, glycam_lib_files, other_lib_files, prep_files, parameter_file);
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
        
        this->input_file_ = pdb_file;
        // std::stringstream out_stream;
        // this->input_file_->PrintOntology(out_stream);
        // std::cout << out_stream.str();
        // int testPoly = this->input_file_->GetMasterCard()->GetNumRemark();
        // std::cout << testPoly << std::endl << std::endl;


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
                float atom_b_factor = atom->GetAtomTempretureFactor();
                new_atom->SetBFactor(atom_b_factor);
                // std::stringstream test;
                // test << atom_b_factor;
                //gmml::log(__LINE__, __FILE__, gmml::INF, test.str());
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
