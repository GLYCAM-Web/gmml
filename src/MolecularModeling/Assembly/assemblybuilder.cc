// Please keep these calls in alphabetical order to avoid duplication
#include <cstdlib> //Yao: to use the exit() function. Right now, rather than throwing exceptions, I use exit().
#include <errno.h>
#include <fstream>
#include <iostream>
#include <iterator> //Added by Yao 07/23/2018
#include <map>      //Added by Yao 01/29/2020
#include <algorithm> //Added by Yao 01/29/2020, for std::find()
#include <math.h>
#include <queue>
#include <regex> // Added by OG 2019.09.05
#include <set>
#include <sstream>
#include <stack>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <string.h>
#include <unistd.h>

// Please keep these calls in directory & alphabetical order to avoid duplication
//    Please put files not ib a subdirectory first
//    Then add subdirectories in alphabetical order
//    Apply those two rules recursively
#include "../../../includes/utils.hpp"
#include "../../../includes/common.hpp"
#include "../../../includes/GeometryTopology/cell.hpp"
#include "../../../includes/GeometryTopology/grid.hpp"
#include "../../../includes/GeometryTopology/rotation.hpp"
#include "../../../includes/GeometryTopology/ResidueLinkages/residue_linkage.hpp"
#include "../../../includes/InputSet/CondensedSequenceSpace/condensedsequence.hpp"
#include "../../../includes/InputSet/CondensedSequenceSpace/condensedsequenceglycam06residue.hpp"
#include "../../../includes/InputSet/CondensedSequenceSpace/condensedsequenceresidue.hpp"
#include "../../../includes/InputSet/CoordinateFileSpace/coordinatefile.hpp"
#include "../../../includes/InputSet/CoordinateFileSpace/coordinatefileprocessingexception.hpp"
#include "../../../includes/InputSet/PdbFileSpace/inputfile.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbatomcard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbatomsection.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbconnectsection.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbfile.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbfileprocessingexception.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbheterogenatomsection.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdblinkcard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdblinkcardresidue.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdblinksection.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbmastercard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbmodelcard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbmodelresidueset.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbmodelsection.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbremarksection.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbtitlesection.hpp"
#include "../../../includes/InputSet/PdbqtFileSpace/pdbqtatom.hpp"
#include "../../../includes/InputSet/PdbqtFileSpace/pdbqtatomcard.hpp"
#include "../../../includes/InputSet/PdbqtFileSpace/pdbqtfile.hpp"
#include "../../../includes/InputSet/PdbqtFileSpace/pdbqtfileprocessingexception.hpp"
#include "../../../includes/InputSet/PdbqtFileSpace/pdbqtmodel.hpp"
#include "../../../includes/InputSet/PdbqtFileSpace/pdbqtmodelcard.hpp"
#include "../../../includes/InputSet/PdbqtFileSpace/pdbqtmodelresidueset.hpp"
#include "../../../includes/InputSet/TopologyFileSpace/topologyangle.hpp"
#include "../../../includes/InputSet/TopologyFileSpace/topologyangletype.hpp"
#include "../../../includes/InputSet/TopologyFileSpace/topologyassembly.hpp"
#include "../../../includes/InputSet/TopologyFileSpace/topologyatom.hpp"
#include "../../../includes/InputSet/TopologyFileSpace/topologyatompair.hpp"
#include "../../../includes/InputSet/TopologyFileSpace/topologybond.hpp"
#include "../../../includes/InputSet/TopologyFileSpace/topologybondtype.hpp"
#include "../../../includes/InputSet/TopologyFileSpace/topologydihedral.hpp"
#include "../../../includes/InputSet/TopologyFileSpace/topologydihedraltype.hpp"
#include "../../../includes/InputSet/TopologyFileSpace/topologyfile.hpp"
#include "../../../includes/InputSet/TopologyFileSpace/topologyfileprocessingexception.hpp"
#include "../../../includes/InputSet/TopologyFileSpace/topologyresidue.hpp"
#include "../../../includes/MolecularMetadata/GLYCAM/glycam06residueinfo.hpp"
#include "../../../includes/MolecularModeling/assembly.hpp"
#include "../../../includes/MolecularModeling/atom.hpp"
#include "../../../includes/MolecularModeling/atomnode.hpp"
#include "../../../includes/MolecularModeling/overlaps.hpp"  //Added by Yao 07/27/2018
#include "../../../includes/MolecularModeling/residue.hpp"
#include "../../../includes/MolecularModeling/residuenode.hpp"
#include "../../../includes/ParameterSet/LibraryFileSpace/libraryfile.hpp"
#include "../../../includes/ParameterSet/LibraryFileSpace/libraryfileatom.hpp"
#include "../../../includes/ParameterSet/LibraryFileSpace/libraryfileprocessingexception.hpp"
#include "../../../includes/ParameterSet/LibraryFileSpace/libraryfileresidue.hpp"
#include "../../../includes/ParameterSet/PrepFileSpace/prepfile.hpp"
#include "../../../includes/ParameterSet/PrepFileSpace/prepfileatom.hpp"
#include "../../../includes/ParameterSet/PrepFileSpace/prepfileprocessingexception.hpp"
#include "../../../includes/ParameterSet/PrepFileSpace/prepfileresidue.hpp"
#include "../../../includes/ParameterSet/ParameterFileSpace/parameterfile.hpp"
#include "../../../includes/ParameterSet/ParameterFileSpace/parameterfileangle.hpp"
#include "../../../includes/ParameterSet/ParameterFileSpace/parameterfileatom.hpp"
#include "../../../includes/ParameterSet/ParameterFileSpace/parameterfilebond.hpp"
#include "../../../includes/ParameterSet/ParameterFileSpace/parameterfiledihedral.hpp"
#include "../../../includes/ParameterSet/ParameterFileSpace/parameterfiledihedralterm.hpp"


using MolecularModeling::Assembly;
using MolecularModeling::AssemblyVector;
using MolecularModeling::Residue;
using MolecularModeling::ResidueVector;
using MolecularModeling::Atom;
using MolecularModeling::AtomVector;

//////////////////////////////////////////////////////////
//                       FUNCTIONS                      //
//////////////////////////////////////////////////////////
bool Assembly::CheckCondensedSequenceSanity(std::string sequence, CondensedSequenceSpace::CondensedSequence::CondensedSequenceGlycam06ResidueTree& prep_residues)
{
    if (sequence.empty()) // Is sequence an empty string?
    {
//        std::cout << "The input sequence (" << sequence << ") is empty" << std::endl;
        return false;
    }
    std::vector<char> badChars = {'\'', '(', ')'}; // Does sequence contain any characters it doesn't? Can expand this list as needed.
    for (auto badChar : badChars)
    {
        if (sequence.find(badChar) != std::string::npos)
        {
//            std::cout << "The input sequence (" << sequence << ") is contains a bad character: " << badChar << std::endl;
            return false;
        }
    }
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
//                std::cout << "The input sequence (" << sequence << ") is not valid" << std::endl;
                return false;
            }
        }
    }
    catch(std::exception ex)
    {
//        std::cout << "The input sequence (" << sequence << ") is not valid" << std::endl;
        return false;
    }

//    std::cout << "The input sequence (" << sequence << ") is valid" << std::endl;
    return true;
}

Assembly::TemplateAssembly* Assembly::BuildTemplateAssemblyFromPrepFile (CondensedSequenceSpace::CondensedSequence::CondensedSequenceGlycam06ResidueTree& glycam06_residue_tree,
                                                                         PrepFileSpace::PrepFile* prep_file)
{
    std::vector<std::string> query_residue_names = std::vector<std::string>();
    for (CondensedSequenceSpace::CondensedSequence::CondensedSequenceGlycam06ResidueTree::iterator it = glycam06_residue_tree.begin(); it != glycam06_residue_tree.end(); it++)
    {
        CondensedSequenceSpace::CondensedSequenceGlycam06Residue* glycam06_residue = *it;
        std::string residue_name = glycam06_residue->GetName();
        query_residue_names.push_back(residue_name);
    }
    return this->BuildTemplateAssemblyFromPrepFile (query_residue_names, prep_file);

}

Assembly::TemplateAssembly* Assembly::BuildTemplateAssemblyFromPrepFile (std::vector<std::string>& query_residue_names, PrepFileSpace::PrepFile* prep_file)
{
    //Sorting query_residue_names to remove duplicate residue names.
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
        else
        {
//            std::cout << "Warning: Cannot find prep file residue with the name: " << *it << std::endl;
//            std::cout << "Either the prep file doesn't have it, or this is not a regular residue at all. For example,\"Deoxy\" means removing parent oxygen it attaches to. " << std::endl;
        }
    }

    template_assembly->SetResidues(template_assembly_residues);

    //Tag cycle atoms and sidechain atoms based on MolecularMetadata lookup map
    for( ResidueVector::iterator it = template_assembly_residues.begin(); it != template_assembly_residues.end(); it++ )
    {
        Residue* residue = (*it);
        std::string residue_name = residue->GetName();
        AtomVector all_atoms = residue->GetAtoms();
        gmml::MolecularMetadata::GLYCAM::Glycam06NamesToTypesLookupContainer glycam06_NamesToTypesMetadata;
        std::vector<std::string> all_types = glycam06_NamesToTypesMetadata.GetTypesForResidue(residue_name);
        //      std::pair<std::multimap<std::string, std::string>::const_iterator, std::multimap<std::string, std::string>::const_iterator> key_range = gmml::MolecularMetadata::GLYCAM::Glycam06NamesToTypesLookupMap.equal_range(residue_name);
        //        if (key_range.first == key_range.second)
        //        {
        //            std::cout << "Warning: no match exists in metadata map for template residue " << residue_name << std::endl;
        //            std::cout << "Cannot tag ring/sidechain atoms for this residue. This might affect accuracy of setting geometry." << std::endl;
        //        }
        //        else
        //        {
        //            std::vector<std::string> all_types = std::vector<std::string>();
        //            for (std::multimap<std::string, std::string>::const_iterator it = key_range.first; it != key_range.second; it++)
        //            {
        //                std::string type =  it->second;
        //                all_types.push_back(type);
        //            }
        if (all_types.empty())
        {
//            std::cout << "Warning: no match exists in metadata map for template residue " << residue_name << "\n"
//                      << "Cannot tag ring/sidechain atoms for this residue. This might affect accuracy of setting geometry." << std::endl;
        }
        else
        {
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
            for (AtomVector::iterator it2 = all_atoms.begin(); it2 != all_atoms.end(); it2++){
                Atom* atom = *it2;
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
        AtomVector atoms = (*it3)->GetAtoms();
        for (AtomVector::iterator it4 = atoms.begin(); it4 != atoms.end();it4++){
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

std::map<int, std::pair<CondensedSequenceSpace::CondensedSequenceGlycam06Residue*, Residue*> >
Assembly::ConvertCondensedSequence2AssemblyResidues(CondensedSequenceSpace::CondensedSequence::CondensedSequenceGlycam06ResidueTree& glycam06_residue_tree, TemplateAssembly* template_assembly)
{
    
    ResidueVector all_template_residues = template_assembly->GetResidues();
    ResidueVector newly_added_residues = ResidueVector();
    std::map<CondensedSequenceSpace::CondensedSequenceGlycam06Residue*, CondensedSequenceSpace::CondensedSequenceGlycam06Residue*> condensed_sequence_child_parent_map =
            std::map<CondensedSequenceSpace::CondensedSequenceGlycam06Residue*, CondensedSequenceSpace::CondensedSequenceGlycam06Residue*>();

    std::map<int, std::pair<CondensedSequenceSpace::CondensedSequenceGlycam06Residue*, Residue*> > index_condensed_sequence_assembly_residue_map =
            std::map<int, std::pair<CondensedSequenceSpace::CondensedSequenceGlycam06Residue*, Residue*> >();
    std::map<Atom*, AtomNode*> new_atom_template_node_map = std::map<Atom*, AtomNode*>();

    int atom_serial_number = 1;
    for (unsigned int i = 0; i< glycam06_residue_tree.size(); i++)
    {
        CondensedSequenceSpace::CondensedSequenceGlycam06Residue* condensed_sequence_residue = glycam06_residue_tree[i];
        std::string condensed_sequence_residue_name = condensed_sequence_residue->GetName();
        if (condensed_sequence_residue_name != "Deoxy"){
            int residue_serial_number = i;
            for (unsigned int j = 0; j < all_template_residues.size(); j++){
                //Copy residues
                if (condensed_sequence_residue_name == all_template_residues[j]->GetName())
                {
                    Residue* template_residue = all_template_residues[j];
                    Residue* assembly_residue = new Residue();
                    //Set the first residue in the tree to be the aglycon.
                    if (i == 0){
                        assembly_residue->SetIsAglycon(true);
                    }
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
                    AtomVector all_template_atoms = template_residue->GetAtoms();
                    for (unsigned int k = 0; k < all_template_atoms.size(); k++)
                    {
                        Atom* template_atom = all_template_atoms[k];
                        Atom* template_atom_copy = new Atom();
                        assembly_residue->AddAtom(template_atom_copy);
                        template_atom_copy->SetResidue(assembly_residue);
                        template_atom_copy->SetName(template_atom->GetName());
                        template_atom_copy->SetNaming(template_atom->GetNaming());
                        template_atom_copy->SetElementSymbol(template_atom->GetElementSymbol());
                        //Attention: SetAtomType()function is overloaded as Atom::SetAtomType() and MolecularModeling::MolecularDynamicAtom::SetAtomType(). You don't really know which one to use.
                        //Likiwise, GetAtomType() is also overloaded.
                        //In my situation, I called MolecularModeling::MolecularModelingAtom::SetAtomType(), but later called Atom::GetAtomType(). The result is empty.
                        //We need to talk about this later
                        template_atom_copy->MolecularDynamicAtom::SetAtomType(template_atom->MolecularDynamicAtom::GetAtomType());
                        template_atom_copy->SetCharge(template_atom->GetCharge());
                        //Create coordinate object for new atom
                        GeometryTopology::Coordinate* copy_coordinate = new GeometryTopology::Coordinate(template_atom->GetCoordinates().at(0));
                        template_atom_copy->AddCoordinate(copy_coordinate);
                        std::stringstream atom_id_stream;
                        atom_id_stream << template_atom->GetName() << "_" << atom_serial_number << "_" << template_residue->GetName() << "_" << "A" << "_" << " " << "_" << "?_" << "?_" << " " << std::endl;
                        template_atom_copy->SetId(atom_id_stream.str());
                        //Tag cycle/sidechain atoms
                        template_atom_copy->MolecularModeling::OligoSaccharideDetectionAtom::SetIsCycle(template_atom->GetIsCycle());
                        template_atom_copy->MolecularModeling::OligoSaccharideDetectionAtom::SetIsSideChain(template_atom->GetIsSideChain());

                        AtomVector template_head_atoms = template_residue->GetHeadAtoms();
                        if (std::find(template_head_atoms.begin(), template_head_atoms.end(), template_atom) != template_head_atoms.end() ){
                            assembly_residue->AddHeadAtom(template_atom_copy);
                        }
                        AtomVector template_tail_atoms = template_residue->GetTailAtoms();
                        if (std::find(template_tail_atoms.begin(), template_tail_atoms.end(), template_atom) != template_tail_atoms.end() ){
                            assembly_residue->AddTailAtom(template_atom_copy);
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
        }//if not "Deoxy" false residue
    }
    //Dealing with deoxy residue
    for (unsigned int i = 0; i< glycam06_residue_tree.size(); i++){
        if (glycam06_residue_tree[i]->GetName() == "Deoxy"){
            CondensedSequenceSpace::CondensedSequenceGlycam06Residue* deoxy_derivative = glycam06_residue_tree[i];
            int parent_index = deoxy_derivative->GetParentId();
            std::string parent_oxygen_name = deoxy_derivative->GetParentOxygen();
            Residue* assembly_residue_parent = index_condensed_sequence_assembly_residue_map[parent_index].second;
            AtomVector parent_atoms = assembly_residue_parent->GetAtoms();

            for (unsigned int j = 0; j< parent_atoms.size(); j++){
                if (parent_atoms[j]->GetName() == parent_oxygen_name){
                    Atom* parent_oxygen_to_remove = parent_atoms[j];
                    AtomVector parent_neighbors = parent_oxygen_to_remove->GetNode()->GetNodeNeighbors();
                    Atom* hydrogen_to_add = new Atom();
                    AtomNode* new_hydrogen_node = new AtomNode();
                    hydrogen_to_add->SetNode(new_hydrogen_node);
                    new_hydrogen_node->SetAtom(hydrogen_to_add);
                    hydrogen_to_add->SetResidue(assembly_residue_parent);
                    assembly_residue_parent->AddAtom(hydrogen_to_add);
                    Atom* parent_neighbor_hydrogen_neighbor = NULL;
                    Atom* parent_neighbor_ring_carbon = NULL;

                    for (unsigned int k = 0; k < parent_neighbors.size(); k++){
                        if (parent_neighbors[k]->GetIsCycle()){  //When the neighbor of the oxygen to be removed is the ring carbon
                            parent_neighbor_ring_carbon = parent_neighbors[k];
                            parent_neighbor_ring_carbon->GetNode()->AddNodeNeighbor(hydrogen_to_add);
                            parent_neighbor_ring_carbon->SetCharge(parent_neighbor_ring_carbon->GetCharge() + parent_oxygen_to_remove->GetCharge()); //Add oxygen charge to ring carbon
                        }
                        else {  //When the neighbor of the oxyben to be removed is the attached hydroxyl hydrogen
                            parent_neighbor_hydrogen_neighbor = parent_neighbors[k];
                            parent_neighbor_ring_carbon->SetCharge(parent_neighbor_ring_carbon->GetCharge() + parent_neighbor_hydrogen_neighbor->GetCharge()); //Add side group H charge to ring carbon
                            assembly_residue_parent->RemoveAtom(parent_neighbor_hydrogen_neighbor);
                        }
                    }

                    hydrogen_to_add->AddCoordinate(parent_oxygen_to_remove->GetCoordinates().at(0));
                    hydrogen_to_add->SetName("Hd");
                    hydrogen_to_add->MolecularDynamicAtom::SetAtomType("H1");
                    hydrogen_to_add->SetCharge(0.0000);  //According to GLYCAM06, the charge of aliphatic oxygen is zero.
                    hydrogen_to_add->SetElementSymbol("H");
                    std::string new_hydrogen_atom_id = parent_oxygen_to_remove->GetId();
                    new_hydrogen_atom_id.replace(0,parent_oxygen_name.size(),"Hd");
                    hydrogen_to_add->SetId(new_hydrogen_atom_id);
                    AtomVector original_parent_tail_atoms = assembly_residue_parent->GetTailAtoms();
                    AtomVector updated_parent_tail_atoms = AtomVector();
                    for (unsigned int k = 0; k < original_parent_tail_atoms.size(); k++){
                        if (original_parent_tail_atoms[k] != parent_oxygen_to_remove){
                            Atom* remaining_tail_atom = original_parent_tail_atoms[k];
                            updated_parent_tail_atoms.push_back(remaining_tail_atom);
                        }
                    }
                    assembly_residue_parent->SetTailAtoms(updated_parent_tail_atoms);
                    assembly_residue_parent->RemoveAtom(parent_oxygen_to_remove);
                }
            }
        }
    }

    return index_condensed_sequence_assembly_residue_map;

}//ConvertCondensedSequence2AssemblyResidues

void Assembly::SetGlycam06ResidueBonding (std::map<int, std::pair<CondensedSequenceSpace::CondensedSequenceGlycam06Residue*, Residue*> >& index_condensed_sequence_assembly_residue_map)
{
    //index_condensed_sequence_assembly_residue_map : map <int ,pair <06 residue, assembly residue> >. The int key is to maintain the order of residue in 06 residue tree.
    std::map<CondensedSequenceSpace::CondensedSequenceGlycam06Residue*, Residue*> condensed_assembly_residue_map =
            std::map<CondensedSequenceSpace::CondensedSequenceGlycam06Residue*, Residue*>();

    std::map<CondensedSequenceSpace::CondensedSequenceGlycam06Residue*, CondensedSequenceSpace::CondensedSequence::CondensedSequenceGlycam06ResidueTree> condensed_residue_chilren_map =
            std::map<CondensedSequenceSpace::CondensedSequenceGlycam06Residue*, CondensedSequenceSpace::CondensedSequence::CondensedSequenceGlycam06ResidueTree> ();

    std::map<CondensedSequenceSpace::CondensedSequenceGlycam06Residue*, CondensedSequenceSpace::CondensedSequence::CondensedSequenceGlycam06ResidueTree> condensed_residue_parent_map =
            std::map<CondensedSequenceSpace::CondensedSequenceGlycam06Residue*, CondensedSequenceSpace::CondensedSequence::CondensedSequenceGlycam06ResidueTree> ();

    //Setting condensed_assembly_residue_map, and initiate the values of the other two maps above to empty.
    for (std::map<int, std::pair<CondensedSequenceSpace::CondensedSequenceGlycam06Residue*, Residue*> >::iterator it = index_condensed_sequence_assembly_residue_map.begin() ;
         it != index_condensed_sequence_assembly_residue_map.end(); it++){
        condensed_assembly_residue_map[it->second.first] = it->second.second;
        condensed_residue_chilren_map[it->second.first] = CondensedSequenceSpace::CondensedSequence::CondensedSequenceGlycam06ResidueTree();
        condensed_residue_parent_map[it->second.first] = CondensedSequenceSpace::CondensedSequence::CondensedSequenceGlycam06ResidueTree();
    }
    //Associating a condensed sequence residue with its children and/or parent, using two nested iterations through index_condensed_sequence_assembly_residue_map;
    for (std::map<int, std::pair<CondensedSequenceSpace::CondensedSequenceGlycam06Residue*, Residue*> >::iterator it = index_condensed_sequence_assembly_residue_map.begin() ;
         it != index_condensed_sequence_assembly_residue_map.end(); it++){
        int index1 = it -> first;
        CondensedSequenceSpace::CondensedSequenceGlycam06Residue* glycam_06_residue = it->second.first;
        for (std::map<int, std::pair<CondensedSequenceSpace::CondensedSequenceGlycam06Residue*, Residue*> >::iterator it2 = index_condensed_sequence_assembly_residue_map.begin() ;
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
    //Moreover, prep file defaults can sometimes be WRONG!!!!!!!!!!!!!!!!!!!!!
    for (std::map<CondensedSequenceSpace::CondensedSequenceGlycam06Residue*, Residue*>::iterator map_it = condensed_assembly_residue_map.begin();
         map_it != condensed_assembly_residue_map.end(); map_it++){
        Residue* residue = map_it->second;
        residue->SetHeadAtoms(AtomVector());
        residue->SetTailAtoms(AtomVector());
    }

    //Set Assembly Residue Nodes based on information extracted from condensed sequence.
    for (std::map<CondensedSequenceSpace::CondensedSequenceGlycam06Residue*, Residue*> ::iterator it = condensed_assembly_residue_map.begin();
         it != condensed_assembly_residue_map.end(); it++){
        CondensedSequenceSpace::CondensedSequenceGlycam06Residue* condensed_residue = it->first;
        Residue* corresponding_assembly_residue = it->second;
        ResidueVector corresponding_assembly_residue_neighbors = ResidueVector();
        //find all assembly residue children for "corresponding_assembly_residue"
        CondensedSequenceSpace::CondensedSequence::CondensedSequenceGlycam06ResidueTree condensed_residue_children = condensed_residue_chilren_map[condensed_residue];
        for (unsigned int i = 0; i < condensed_residue_children.size(); i++){

            Residue* corresponding_assembly_residue_child = condensed_assembly_residue_map[condensed_residue_children[i] ];
            corresponding_assembly_residue_neighbors.push_back(corresponding_assembly_residue_child);
        }

        //find all assembly residue parents for "corresponing_assembly_residue"
        CondensedSequenceSpace::CondensedSequence::CondensedSequenceGlycam06ResidueTree condensed_residue_parents = condensed_residue_parent_map[condensed_residue];
        for (unsigned int i = 0; i < condensed_residue_parents.size(); i++){
            Residue* corresponding_assembly_residue_parent = condensed_assembly_residue_map[condensed_residue_parents[i] ];
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
            Residue* neighbor = corresponding_assembly_residue_neighbors[i];
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
    for (std::map<CondensedSequenceSpace::CondensedSequenceGlycam06Residue*, Residue*>::iterator it =condensed_assembly_residue_map.begin();
         it != condensed_assembly_residue_map.end(); it++){
        CondensedSequenceSpace::CondensedSequenceGlycam06Residue* condensed_residue = it->first;
        Residue* corresponding_assembly_residue = it->second;
        CondensedSequenceSpace::CondensedSequence::CondensedSequenceGlycam06ResidueTree condensed_residue_children = condensed_residue_chilren_map[condensed_residue];

        std::map<Atom*, Atom*> parent_tail_child_head_map = std::map<Atom*, Atom*>();
        //find out the name of parent oxygen
        for (unsigned int i = 0; i < condensed_residue_children.size(); i++){
            CondensedSequenceSpace::CondensedSequenceGlycam06Residue* condensed_residue_child = condensed_residue_children[i];
            std::string parent_oxygen_name = condensed_residue_child -> GetParentOxygen();
            std::string anomeric_carbon_name = condensed_residue_child-> GetAnomericCarbon();
            Residue* assembly_residue_child = condensed_assembly_residue_map[condensed_residue_child];

            AtomVector all_atoms_in_assembly_residue = corresponding_assembly_residue->GetAtoms();
            AtomVector all_atoms_in_child_assembly_residue = assembly_residue_child->GetAtoms();
            for (unsigned int j = 0; j < all_atoms_in_assembly_residue.size(); j++){
                if (all_atoms_in_assembly_residue[j]->GetName() == parent_oxygen_name){
                    for (unsigned int k = 0; k < all_atoms_in_child_assembly_residue.size(); k++){
                        if (all_atoms_in_child_assembly_residue[k]->GetName() == anomeric_carbon_name){
                            Atom* parent_tail_atom = all_atoms_in_assembly_residue[j];
                            Atom* child_head_atom = all_atoms_in_child_assembly_residue[k];
                            parent_tail_child_head_map[parent_tail_atom] = child_head_atom;
                        }
                    }
                }
            }
        }

        //Also set a bond between child anomeric carbon and linking parent tail atom.
        for (std::map<Atom*, Atom*>::iterator tail_head_mapit = parent_tail_child_head_map.begin();
             tail_head_mapit != parent_tail_child_head_map.end(); tail_head_mapit++){
            Atom* parent_tail_atom = tail_head_mapit->first;
            Atom* child_head_atom = tail_head_mapit->second;
            //Add parent tail to child head's node neighbors.
            AtomNode* head_atom_node = NULL;
            if (child_head_atom->GetNode() == NULL){
                head_atom_node = new AtomNode();
                head_atom_node->SetAtom(child_head_atom);
                child_head_atom->SetNode(head_atom_node);
            }
            else
                head_atom_node = child_head_atom->GetNode();

            AtomVector existing_head_atom_neighbors = head_atom_node->GetNodeNeighbors();
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

            AtomVector existing_tail_atom_neighbors = tail_atom_node->GetNodeNeighbors();
            if (std::find (existing_tail_atom_neighbors.begin(), existing_tail_atom_neighbors.end(), child_head_atom) ==
                    existing_tail_atom_neighbors.end() ){
                tail_atom_node->AddNodeNeighbor(child_head_atom);
            }
        }

        //Add all tail/head atoms to corresponding residue's tail/head atoms, as well as corresponding residueNode connecting atoms
        for (std::map<Atom*, Atom*>::iterator it2 = parent_tail_child_head_map.begin(); it2 != parent_tail_child_head_map.end(); it2++){
            Atom* parent_tail_atom = it2->first;
            Atom* child_head_atom = it2->second;
            Residue* parent_residue = parent_tail_atom->GetResidue();
            Residue* child_residue = child_head_atom->GetResidue();
            parent_residue->AddTailAtom(parent_tail_atom);
            parent_residue->GetNode()->AddResidueNodeConnectingAtom(parent_tail_atom);
            child_residue->AddHeadAtom(child_head_atom);
            child_residue->GetNode()->AddResidueNodeConnectingAtom(child_head_atom);
        }
        //If child residue is derivative, adjust charge
        for (std::map<Atom*, Atom*>::iterator it2 = parent_tail_child_head_map.begin(); it2 != parent_tail_child_head_map.end(); it2++){
            Atom* parent_tail_atom = it2->first;
            Atom* child_head_atom = it2->second;
            Residue* parent_residue = parent_tail_atom->GetResidue();
            Residue* child_residue = child_head_atom->GetResidue();
            if (child_residue->GetIsSugarDerivative()){
                AtomVector all_parent_tail_atoms = parent_residue->GetTailAtoms();
                for (unsigned int i = 0; i < all_parent_tail_atoms.size(); i++){
                    if (all_parent_tail_atoms[i] == parent_tail_atom){
                        unsigned int branch_index = i;
                        this->AdjustCharge(child_residue, parent_residue, branch_index);
                    }
                }
            }
        }

    }//for
    //The for loop above set parent residue tail atoms and child head atoms. But if a residue is either at reducing or non-reducing terminal:
    //If a residue is the reducing terminal/aglycon, at this point it doesn't have head atoms set, since it has no parent. Set head atom equals tail atom.
    //If a residue is the non-reducing terminal, at this point it doesn't have tail atoms set, since it has no child. Set tail atom eqals head atom.
    for (std::map<CondensedSequenceSpace::CondensedSequenceGlycam06Residue*, Residue*>::iterator mapit =
         condensed_assembly_residue_map.begin(); mapit != condensed_assembly_residue_map.end(); mapit++){
        Residue* residue = mapit->second;
        AtomVector head_atoms = residue->GetHeadAtoms();
        AtomVector tail_atoms = residue->GetTailAtoms();
        //If this residue is reducing terminal, head atom is not set. Set below:
        if (head_atoms.empty() && !tail_atoms.empty()){
            residue->SetHeadAtoms(tail_atoms);
        }
        //If this residue is non-reducing terminal, tail atom is not set. Set below:
        if (!head_atoms.empty() && tail_atoms.empty()){
            residue->SetTailAtoms(head_atoms);
        }
        //If both not set, something is wrong. Print warning message:
        if (head_atoms.empty() && tail_atoms.empty()){
//            std::cout << "Warning: " << "Both head and tail atoms are not set(empty) for residue " << residue->GetName() << std::endl;
        }
    }

}//SetGlycam06ResidueBonding

void Assembly::RecursivelyTagDihedrals(Residue* parent_residue, std::multimap<int, std::pair<AtomVector*, std::string> >& index_dihedral_map, int& linkage_index)
{
    AtomVector all_tail_atoms = parent_residue->GetTailAtoms();
    AtomVector all_head_atoms = parent_residue->GetHeadAtoms();
    for (unsigned int i = 0; i < all_tail_atoms.size(); i++){
        Atom* tail_atom = all_tail_atoms[i];
        //For non-reducing terminal/derivative residues, tail atom equals head atom. In this case, the tail atoms aren't really connecting to a child residue.
        //If this is the case, stop recursion from further going down the oligosaccharide tree through such tail atoms.
        //This if statement below makes sure a tail atom is a head atom at the same time.If so, skip any operations.
        if ( parent_residue->GetIsAglycon() || std::find(all_head_atoms.begin(),all_head_atoms.end(),tail_atom) == all_head_atoms.end()){
            AtomVector tail_atom_neighbors = tail_atom->GetNode()->GetNodeNeighbors();
            AtomVector all_atoms_in_residue = parent_residue->GetAtoms();
            for (unsigned int j = 0; j < tail_atom_neighbors.size(); j++){
                Atom* neighbor_atom = tail_atom_neighbors[j];
                //if a neighbor is outside of a parent residue, it must be the head atom of a child residue.There should only exist one such atom, otherwise, something is amiss.
                if (std::find(all_atoms_in_residue.begin(), all_atoms_in_residue.end(), neighbor_atom) == all_atoms_in_residue.end()){
                    Atom* head_atom_of_child_residue = neighbor_atom;
                    Residue* child_residue = head_atom_of_child_residue->GetResidue();
                    AtomVector all_atoms_in_child_residue = child_residue->GetAtoms();

                    //Set C-O(tail atom)-C(head atom ) angle to 120 deg. The first C is the neighbor of tail atom that is not the head atom && not a hydrogen
                    //The only exception is ROH, where you have to use that hydrogen
                    //Right now,rely on the first letter of atom name to determine element type (if hydrogen or not). Better solution is the rule class.
                    Atom* non_hydrogen_tail_atom_neighbor = NULL;
                    for (unsigned int k = 0; k < tail_atom_neighbors.size(); k++){
                        if (parent_residue->GetName() == "ROH" && tail_atom_neighbors[k] != head_atom_of_child_residue)
                            non_hydrogen_tail_atom_neighbor = tail_atom_neighbors[k];

                        else if (tail_atom_neighbors[k] != head_atom_of_child_residue && tail_atom_neighbors[k]->GetName().substr(0,1) != "H"){
                            non_hydrogen_tail_atom_neighbor = tail_atom_neighbors[k];
                        }
                    }

                    //If child is a derivative, set only one dihedral: HX-CX-OX(tail atom)-head atom. Set this to 0 degree, making the derivative head eclipse the hydrogen. This is a general solution for bad derivative angles
                    //This torsion closely resembles the psi torsion for regular sugar-sugar connection.
                    if (child_residue->GetIsSugarDerivative()){
                        //! \todo Add these cout statements to the debugging mechanism once the DebugLevel class (or whatever) is implemented.
                        // std::cout << "Setting derivative psi torsion for residue: " << child_residue->GetName() <<std::endl;
                        Atom* derivative_atom_4 = head_atom_of_child_residue;
                        Atom* derivative_atom_3 = tail_atom;
                        Atom* derivative_atom_2 = non_hydrogen_tail_atom_neighbor;
                        Atom* derivative_atom_1 = NULL;
                        AtomVector derivative_atom_2_neighbors = derivative_atom_2->GetNode()->GetNodeNeighbors();
                        //Derivative atom 1 should be a exocyclic (normally hydrogen) neighbor of the neighbor of tail atom (neighbor of neighbor of tail oxygen), for example: H4-C4-O4-C1
                        for (unsigned int k = 0; k < derivative_atom_2_neighbors.size(); k++){
                            if (!derivative_atom_2_neighbors[k]->GetIsCycle() && derivative_atom_2_neighbors[k] != derivative_atom_3){
                                derivative_atom_1 = derivative_atom_2_neighbors[k];
                            }
                        }
                        if (derivative_atom_1 == NULL || derivative_atom_2 == NULL || derivative_atom_3 == NULL || derivative_atom_4 == NULL ){
                            //std::cout << "SetPsiDihedral: cannot find all four psi atoms. Skipping." << std::endl;
                        }
                        else {
                            //! \todo Add these cout statements to the debugging mechanism once the DebugLevel class (or whatever) is implemented.
                            // std::cout << "Derivative psi: " << derivative_atom_1->GetName() << "-" << derivative_atom_2->GetName() << "-" << derivative_atom_3->GetName() << "-" << derivative_atom_4->GetName() <<std::endl;
                            AtomVector* psi_atoms = new AtomVector();
                            psi_atoms->push_back(derivative_atom_1);
                            psi_atoms->push_back(derivative_atom_2);
                            psi_atoms->push_back(derivative_atom_3);
                            psi_atoms->push_back(derivative_atom_4);
                            std::pair<AtomVector*, std::string> dihedral_type_pair = std::make_pair(psi_atoms, "psi");
                            std::pair<int, std::pair<AtomVector*, std::string> > linkage_index_dihedral_type_pair = std::make_pair(linkage_index, dihedral_type_pair);
                            index_dihedral_map.insert(linkage_index_dihedral_type_pair);
                        }


                    }
                    else{
                        //Set Dihedral phi:
                        Atom* phi_atom_1 = non_hydrogen_tail_atom_neighbor;
                        Atom* phi_atom_2 = tail_atom;
                        Atom* phi_atom_3 = head_atom_of_child_residue;
                        Atom* phi_atom_4 = NULL;

                        AtomVector head_atom_neighbors = head_atom_of_child_residue->GetNode()->GetNodeNeighbors();
                        std::string anomeric_carbon_index_str = head_atom_of_child_residue->GetName().substr(1,1); //The "1" in C1, the "2" in C2
                        std::stringstream s1;
                        s1 << anomeric_carbon_index_str;
                        int anomeric_carbon_index;
                        s1 >> anomeric_carbon_index;
                        for (unsigned int k = 0; k< head_atom_neighbors.size(); k++){
                            Atom* neighbor = head_atom_neighbors[k];

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
                            //std::cout << "SetPhiDihedral: cannot find all four phi atoms. Skipping." << std::endl;
                        }
                        //If phi_atom_4 can be found, set phi.
                        else{
                            AtomVector* phi_atoms = new AtomVector();
                            phi_atoms->push_back(phi_atom_1);
                            phi_atoms->push_back(phi_atom_2);
                            phi_atoms->push_back(phi_atom_3);
                            phi_atoms->push_back(phi_atom_4);
                            std::pair<AtomVector*, std::string> dihedral_type_pair = std::make_pair(phi_atoms, "phi");
                            std::pair<int, std::pair<AtomVector*, std::string> > linkage_index_dihedral_type_pair = std::make_pair(linkage_index, dihedral_type_pair);
                            index_dihedral_map.insert(linkage_index_dihedral_type_pair);
                        }
                        //Set psi, for example: psi H4-C4-O4-C1 to 0 deg
                        Atom* psi_atom_4 = head_atom_of_child_residue;
                        Atom* psi_atom_3 = tail_atom;
                        Atom* psi_atom_2 = non_hydrogen_tail_atom_neighbor;
                        Atom* psi_atom_1 = NULL;
                        AtomVector psi_atom_2_neighbors = psi_atom_2->GetNode()->GetNodeNeighbors();
                        //Psi atom 1 should be a exocyclic (normally hydrogen) neighbor of the neighbor of tail atom (neighbor of neighbor of tail oxygen), for example: psi H4-C4-O4-C1
                        for (unsigned int k = 0; k < psi_atom_2_neighbors.size(); k++){
                            if (!psi_atom_2_neighbors[k]->GetIsCycle() && psi_atom_2_neighbors[k] != psi_atom_3){
                                psi_atom_1 = psi_atom_2_neighbors[k];
                            }
                        }

                        if (psi_atom_1 == NULL || psi_atom_2 == NULL || psi_atom_3 == NULL || psi_atom_4 == NULL ){
                            //std::cout << "SetPsiDihedral: cannot find all four psi atoms. Skipping." << std::endl;
                        }
                        else {
                            AtomVector* psi_atoms = new AtomVector();
                            psi_atoms->push_back(psi_atom_1);
                            psi_atoms->push_back(psi_atom_2);
                            psi_atoms->push_back(psi_atom_3);
                            psi_atoms->push_back(psi_atom_4);
                            std::pair<AtomVector*, std::string> dihedral_type_pair = std::make_pair(psi_atoms, "psi");
                            std::pair<int, std::pair<AtomVector*, std::string> > linkage_index_dihedral_type_pair = std::make_pair(linkage_index, dihedral_type_pair);
                            index_dihedral_map.insert(linkage_index_dihedral_type_pair);
                        }

                        //Set Omega (if exists) , for example, C4-C5-C6-O6. Set this to 180 deg
                        //Set Omega (if exists) , for example, C4-C5-C6-O6. Set this to 180 deg

                        Atom* omega_atom_4 = psi_atom_3;
                        Atom* omega_atom_3 = psi_atom_2;
                        Atom* omega_atom_2 = NULL;
                        Atom* omega_atom_1 = NULL;
                        //omega atom 3 should be a exocyclic non-hydrogen(probably carbon) atom that 's connects to the atoms that connects to tail atom (neighbor of neighbor of tail oxygen)
                        //If such an atom is already on the ring, then there is no omega angle. If it is exocyclic, then omega exists.
                        if (omega_atom_3->GetIsCycle()){
                            //std::cout << "Tail atom is directly attached to ring atom. So omega does not exist. Skip setting omega." << std::endl;
                        }
                        else{
                            //Choose omega atom 2 from the neighbors of omega atom 3. It can't be omega_atom_4, and it shouldn't be a hydrogen
                            AtomVector omega_atom_3_neighbors = omega_atom_3->GetNode()->GetNodeNeighbors();
                            for (AtomVector::iterator atom_it4 = omega_atom_3_neighbors.begin(); atom_it4 != omega_atom_3_neighbors.end(); atom_it4++){
                                Atom* neighbor = *atom_it4;
                                if (neighbor->GetElementSymbol() != "H" && neighbor != omega_atom_4){
                                    omega_atom_2 = neighbor;
                                }
                            }
                        }
                        //Once omega atom 2 is identified, get its non-hydrogen node neighbor, this should be omega atom 1.
                        if (omega_atom_2 !=NULL){       //if there is such an exocyclic atom, then omega atom 2 exists.
                            AtomVector omega_atom_2_neighbors = omega_atom_2->GetNode()->GetNodeNeighbors();
                            for (AtomVector::iterator atom_it5 = omega_atom_2_neighbors.begin(); atom_it5 != omega_atom_2_neighbors.end(); atom_it5++){
                                Atom* neighbor = *atom_it5;
                                if (neighbor->GetElementSymbol() != "H" && neighbor != omega_atom_3){
                                    omega_atom_1 = neighbor;
                                }
                            }
                        }
                        //If all four omega atoms exist, set omega to -60 deg, assuming gt.
                        if (omega_atom_4 == NULL || omega_atom_3 == NULL || omega_atom_2 == NULL || omega_atom_1 == NULL){
                            //std::cout << "SetOmegaDihedral: cannot find all four omega atoms. Skipping." << std::endl;
                        }
                        else {
                            AtomVector* omega_atoms = new AtomVector();
                            omega_atoms->push_back(omega_atom_1);
                            omega_atoms->push_back(omega_atom_2);
                            omega_atoms->push_back(omega_atom_3);
                            omega_atoms->push_back(omega_atom_4);
                            std::pair<AtomVector*, std::string> dihedral_type_pair = std::make_pair(omega_atoms, "omega");
                            std::pair<int, std::pair<AtomVector*, std::string> > linkage_index_dihedral_type_pair = std::make_pair(linkage_index, dihedral_type_pair);
                            index_dihedral_map.insert(linkage_index_dihedral_type_pair);
                        }

                    }//else Done setting phi,psi, omega(if exists)
                    //Start new recursion
                    Residue* new_parent_residue = child_residue;
                    linkage_index++;
                    this->RecursivelyTagDihedrals(new_parent_residue, index_dihedral_map, linkage_index);
                }//if
            }//for
        }//if

    }//for

}

void Assembly::RecursivelySetAngleGeometry (Residue* parent_residue)
{
    AtomVector all_tail_atoms = parent_residue->GetTailAtoms();
    AtomVector all_head_atoms = parent_residue->GetHeadAtoms();
    for (unsigned int i = 0; i < all_tail_atoms.size(); i++){
        Atom* tail_atom = all_tail_atoms[i];
        //For non-reducing terminal/derivative residues, tail atom equals head atom. In this case, the tail atoms aren't really connecting to a child residue.
        //If this is the case, stop recursion from further going down the oligosaccharide tree through such tail atoms.
        //This if statement below makes sure a tail atom is a head atom at the same time.If so, skip any operations.
        if ( parent_residue->GetIsAglycon() || std::find(all_head_atoms.begin(),all_head_atoms.end(),tail_atom) == all_head_atoms.end()){
            AtomVector tail_atom_neighbors = tail_atom->GetNode()->GetNodeNeighbors();
            AtomVector all_atoms_in_residue = parent_residue->GetAtoms();
            for (unsigned int j = 0; j < tail_atom_neighbors.size(); j++){
                Atom* neighbor_atom = tail_atom_neighbors[j];
                //if a neighbor is outside of a parent residue, it must be the head atom of a child residue.There should only exist one such atom, otherwise, something is amiss.
                if (std::find(all_atoms_in_residue.begin(), all_atoms_in_residue.end(), neighbor_atom) == all_atoms_in_residue.end()){
                    Atom* head_atom_of_child_residue = neighbor_atom;
                    Residue* child_residue = head_atom_of_child_residue->GetResidue();
                    //Right now, all residues are at the position of the template residue. That is, they are all around the orgin and stacked upon each other.
                    //SetResidueResidueBondDistance function: takes a pair of parent tail/child head atoms as argument. This function keeps the parent residue intact,but
                    //finds out the new position of child head atom, and move atoms of child residue accordingly.(i.e. grafting)
                    this->SetResidueResidueBondDistance(tail_atom, head_atom_of_child_residue);

                    //Set C-O(tail atom)-C(head atom ) angle to 120 deg. The first C is the neighbor of tail atom that is not the head atom && not a hydrogen
                    //The only exception is ROH, where you have to use that hydrogen
                    //Right now,rely on the first letter of atom name to determine element type (if hydrogen or not). Better solution is the rule class.
                    Atom* non_hydrogen_tail_atom_neighbor = NULL;
                    for (unsigned int k = 0; k < tail_atom_neighbors.size(); k++){
                        if (parent_residue->GetName() == "ROH" && tail_atom_neighbors[k] != head_atom_of_child_residue)
                            non_hydrogen_tail_atom_neighbor = tail_atom_neighbors[k];

                        else if (tail_atom_neighbors[k] != head_atom_of_child_residue && tail_atom_neighbors[k]->GetName().substr(0,1) != "H"){
                            non_hydrogen_tail_atom_neighbor = tail_atom_neighbors[k];
                        }
                    }

                    //Exit if atom 1 cannot be be found
                    if (non_hydrogen_tail_atom_neighbor != NULL){
                        const double angle_to_set = 109.4;	//assuming sp3 tetrahedral
                        this->SetAngle(non_hydrogen_tail_atom_neighbor, tail_atom, head_atom_of_child_residue, angle_to_set);
                    }
                    //Start new recursion
                    Residue* new_parent_residue = child_residue;
                    this->RecursivelySetAngleGeometry(new_parent_residue);
                }
            }
        }
    }
}

// OG These two functions are moving to carbohydrateBuilder
void Assembly::FigureOutResidueLinkagesInGlycan(Residue *from_this_residue1, Residue *to_this_residue2, ResidueLinkageVector *residue_linkages)
{
    ResidueVector neighbors = to_this_residue2->GetNode()->GetResidueNeighbors();
    for(auto &neighbor : neighbors)
    {
        if(neighbor->GetIndex() != from_this_residue1->GetIndex()) // If not the previous residue
        {
            residue_linkages->emplace_back(neighbor, to_this_residue2);
        }
    }
    for(auto &neighbor : neighbors)
    {
        if(neighbor->GetIndex() != from_this_residue1->GetIndex())
        {
            this->FigureOutResidueLinkagesInGlycan(to_this_residue2, neighbor, residue_linkages);
        }
    }
    return;
}

// OG. This is agnostic of it being generated from a condensed sequence. Could be any assembly that is a glycan with good connectivites and GLYCAM nomenclature for atoms and residues.
void Assembly::SetDihedralAngleGeometryWithMetadata()
{
    ResidueLinkageVector all_residue_linkages;
    this->FigureOutResidueLinkagesInGlycan(this->GetResidues().at(0), this->GetResidues().at(0), &all_residue_linkages);
    for(auto &linkage : all_residue_linkages)
    {
        linkage.SetDefaultShapeUsingMetadata();

    }
    // Resovlving overlaps should be a separated function, but I don't want assembly to have a ResidueLinkageVector member. Need new, seperate class.
    for(auto &linkage : all_residue_linkages)
    {
        //AtomVector overlapAtomSet1, AtomVector overlapAtomSet2, double overlapTolerance, int angleIncrement
        linkage.SimpleWiggle(this->GetAllAtomsOfAssembly(), this->GetAllAtomsOfAssembly(), 0.1, 5);
    }
}

/* Oliver needs to clarify what he is doing here:
 * I don't currently have the ability to "select" which rotamers I will write out. I need to add default and overwritable indexes for residue_linkage class for that to work.
 * I want this "selection" to work at the level of AddMetadata. So I can only add relevant metadata while linkages or rotamers that are deselected are not used.
 *  Having flags active and inactive metadata in that class might be another option.
 * I've updated the pdb writer to allow me to pass in a model_index number X, which will trigger writing out of coordinate set X from coordinateVector in Atom.
 * I can create an Assembly::saveCurrentCoordinates function that uses Atom::AddCoordinate to push_back a coordinate set related to a rotamer.
 * Then I'll create a wrapper function that allows me to write out each coordinate set.... Either I can pass in a name (e.g. residue_linkage_index_ gg) or
 *  I just write the PDB file as I'm creating all the possible combos... But same issue with naming.
 * The following functions have been added by me already:
 * Assembly::FigureOutResidueLinkagesInGlycan
 * Assembly::SetDihedralAngleGeometryWithMetadata
 * In Residue_Linkage I have a simple wiggle function, I have implemented it here. The problem is that it explores all rotamers. That won't work.
 * Functionality is there for it to only explore within one rotamer, but need higher level control of that.


  */


// By Yao, Oliver wishes to replace  with SetDihedralAngleGeometryWithMetadata
void Assembly::RecursivelySetDihedralAngleGeometry (Residue* parent_residue)
{
    AtomVector all_tail_atoms = parent_residue->GetTailAtoms();
    AtomVector all_head_atoms = parent_residue->GetHeadAtoms();
    for (unsigned int i = 0; i < all_tail_atoms.size(); i++){
        Atom* tail_atom = all_tail_atoms[i];
        //For non-reducing terminal/derivative residues, tail atom equals head atom. In this case, the tail atoms aren't really connecting to a child residue.
        //If this is the case, stop recursion from further going down the oligosaccharide tree through such tail atoms.
        //This if statement below makes sure a tail atom is a head atom at the same time.If so, skip any operations.
        if ( parent_residue->GetIsAglycon() || std::find(all_head_atoms.begin(),all_head_atoms.end(),tail_atom) == all_head_atoms.end()){
            AtomVector tail_atom_neighbors = tail_atom->GetNode()->GetNodeNeighbors();
            AtomVector all_atoms_in_residue = parent_residue->GetAtoms();
            for (unsigned int j = 0; j < tail_atom_neighbors.size(); j++){
                Atom* neighbor_atom = tail_atom_neighbors[j];
                //if a neighbor is outside of a parent residue, it must be the head atom of a child residue.There should only exist one such atom, otherwise, something is amiss.
                if (std::find(all_atoms_in_residue.begin(), all_atoms_in_residue.end(), neighbor_atom) == all_atoms_in_residue.end()){
                    Atom* head_atom_of_child_residue = neighbor_atom;
                    Residue* child_residue = head_atom_of_child_residue->GetResidue();
                    //If child is a derivative, set only one dihedral: HX-CX-OX(tail atom)-head atom. Set this to 0 degree, making the derivative head eclipse the hydrogen. This is a general solution for bad derivative angles
                    //This torsion closely resembles the psi torsion for regular sugar-sugar connection.
                    Atom* non_hydrogen_tail_atom_neighbor = NULL;
                    for (unsigned int k = 0; k < tail_atom_neighbors.size(); k++){
                        if (parent_residue->GetName() == "ROH" && tail_atom_neighbors[k] != head_atom_of_child_residue)
                            non_hydrogen_tail_atom_neighbor = tail_atom_neighbors[k];

                        else if (tail_atom_neighbors[k] != head_atom_of_child_residue && tail_atom_neighbors[k]->GetName().substr(0,1) != "H"){
                            non_hydrogen_tail_atom_neighbor = tail_atom_neighbors[k];
                        }
                    }
                    if (child_residue->GetIsSugarDerivative()){
                        Atom* derivative_atom_4 = head_atom_of_child_residue;
                        Atom* derivative_atom_3 = tail_atom;
                        Atom* derivative_atom_2 = non_hydrogen_tail_atom_neighbor;
                        Atom* derivative_atom_1 = NULL;
                        AtomVector derivative_atom_2_neighbors = derivative_atom_2->GetNode()->GetNodeNeighbors();
                        //Derivative atom 1 should be a exocyclic (normally hydrogen) neighbor of the neighbor of tail atom (neighbor of neighbor of tail oxygen), for example: H4-C4-O4-C1
                        for (unsigned int k = 0; k < derivative_atom_2_neighbors.size(); k++){
                            if (!derivative_atom_2_neighbors[k]->GetIsCycle() && derivative_atom_2_neighbors[k] != derivative_atom_3){
                                derivative_atom_1 = derivative_atom_2_neighbors[k];
                            }
                        }
                        if (derivative_atom_1 != NULL && derivative_atom_2 != NULL && derivative_atom_3 != NULL && derivative_atom_4 != NULL ){
                            const double derivative_psi = 0.0;
                            this->SetDihedral(derivative_atom_1, derivative_atom_2, derivative_atom_3, derivative_atom_4, derivative_psi);
                        }
                    }
                    //Otherwise,attempt setting dihedrals
                    else
                    {
                        //Set Dihedral phi:
                        Atom* phi_atom_1 = non_hydrogen_tail_atom_neighbor;
                        Atom* phi_atom_2 = tail_atom;
                        Atom* phi_atom_3 = head_atom_of_child_residue;
                        Atom* phi_atom_4 = NULL;

                        AtomVector head_atom_neighbors = head_atom_of_child_residue->GetNode()->GetNodeNeighbors();
                        std::string anomeric_carbon_index_str = head_atom_of_child_residue->GetName().substr(1,1); //The "1" in C1, the "2" in C2
                        std::stringstream s1;
                        s1 << anomeric_carbon_index_str;
                        int anomeric_carbon_index;
                        s1 >> anomeric_carbon_index;
                        for (unsigned int k = 0; k< head_atom_neighbors.size(); k++){
                            Atom* neighbor = head_atom_neighbors[k];

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
                        if (phi_atom_4 != NULL && phi_atom_3 != NULL && phi_atom_2 != NULL && phi_atom_1 != NULL ){
                            const double dihedral_phi = 180.0;
                            this->SetDihedral(phi_atom_1, phi_atom_2, phi_atom_3, phi_atom_4, dihedral_phi);
                        }
                        //Set psi, for example: psi H4-C4-O4-C1 to 0 deg
                        Atom* psi_atom_4 = head_atom_of_child_residue;
                        Atom* psi_atom_3 = tail_atom;
                        Atom* psi_atom_2 = non_hydrogen_tail_atom_neighbor;
                        Atom* psi_atom_1 = NULL;
                        AtomVector psi_atom_2_neighbors = psi_atom_2->GetNode()->GetNodeNeighbors();
                        //Psi atom 1 should be a exocyclic (normally hydrogen) neighbor of the neighbor of tail atom (neighbor of neighbor of tail oxygen), for example: psi H4-C4-O4-C1
                        for (unsigned int k = 0; k < psi_atom_2_neighbors.size(); k++){
                            if (!psi_atom_2_neighbors[k]->GetIsCycle() && psi_atom_2_neighbors[k] != psi_atom_3){
                                psi_atom_1 = psi_atom_2_neighbors[k];
                            }
                        }

                        if (psi_atom_1 != NULL && psi_atom_2 != NULL && psi_atom_3 != NULL && psi_atom_4 != NULL )
                        {
                            // OG 2019.09.05 quickfix:
                            // The psi angle is ok for 2 bond linkages (e.g. 1-2, 1-4 etc) but for 3 bond and longer we use
                            // different atoms and different psi values. First I get the correct atom, and then I set dihedral to 180 instead of 0.
                            double dihedral_psi = 0.0; //default value for most linkages
                            // Check if this is a 3-bond or longer linkage:
                            std::regex regex_query("O[6-9]", std::regex_constants::ECMAScript);
                            if (std::regex_search(tail_atom->GetName(), regex_query) ) // if tail atom name matches "O[6-9]"
                            {
                                dihedral_psi = 180.0;
                                //Psi atom 1 should be a non-hydrogen neighbor of the neighbor of tail atom (neighbor of neighbor of tail oxygen), for example: psi C5-C6-O6-C1
                                for (unsigned int k = 0; k < psi_atom_2_neighbors.size(); k++)
                                {   // If not a hydrogen and not psi_atom_3
                                    if((psi_atom_2_neighbors[k]->GetName().substr(0,1) != "H") && (psi_atom_2_neighbors[k] != psi_atom_3))
                                    {
                                        psi_atom_1 = psi_atom_2_neighbors[k];
                                    }
                                }
                            }
                            // End OG 2019.09.05 quickfix
                            this->SetDihedral(psi_atom_1, psi_atom_2, psi_atom_3, psi_atom_4, dihedral_psi);
                        }

                        //Set Omega (if exists) , for example, C4-C5-C6-O6. Set this to 180 deg

                        Atom* omega_atom_4 = psi_atom_3;
                        Atom* omega_atom_3 = psi_atom_2;
                        Atom* omega_atom_2 = NULL;
                        Atom* omega_atom_1 = NULL;
                        //omega atom 3 should be a exocyclic non-hydrogen(probably carbon) atom that 's connects to the atoms that connects to tail atom (neighbor of neighbor of tail oxygen)
                        //If such an atom is already on the ring, then there is no omega angle. If it is exocyclic, then omega exists.
                        if (omega_atom_3->GetIsCycle()){
                            //std::cout << "Tail atom is directly attached to ring atom. So omega does not exist. Skip setting omega." << std::endl;
                        }
                        else{
                            //Choose omega atom 2 from the neighbors of omega atom 3. It can't be omega_atom_4, and it shouldn't be a hydrogen
                            AtomVector omega_atom_3_neighbors = omega_atom_3->GetNode()->GetNodeNeighbors();
                            for (AtomVector::iterator atom_it4 = omega_atom_3_neighbors.begin(); atom_it4 != omega_atom_3_neighbors.end(); atom_it4++){
                                Atom* neighbor = *atom_it4;
                                if (neighbor->GetElementSymbol() != "H" && neighbor != omega_atom_4){
                                    omega_atom_2 = neighbor;
                                }
                            }
                        }
                        //Once omega atom 2 is identified, get its non-hydrogen node neighbor, this should be omega atom 1.
                        if (omega_atom_2 !=NULL){	//if there is such an exocyclic atom, then omega atom 2 exists.
                            AtomVector omega_atom_2_neighbors = omega_atom_2->GetNode()->GetNodeNeighbors();
                            for (AtomVector::iterator atom_it5 = omega_atom_2_neighbors.begin(); atom_it5 != omega_atom_2_neighbors.end(); atom_it5++){
                                Atom* neighbor = *atom_it5;
                                if (neighbor->GetElementSymbol() != "H" && neighbor != omega_atom_3){
                                    omega_atom_1 = neighbor;
                                }
                            }
                        }
                        if ( omega_atom_4 != NULL && omega_atom_3 != NULL && omega_atom_2 != NULL && omega_atom_1 != NULL ) {
                            const double dihedral_omega = -60.0;
                            this->SetDihedral(omega_atom_1, omega_atom_2, omega_atom_3, omega_atom_4, dihedral_omega);
                        }
                    }//else Done setting phi,psi, omega(if exists)

                    //Start new recursion
                    Residue* new_parent_residue = child_residue;
                    this->RecursivelySetDihedralAngleGeometry(new_parent_residue);
                }//if
            }//for
        }//if

    }//for
}

ResidueVector Assembly::FindClashingResidues()
{
    //TODO: to make the code run faster, for the nested for loops that iterate througha atoms, it1 = it+1
    ResidueVector clashing_residues;
    AtomVector all_atoms_of_assembly = this->GetAllAtomsOfAssembly();
    //Exhaustively compare two different atoms in assembly, using two nested for loops
    for (AtomVector::iterator it = all_atoms_of_assembly.begin(); it != all_atoms_of_assembly.end() ; it++){
        Atom* current_atom = *it;
        unsigned int bond_by_distance_count = 0;  //How many bonds does a particular atom have, according to bond by distance
        for (AtomVector::iterator it1 = it + 1; it1 != all_atoms_of_assembly.end(); it1++)
        {
            Atom* another_atom = *it1;
            if(current_atom->CheckIfOtherAtomIsWithinBondingDistance(another_atom))
            {
                ++bond_by_distance_count;
            }
        }
        unsigned int actual_bond_count = current_atom->GetNode()->GetNodeNeighbors().size();
        //If bond by distance gives more bonds than the number of bonds in atom node, then this atom is clashing with other atoms. This residue is a clashing residue
        if (bond_by_distance_count > actual_bond_count){
            //If multiple atoms in a residue is clashing, prevent duplicate tagging of the corresponding residue as clashing
            if (std::find (clashing_residues.begin(), clashing_residues.end(), current_atom->GetResidue()) == clashing_residues.end())
            {
                clashing_residues.push_back(current_atom->GetResidue());
            }
        }
    }
    //! \todo Add these cout statements to the debugging mechanism once the DebugLevel class (or whatever) is implemented.
    // std::cout << "All clashing residues were identified in this order:" << std::endl;
    // for (ResidueVector::iterator it = clashing_residues.begin(); it != clashing_residues.end(); it++){
    // std::cout << (*it)->GetIndex() << "---" << (*it)->GetName() << std::endl;
    // }
    // std::cout << "\n";
    return clashing_residues;
}

std::vector<ResidueVector> Assembly::FindPathToCommonAncestors(ResidueVector& all_clashing_residues)
{
    std::vector<ResidueVector> all_clashing_residue_parent_paths = std::vector<ResidueVector> ();
    ResidueVector visited_residues = ResidueVector();
    //Starting from each residue, construct a pathway until a branching point(a residue with multiple tail atoms), add residue in this pathway to a ResidueVector
    for (ResidueVector::iterator it = all_clashing_residues.begin(); it != all_clashing_residues.end(); it++){
        Residue* clashing_residue = *it;
        if (!clashing_residue->GetIsAglycon() && std::find(visited_residues.begin(), visited_residues.end(), clashing_residue) == visited_residues.end() ){
            Residue* current_residue = clashing_residue;
            ResidueVector path = ResidueVector();
            //Going from the non-reducing end(starting from clashing residue) towards the reducing end using head atom -> parent tail atom relationship
            while (true){
                //If the code has reached the aglycon, or has reached a residue that's already visited. Then this residue is the endpoint of a clashing pathway.
                //It's either the aglycon, or a residue with multiple branches. This is the last element is a clashing pathway residue vector, and is called
                //"common ancestor" in another function somewhere downstream. Termination condition is reached, break out of while loop
                if (current_residue->GetIsAglycon() || std::find(visited_residues.begin(), visited_residues.end(), current_residue) != visited_residues.end()){
                    visited_residues.push_back(current_residue);
                    path.push_back(current_residue);
                    all_clashing_residue_parent_paths.push_back(path);
                    break;
                }
                //Otherwise, add the current residue to clashing pathway and keep going upward.
                else{
                    visited_residues.push_back(current_residue);
                    path.push_back(current_residue);
                    AtomVector head_atom_neighbors = current_residue->GetHeadAtoms().at(0)->GetNode()->GetNodeNeighbors();
                    for (unsigned int i = 0; i < head_atom_neighbors.size(); i++){
                        Atom* neighbor = head_atom_neighbors[i];
                        if (neighbor->GetResidue() != current_residue){
                            Residue* parent_residue = neighbor->GetResidue();
                            current_residue = parent_residue;
                        }
                    }
                }
            }

        }
    }
    return all_clashing_residue_parent_paths;
}

void Assembly::ResolveClashes(std::vector<ResidueVector>& fused_clashing_paths,
                              std::multimap<int, std::pair<AtomVector*, std::string> >& index_dihedral_map)
{
    //For each common ancestor, go through all its clashing pathways one by one.
    for (std::vector<ResidueVector>::iterator it = fused_clashing_paths.begin(); it != fused_clashing_paths.end(); it++){
        std::vector <AtomVector*> all_omega_dihedrals = this->FindAllOmegaTorsionsInPathway(*it, index_dihedral_map);
        //If availble dihedrals are found, initiate clash resolution process
        if (!all_omega_dihedrals.empty()){
            //By limited grid search, find the set of coordinate resulting in least clash
            GeometryTopology::CoordinateVector least_clash_coordinates_for_this_pathway = this->FindBestSetOfTorsions(all_omega_dihedrals);
            //For each atom in assembly, set coordinate according to the best set of coordiante found. This will crudely resolve clashes.
            AtomVector all_atoms_in_assembly = this->GetAllAtomsOfAssembly();
            for (unsigned int j = 0; j < all_atoms_in_assembly.size(); j++){
                GeometryTopology::Coordinate* new_coordinate = least_clash_coordinates_for_this_pathway[j];
                GeometryTopology::CoordinateVector new_coordinate_set = GeometryTopology::CoordinateVector();
                new_coordinate_set.push_back(new_coordinate);
                all_atoms_in_assembly[j]->SetCoordinates(new_coordinate_set);
            }
        }
    }
}

std::vector< AtomVector* > Assembly::FindAllOmegaTorsionsInPathway (ResidueVector& pathway, std::multimap<int, std::pair<AtomVector*, std::string> >&
                                                                    index_dihedral_map)
{
    //Identify all head atoms present in pathway
    AtomVector all_head_atoms_in_pathway = AtomVector();
    for (ResidueVector::reverse_iterator it = pathway.rbegin(); it != pathway.rend(); it++){
        AtomVector head_atoms_in_residue = (*it)->GetHeadAtoms();
        for (AtomVector::iterator it2 = head_atoms_in_residue.begin(); it2 != head_atoms_in_residue.end(); it2++){
            all_head_atoms_in_pathway.push_back(*it2);
        }
    }
    //Identify all tail atoms connected to the head atoms in pathway.
    AtomVector all_tail_atoms_in_pathway = AtomVector();
    for (AtomVector::iterator it = all_head_atoms_in_pathway.begin(); it != all_head_atoms_in_pathway.end(); it++){
        Atom* head_atom = *it;
        AtomVector head_atom_neighbors = head_atom->GetNode()->GetNodeNeighbors();
        for (AtomVector::iterator it2 = head_atom_neighbors.begin(); it2 != head_atom_neighbors.end(); it2++){
            Atom* neighbor = *it2;
            if (neighbor->GetResidue() != head_atom->GetResidue()){
                Atom* connected_tail_atom = neighbor;
                all_tail_atoms_in_pathway.push_back(connected_tail_atom);
            }
        }
    }
    //Phi,psi, or omega torsion angles all contain the tail atom. So, look at each tail atom:
    std::vector<AtomVector*>  omega_torsions_in_pathway = std::vector<AtomVector*>();
    for (AtomVector::iterator it = all_tail_atoms_in_pathway.begin(); it != all_tail_atoms_in_pathway.end(); it++){
        Atom* tail_atom = *it;
        //From all available dihedrals:
        for (std::multimap<int, std::pair<AtomVector*, std::string> >::iterator it2 = index_dihedral_map.begin(); it2 != index_dihedral_map.end(); it2++){
            AtomVector* dihedral_atoms = it2->second.first;
            std::string dihedral_type = it2->second.second;
            //If a dihedral contains a tail atom in the pathway, and has type "psi" or "omega", this is the dihedral we want. Later perform grid search on these torsion to alleviate clashes.
            if (std::find(dihedral_atoms->begin(), dihedral_atoms->end(), tail_atom) != dihedral_atoms->end() && (dihedral_type == "psi" || dihedral_type == "omega")){
                omega_torsions_in_pathway.push_back(dihedral_atoms);
            }
        }
    }
    return omega_torsions_in_pathway;
}

GeometryTopology::CoordinateVector Assembly::FindBestSetOfTorsions(std::vector<AtomVector*>& available_dihedrals)
{
    //For each dihedral, enable them to rotate -5.0, 0, and 5.0 degrees.
    std::vector<double> rotation_values = std::vector<double>();
    rotation_values.push_back(0.0);
    rotation_values.push_back(-5.0);
    rotation_values.push_back(5.0);
    //Start generating a relationship betwewn each dihedral and all rotation values, in the form of a vector of pairs. Each pair represents one dihedral
    //Pair.first is an AtomVector* containing the four atoms of a dihedral. Pair.second is a vector of double. In this function, it's always "rotation_values"
    std::vector<std::pair<AtomVector*, std::vector<double> > > all_dihedral_rotation_values = std::vector<std::pair<AtomVector*, std::vector<double> > >();
    for (unsigned int i = 0; i < available_dihedrals.size(); i++){
        AtomVector* dihedral_atoms = available_dihedrals[i];
        std::vector<double> possible_rotation_values = rotation_values;
        all_dihedral_rotation_values.push_back(std::make_pair(dihedral_atoms, possible_rotation_values));
    }
    //For the vector of pair generated above, go through each dihedral, expand their rotation angles in the form of recursive for loops. Exhaustively list each possibility as a "combination"
    typedef std::vector<std::pair<AtomVector*, double> > combination;
    std::vector<combination> all_combinations = std::vector<combination>();
    std::vector<double> angle_index_per_dihedral = std::vector<double> (available_dihedrals.size(), gmml::dNotSet);
    this ->GenerateAllTorsionCombinations(all_dihedral_rotation_values, 0, all_combinations, angle_index_per_dihedral);
    //For each combination, rotate coordinates accordingly, then compute clash score. If clash score is lower than the current minimum, record the current set of coordinate as the least-clashing
    //coordinate set
    //Make a pair. Pair.first is clash score, pair.second is coordinate set
    std::pair <double, GeometryTopology::CoordinateVector> least_clash_coordinate = std::pair <double, GeometryTopology::CoordinateVector>();
    AtomVector all_atoms_in_assembly = this->GetAllAtomsOfAssembly();
    //Compute rotation score using the originial coordinate. Initiate least_clash_coordinate to contain this original state.
    double initial_clash_score = gmml::CalculateAtomicOverlaps(all_atoms_in_assembly, all_atoms_in_assembly);
    least_clash_coordinate.first = initial_clash_score;
    for (unsigned int i = 0; i < all_atoms_in_assembly.size(); i++){
        Atom* atom = all_atoms_in_assembly[i];
        least_clash_coordinate.second.push_back( new GeometryTopology::Coordinate(atom->GetCoordinates().at(0)));
    }
    //For each combination, rotate accordingly
    for (unsigned int i = 0; i < all_combinations.size(); i++){
        combination& rotation_set = all_combinations[i];
        for (combination::iterator it = rotation_set.begin(); it != rotation_set.end(); it++){
            AtomVector* dihedral_atoms = it->first;
            GeometryTopology::Coordinate* pivot_point = dihedral_atoms->at(1)->GetCoordinates().at(0);
            GeometryTopology::Coordinate* rotation_axis = new GeometryTopology::Coordinate();
            rotation_axis->SetX(dihedral_atoms->at(2)->GetCoordinates().at(0)->GetX() - dihedral_atoms->at(1)->GetCoordinates().at(0)->GetX());
            rotation_axis->SetY(dihedral_atoms->at(2)->GetCoordinates().at(0)->GetY() - dihedral_atoms->at(1)->GetCoordinates().at(0)->GetY());
            rotation_axis->SetZ(dihedral_atoms->at(2)->GetCoordinates().at(0)->GetZ() - dihedral_atoms->at(1)->GetCoordinates().at(0)->GetZ());
            AtomVector atoms_to_rotate = AtomVector();
            atoms_to_rotate.push_back(dihedral_atoms->at(1));
            dihedral_atoms->at(2)->FindConnectedAtoms(atoms_to_rotate);
            GeometryTopology::CoordinateVector original_coordinates = GeometryTopology::CoordinateVector();
            for (unsigned int j = 0; j < atoms_to_rotate.size(); j++){
                original_coordinates.push_back(atoms_to_rotate[j]->GetCoordinates().at(0));
            }
            double rotation_value = it->second;
            GeometryTopology::Rotation rotation_operation = GeometryTopology::Rotation();
            //Perform rotation, but return rotated coordinates rather than actually perform rotation
            GeometryTopology::CoordinateVector rotated_coordinates = rotation_operation.RotateCoordinates(pivot_point, rotation_axis, rotation_value, original_coordinates);
            //Reset atom coordinates to the rotated coordinated set
            for (unsigned int j = 0; j < atoms_to_rotate.size(); j++){
                GeometryTopology::CoordinateVector new_coordinates = GeometryTopology::CoordinateVector();
                new_coordinates.push_back(rotated_coordinates[j]);
                atoms_to_rotate[j]->SetCoordinates(new_coordinates);

            }

        }//rotate each combination
        //After setting coordinates, compute clash score.
        double clash_score = gmml::CalculateAtomicOverlaps(all_atoms_in_assembly, all_atoms_in_assembly);
        //If current clash score is lower than current minimum, overwrite lowest clash score, and best coordinate set.
        if (clash_score < least_clash_coordinate.first){
            least_clash_coordinate.first = clash_score;
            for (unsigned int j = 0; j < all_atoms_in_assembly.size(); j++){
                least_clash_coordinate.second[j]->SetX(all_atoms_in_assembly[j]->GetCoordinates().at(0)->GetX());
                least_clash_coordinate.second[j]->SetY(all_atoms_in_assembly[j]->GetCoordinates().at(0)->GetY());
                least_clash_coordinate.second[j]->SetZ(all_atoms_in_assembly[j]->GetCoordinates().at(0)->GetZ());
            }
        }

    }
    //Return best set of coordinate
    return least_clash_coordinate.second;
}

void Assembly::GenerateAllTorsionCombinations(std::vector<std::pair<AtomVector*, std::vector<double> > >& all_dihedral_rotation_values, unsigned int current_dihedral_index ,
                                              std::vector<std::vector<std::pair<AtomVector*, double> > >& container_for_combinations, std::vector<double>& angle_index_per_dihedral)
{
    //For each dihedral to populate, get all its allowed rotation values
    std::vector<double>& rotation_values = all_dihedral_rotation_values[current_dihedral_index].second;
    //For each rotation value
    for (unsigned int j = 0; j < rotation_values.size(); j++){
        //Record the current rotation value of the current dihedral to the temporary container "angle_index_per_dihedral"
        angle_index_per_dihedral[current_dihedral_index] = rotation_values[j];
        //If the current dihedral is not the last one to rotate, move the recursive process to the next dihedral
        if (current_dihedral_index != all_dihedral_rotation_values.size() -1){
            unsigned int next_dihedral_index = current_dihedral_index + 1;
            this->GenerateAllTorsionCombinations(all_dihedral_rotation_values, next_dihedral_index, container_for_combinations, angle_index_per_dihedral);
        }
        //If the current dihedral is already the last dihedral, the code has reached the "dead end" of nested for loop.Now look at the temporary container for the rotation value at each
        //dihedral. These become a new combination.
        else{
            std::vector<std::pair<AtomVector*, double> > new_combination = std::vector<std::pair<AtomVector*, double> >();
            for (unsigned int dihedral_position = 0; dihedral_position < angle_index_per_dihedral.size(); dihedral_position++){
                AtomVector* dihedral = all_dihedral_rotation_values[dihedral_position].first;
                double torsion_value_at_this_dihedral = angle_index_per_dihedral[dihedral_position];
                new_combination.push_back(std::make_pair(dihedral, torsion_value_at_this_dihedral));
            }
            container_for_combinations.push_back(new_combination);
        }
    }
}

//This is a wrapper for Python to call.
void Assembly::BuildAssemblyFromCondensedSequence(std::string condensed_sequence, std::string prep_file_path)
{
    PrepFileSpace::PrepFile* prepfile = new PrepFileSpace::PrepFile(prep_file_path);
    this->BuildAssemblyFromCondensedSequence(condensed_sequence, prepfile);
}

//New BuildAssemblyFromCondensedSequence() created by Yao on 06/25/2018. This will replace the old version below.
void Assembly::BuildAssemblyFromCondensedSequence(std::string condensed_sequence, PrepFileSpace::PrepFile* prep_file)
{
//    std::cout << "Building Assembly From Condensed Sequence......" << std::endl;
    CondensedSequenceSpace::CondensedSequence sequence (condensed_sequence);
    CondensedSequenceSpace::CondensedSequence::CondensedSequenceGlycam06ResidueTree glycam06_residues = sequence.GetCondensedSequenceGlycam06ResidueTree();
    
    //    CondensedSequenceSpace::CondensedSequence::CondensedSequenceResidueTree res_tree = sequence.GetCondensedSequenceResidueTree();
    //    CondensedSequenceSpace::CondensedSequence::CondensedSequenceRotamersAndGlycosidicAnglesInfo info = sequence.GetCondensedSequenceRotamersAndGlycosidicAnglesInfo(res_tree);

    Assembly::TemplateAssembly* template_assembly = this-> BuildTemplateAssemblyFromPrepFile (glycam06_residues, prep_file);

    std::map<int, std::pair<CondensedSequenceSpace::CondensedSequenceGlycam06Residue*, Residue*> > glycam06_assembly_residue_map =
            this -> ConvertCondensedSequence2AssemblyResidues (glycam06_residues, template_assembly);

    this -> SetGlycam06ResidueBonding (glycam06_assembly_residue_map);

    std::multimap<int, std::pair<AtomVector*, std::string> > index_dihedral_map = std::multimap<int, std::pair<AtomVector*, std::string> >();
    for (std::map<int,std::pair<CondensedSequenceSpace::CondensedSequenceGlycam06Residue*, Residue*> >::iterator it =
         glycam06_assembly_residue_map.begin(); it != glycam06_assembly_residue_map.end(); it++){
        CondensedSequenceSpace::CondensedSequenceGlycam06Residue* glycam_06_res = it->second.first;
        Residue* corresponding_assembly_residue = it->second.second;
        //The atom and the only atom without a parent is the absolute parent(terminal).
        if (glycam_06_res->GetParentId() == gmml::iNotSet && glycam_06_res->GetName() != "Deoxy"){
            Residue* root = corresponding_assembly_residue;
            //  TURN OFF GEOMETRY OPS
            this->RecursivelySetAngleGeometry(root);
            this->SetDihedralAngleGeometryWithMetadata(); // Removed for current push. Need to tag dihedrals. Also need to resolve clashes.
            //this->RecursivelySetDihedralAngleGeometry(root); // replaced by SetDihedralAngleGeometryWithMetadata. OG 2019.09.
            //          The Recursive function below needs to number all dihedrals, so it needs to know the linkage index at the beginning.
            //          Linkage index is incremented inside function once a linkage has been processed.


            // OG 2019.09:
            // This next part should be replaced. The metadata contains the tags, and should be edited to include a linkage_index
            // However, the ResolveClashes function uses the index_dihedral_map, so that must be replaced at the same time.


//            int linkage_index = 0;
//            this->RecursivelyTagDihedrals(root, index_dihedral_map, linkage_index);
            break;
        }
    }

//    //Find and resolve clashes below(crudely)
//    ResidueVector clashing_residues = this->FindClashingResidues();
//    std::vector<ResidueVector> clashing_residue_parent_paths = this -> FindPathToCommonAncestors(clashing_residues);
//    this->ResolveClashes(clashing_residue_parent_paths, index_dihedral_map);
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
                GeometryTopology::CoordinateVector cartesian_coordinate_list = GeometryTopology::CoordinateVector();
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
                PrepFileSpace::PrepFileAtomVector prep_atoms = prep_residue->GetAtoms();
                for(PrepFileSpace::PrepFileAtomVector::iterator it1 = prep_atoms.begin(); it1 != prep_atoms.end(); it1++)
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
                        GeometryTopology::CoordinateVector coordinate_list;
                        //std::vector<GeometryTopology::Coordinate*> coordinate_list = std::vector<GeometryTopology::Coordinate*>();
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
                        double bond_length = prep_atom->GetBondLength();
                        double angle_value = prep_atom->GetAngle();
                        double dihedral_value = prep_atom->GetDihedral();
                        coordinate = coordinate->ConvertInternalCoordinate2CartesianCoordinate(
                                    coordinate_list, bond_length, angle_value, dihedral_value);
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
//                std::cout << "Residue " << glycam06_residue_name << " has not been found in the database" << std::endl;
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
//        std::cout << "Building assembly from " << sequence << " failed." << std::endl;
    }


}

void Assembly::GenerateRotamersForCondensedSequence (Assembly* working_assembly, CondensedSequenceSpace::CondensedSequence::CondensedSequenceRotamersAndGlycosidicAnglesInfo 
                                                     rotamers_glycosidic_angles_info, std::multimap<int, std::pair<AtomVector*, std::string> >& index_dihedral_map)
{
    std::vector<std::pair<AtomVector*, std::vector<double> > > all_dihedral_rotation_values = std::vector<std::pair<AtomVector*, std::vector<double> > >();
    typedef std::multimap<int, std::pair<AtomVector*, std::string> > index_torsion_map;
    //Go through each linkage, look at its phi,psi and omega(if exists)
    for(unsigned int i = 0; i < rotamers_glycosidic_angles_info.size(); i++){
        CondensedSequenceSpace::RotamersAndGlycosidicAnglesInfo* linkage_info = rotamers_glycosidic_angles_info[i].second;
        int linkage_index = linkage_info->GetLinkageIndex();
        std::pair<index_torsion_map::iterator, index_torsion_map::iterator > current_linkage_dihedral_range = index_dihedral_map.equal_range(linkage_index);
        //First, check for torsions explicitly set by the user. Apply user settting instead of default setting.
        std::vector<std::pair<std::string, double> > user_defined_angles = linkage_info->GetEnabledGlycosidicAngles();
        //Use a string to record what type of torsions are being set. For example, if phi and psi are explicitly set, user_defined_angle_types = "phipsi"
        std::string user_defined_angle_types = std::string();
        for (unsigned int j = 0; j < user_defined_angles.size(); j++){
            std::string& dihedral_type = user_defined_angles[j].first;
            double& rotation_value = user_defined_angles[j].second;
            for (index_torsion_map::iterator it = current_linkage_dihedral_range.first; it != current_linkage_dihedral_range.second; it++){
                std::string& map_dihedral_type = it->second.second;
                if (map_dihedral_type == dihedral_type){
                    AtomVector* dihedral_atoms = it->second.first;
                    std::vector<double> rotation_values = std::vector<double>(1,rotation_value);
                    all_dihedral_rotation_values.push_back(std::make_pair(dihedral_atoms, rotation_values));
                    user_defined_angle_types += map_dihedral_type;
                }
            }
        }
        //Then, go through all rotamers, including user-defined and default torsions
        std::vector<std::pair<std::string,std::vector<std::string> > > selected_rotamers = linkage_info->GetSelectedRotamers();
        for (unsigned int j = 0; j < selected_rotamers.size(); j++){
            std::string& dihedral_type = selected_rotamers[j].first;
            //If a torsion is not being explicitly set:
            if (user_defined_angle_types.find(dihedral_type) == std::string::npos){
                std::vector<std::string>& rotamer_names = selected_rotamers[j].second;
                std::vector<double> rotation_values = std::vector<double>();
                //This should be put into metadata one day.
                for (unsigned int k = 0; k < rotamer_names.size(); k++){
                    if (rotamer_names[k] == "gg"){
                        rotation_values.push_back(-60.0);
                    }
                    if (rotamer_names[k] == "gt"){
                        rotation_values.push_back(60.0);
                    }
                    if (rotamer_names[k] == "tg"){
                        rotation_values.push_back(180.0);
                    }
                    if (rotamer_names[k] == "g"){
                        rotation_values.push_back(60.0);
                    }
                    if (rotamer_names[k] == "t"){
                        rotation_values.push_back(180.0);
                    }
                    if (rotamer_names[k] == "-g"){
                        rotation_values.push_back(-60.0);
                    }
                }
                for (index_torsion_map::iterator it = current_linkage_dihedral_range.first; it != current_linkage_dihedral_range.second; it++){
                    std::string& map_dihedral_type = it->second.second;
                    if (map_dihedral_type == dihedral_type){
                        AtomVector* dihedral_atoms = it->second.first;
                        all_dihedral_rotation_values.push_back(std::make_pair(dihedral_atoms,rotation_values));
                    }
                }
            }
        }

    }
    //Generate All Rotation Combinations
    typedef std::vector<std::pair<AtomVector*, double> > combination;
    std::vector<combination> all_rotation_combinations = std::vector<combination>();
    std::vector<double> angle_index_per_dihedral = std::vector<double> (all_dihedral_rotation_values.size(), gmml::dNotSet);
    working_assembly ->GenerateAllTorsionCombinations(all_dihedral_rotation_values, 0, all_rotation_combinations, angle_index_per_dihedral);
    //Rotate according to each combination
    std::vector<GeometryTopology::CoordinateVector> rotamer_coordinate_sets = std::vector<GeometryTopology::CoordinateVector>();
    for (unsigned int i = 0; i < all_rotation_combinations.size(); i++){
        combination& rotamer_rotation_set = all_rotation_combinations[i];
        for (unsigned int j = 0; j < rotamer_rotation_set.size(); j++){
            std::pair<AtomVector*, double>& dihedral_rotation_value_pair = rotamer_rotation_set[j];
            AtomVector* dihedral_atoms = dihedral_rotation_value_pair.first;
            double rotation_value = dihedral_rotation_value_pair.second;
            working_assembly ->SetDihedral(dihedral_atoms->at(0), dihedral_atoms->at(1), dihedral_atoms->at(2), dihedral_atoms->at(3), rotation_value);
        }
        //Copy current coordinate objects
        GeometryTopology::CoordinateVector new_rotamer_coordinate_set = GeometryTopology::CoordinateVector();
        AtomVector all_atoms_of_assembly = working_assembly->GetAllAtomsOfAssembly();
        for (unsigned int j = 0; j < all_atoms_of_assembly.size(); j++){
            GeometryTopology::Coordinate* original_coordinate = all_atoms_of_assembly[j]->GetCoordinates().at(0);
            GeometryTopology::Coordinate* copied_coordinate = new GeometryTopology::Coordinate(original_coordinate);
            new_rotamer_coordinate_set.push_back(copied_coordinate);
        }
        rotamer_coordinate_sets.push_back(new_rotamer_coordinate_set);
    }
    
    AtomVector all_atoms_of_assembly = working_assembly->GetAllAtomsOfAssembly();
    //Empty current atoms coordinates
    for (unsigned int i = 0; i < all_atoms_of_assembly.size(); i++){
        all_atoms_of_assembly[i]->SetCoordinates(GeometryTopology::CoordinateVector());
    }
    //Add each rotamer coordinate set one by one to atoms in assembly
    for (unsigned int i = 0; i < rotamer_coordinate_sets.size(); i++){
        GeometryTopology::CoordinateVector& rotamer_set = rotamer_coordinate_sets[i];
        for (unsigned int j = 0; j < rotamer_set.size(); j++){
            all_atoms_of_assembly[j]->AddCoordinate(rotamer_set[j]);
        }
    }
    
    
}

/** ***************************************************************************
 *  *          Build from PDB files
 *   * ************************************************************************* **/

void Assembly::BuildAssemblyFromPdbFile(std::string pdb_file_path, std::vector<std::string> amino_lib_files, std::vector<std::string> glycam_lib_files,
                                        std::vector<std::string> other_lib_files, std::vector<std::string> prep_files, std::string parameter_file)
{
    // std::cout << "Building assembly from pdb file ..." << std::endl;
    gmml::log(__LINE__, __FILE__, gmml::INF, "Building assembly from pdb file ...");
    //    std::cout << "Reading PDB file into PdbFileSpace::PdbFile structure." << std::endl;
    PdbFileSpace::PdbFile* pdb_file=NULL;
    try
    {
        gmml::log(__LINE__, __FILE__, gmml::INF, "Reading PDB file into PdbFileSpace::PdbFile structure ...");
        pdb_file = new PdbFileSpace::PdbFile(pdb_file_path);
    }
    catch(PdbFileSpace::PdbFileProcessingException &ex)
    {
//        std::cout << "Generating PdbFileSpace::PdbFile structure from " << pdb_file_path << "failed." << std::endl;
    }
    this->BuildAssemblyFromPdbFile(pdb_file, amino_lib_files, glycam_lib_files, other_lib_files, prep_files, parameter_file);
}

void Assembly::BuildAssemblyFromAtomStream(std::stringstream& atomStream)
{
//    std::cout << "Building assembly from stringstream ..." << std::endl;
    //    std::cout << "Reading PDB file into PdbFileSpace::PdbFile structure." << std::endl;
    PdbFileSpace::PdbFile* pdb_file=NULL;
    try
    {
        gmml::log(__LINE__, __FILE__, gmml::INF, "Reading stringstream into PdbFileSpace::PdbFile structure ...");
        pdb_file = new PdbFileSpace::PdbFile(atomStream);
    }
    catch(PdbFileSpace::PdbFileProcessingException &ex)
    {
//        std::cout << "Generating PdbFileSpace::PdbFile structure from stringstream failed." << std::endl;
    }
    std::vector<std::string> amino_lib_files;
    std::vector<std::string> glycam_lib_files;
    std::vector<std::string> other_lib_files;
    std::vector<std::string> prep_files;
    std::string parameter_file;
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
                std::string atom_element = atom->GetAtomElementSymbol();
                new_atom->SetElementSymbol(atom_element);
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
            if( (residue->GetName() == "ALA") ||
                    (residue->GetName() == "ARG") ||
                    (residue->GetName() == "ASN") ||
                    (residue->GetName() == "ASP") ||
                    (residue->GetName() == "CYS") ||
                    (residue->GetName() == "GLU") ||
                    (residue->GetName() == "GLN") ||
                    (residue->GetName() == "GLY") ||
                    (residue->GetName() == "HIS") ||
                    (residue->GetName() == "ILE") ||
                    (residue->GetName() == "LEU") ||
                    (residue->GetName() == "LYS") ||
                    (residue->GetName() == "MET") ||
                    (residue->GetName() == "PHE") ||
                    (residue->GetName() == "PRO") ||
                    (residue->GetName() == "SER") ||
                    (residue->GetName() == "THR") ||
                    (residue->GetName() == "TRP") ||
                    (residue->GetName() == "TYR") ||
                    (residue->GetName() == "VAL") )
            {
                residue->SetChemicalType("Amino Acid");
            }
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
//    std::cout << "Building assembly from pdbqt file ..." << std::endl;
//    std::cout << "Reading PDBQT file into PdbqtFileSpace::PdbqtFile structure." << std::endl;
    PdbqtFileSpace::PdbqtFile pdbqt_file;
    try
    {
        gmml::log(__LINE__, __FILE__, gmml::INF, "Reading PDBQT file into PdbqtFileSpace::PdbqtFile structure ...");
        pdbqt_file = PdbqtFileSpace::PdbqtFile(pdbqt_file_path);
    }
    catch(PdbqtFileSpace::PdbqtFileProcessingException &ex)
    {
//        std::cout << "Generating PdbqtFileSpace::PdbqtFile structure from " << pdbqt_file_path << "failed." << std::endl;
    }
    this->BuildAssemblyFromPdbqtFile(&pdbqt_file, parameter_file);
}


void Assembly::BuildAssemblyFromPdbqtFile(PdbqtFileSpace::PdbqtFile *pdbqt_file, std::string parameter_file)
{
//    std::cout << "Building assembly from pdbqt file ..." << std::endl;
    gmml::log(__LINE__, __FILE__, gmml::INF, "Building assembly from pdbqt file ...");
    //Source of metadata: http://autodock.scripps.edu/faqs-help/manual/autodock-4-2-user-guide/AutoDock4.2_UserGuide.pdf
    std::map<std::string, std::string> ad_type_element_symbol_map = {
            {"C","C"}, {"A","C"}, {"HD","H"}, {"H","H"}, {"HS","H"}, {"N","N"}, {"NA","N"}, {"NS","N"}, {"NR", "N"}, {"OA","O"}, {"OS","O"}, {"F","F"}, {"Mg","Mg"}, {"MG","Mg"}, {"P","P"}, {"SA","S"}, 
	    {"S","S"}, {"Cl","Cl"}, {"CL","Cl"}, {"CA","Ca"}, {"Ca","Ca"}, {"MN","Mn"}, {"Mn","Mn"}, {"FE","Fe"}, {"Fe","Fe"}, {"ZN","Zn"}, {"Zn","Zn"}, {"BR","Br"}, {"Br","Br"}, {"I","I"}, {"Z","Z"}, 
	    {"G","C"}, {"GA","C"}, {"J","C"}, {"Q","C"}
    };
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
		std::string autodock_type = atom->GetAtomType();
                new_atom->MolecularDynamicAtom::SetAtomType(autodock_type);
		if (ad_type_element_symbol_map.find(autodock_type) != ad_type_element_symbol_map.end()){
		    new_atom->SetElementSymbol(ad_type_element_symbol_map[autodock_type]);
		}
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
//    std::cout << "Building assembly from an AMBER parameter-topology file ..." << std::endl;
//    std::cout << "Reading AMBER parameter-topology file into TopologyFileSpace::TopologyFile structure." << std::endl;
    gmml::log(__LINE__, __FILE__, gmml::INF, "Reading AMBER parameter-topology file into TopologyFileSpace::TopologyFile structure.");
    TopologyFileSpace::TopologyFile topology_file;
    try
    {
        topology_file = TopologyFileSpace::TopologyFile(topology_file_path);
    }
    catch(TopologyFileSpace::TopologyFileProcessingException &ex)
    {
//        std::cout << "Generating TopologyFileSpace::TopologyFile structure from " << topology_file_path << "failed." << std::endl;
    }
    this->BuildAssemblyFromTopologyFile(&topology_file, parameter_file);
}

void Assembly::BuildAssemblyFromTopologyFile(TopologyFileSpace::TopologyFile *topology_file, std::string parameter_file)
{
//    std::cout << "Building assembly from topology file ..." << std::endl;
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
//    std::cout << "Building assembly from AMBER library/off file ..." << std::endl;
//   std::cout << "Reading AMBER library/off file into LibraryFileSpace::LibraryFile structure." << std::endl;
    gmml::log(__LINE__, __FILE__, gmml::INF, "Reading AMBER library/off file into LibraryFileSpace::LibraryFile structure ...");
    LibraryFileSpace::LibraryFile library_file;
    try
    {
        library_file = LibraryFileSpace::LibraryFile(library_file_path);
    }
    catch(LibraryFileSpace::LibraryFileProcessingException &ex)
    {
//        std::cout << "Generating LibraryFileSpace::LibraryFile structure from " << library_file_path << "failed." << std::endl;
    }
    this->BuildAssemblyFromLibraryFile(&library_file, parameter_file);
}



void Assembly::BuildAssemblyFromLibraryFile(LibraryFileSpace::LibraryFile *library_file, std::string parameter_file)
{
//    std::cout << "Building assembly from library file ..." << std::endl;
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
//    std::cout << "Building assembly from AMBER parameter-topology and unknown-style coordinate files ..." << std::endl;
//    std::cout << "Reading AMBER parameter-topology and unknown-style coordinate files into their file structure." << std::endl;
    gmml::log(__LINE__, __FILE__, gmml::INF, "Reading AMBER parameter-topology and unknown-style coordinate files into their file structure.");
    TopologyFileSpace::TopologyFile topology_file;
    CoordinateFileSpace::CoordinateFile coordinate_file;
    try
    {
        topology_file = TopologyFileSpace::TopologyFile(topology_file_path);
    }
    catch(TopologyFileSpace::TopologyFileProcessingException &ex)
    {
//        std::cout << "Generating TopologyFileSpace::TopologyFile structure from " << topology_file_path << "failed." << std::endl;
    }
    try
    {
        coordinate_file = CoordinateFileSpace::CoordinateFile(coordinate_file_path);
    }
    catch(CoordinateFileSpace::CoordinateFileProcessingException &ex)
    {
//        std::cout << "Generating CoordinateFileSpace::CoordinateFile structure from " << coordinate_file_path << "failed." << std::endl;
    }
    this->BuildAssemblyFromTopologyCoordinateFile(&topology_file, &coordinate_file, parameter_file);
}

void Assembly::BuildAssemblyFromTopologyCoordinateFile(TopologyFileSpace::TopologyFile *topology_file, CoordinateFileSpace::CoordinateFile *coordinate_file, std::string parameter_file)
{
//    std::cout << "Building assembly from topology and coordinate files ..." << std::endl;
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

            GeometryTopology::CoordinateVector coord_file_coordinates = coordinate_file->GetCoordinates();
            //std::vector<GeometryTopology::Coordinate*> coord_file_coordinates = coordinate_file->GetCoordinates();
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
//    std::cout << "Building assembly from prep file ..." << std::endl;
//    std::cout << "Reading Prep file into PrepFileSpace::PrepFile structure." << std::endl;
    gmml::log(__LINE__, __FILE__, gmml::INF, "Reading Prep file into PrepFileSpace::PrepFile structure ...");
    PrepFileSpace::PrepFile prep_file;
    try
    {
        prep_file = PrepFileSpace::PrepFile(prep_file_path);
    }
    catch(PrepFileSpace::PrepFileProcessingException &ex)
    {
//        std::cout << "Generating PrepFileSpace::PrepFile structure from " << prep_file_path << "failed." << std::endl;
    }
    this->BuildAssemblyFromPrepFile(&prep_file, parameter_file);
}



void Assembly::BuildAssemblyFromPrepFile(PrepFileSpace::PrepFile *prep_file, std::string parameter_file)
{
//    std::cout << "Building assembly from prep file ..." << std::endl;
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
        GeometryTopology::CoordinateVector cartesian_coordinate_list = GeometryTopology::CoordinateVector();
        int head_atom_index = (int) INFINITY;
        int tail_atom_index = (int) -INFINITY;
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
        PrepFileSpace::PrepFileAtomVector prep_atoms = prep_residue->GetAtoms();
        PrepFileSpace::PrepFileAtomVector parent_atoms = prep_residue->GetAtomsParentVector();

        for(PrepFileSpace::PrepFileAtomVector::iterator it1 = prep_atoms.begin(); it1 != prep_atoms.end(); it1++)
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
                GeometryTopology::CoordinateVector coordinate_list;
                //std::vector<GeometryTopology::Coordinate*> coordinate_list = std::vector<GeometryTopology::Coordinate*>();
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
                GeometryTopology::Coordinate* coordinate = new GeometryTopology::Coordinate();
                double bond_length = prep_atom->GetBondLength();
                double angle_value = prep_atom->GetAngle();
                double dihedral_value = prep_atom->GetDihedral();
                coordinate = coordinate->ConvertInternalCoordinate2CartesianCoordinate(
                            coordinate_list, bond_length, angle_value, dihedral_value);
                cartesian_coordinate_list.push_back(coordinate);
                //original
                //GeometryTopology::Coordinate* coordinate = GeometryTopology::Coordinate::ConvertInternalCoordinate2CartesianCoordinate(coordinate_list, prep_atom->GetBondLength(),
                //prep_atom->GetAngle(), prep_atom->GetDihedral());
                //cartesian_coordinate_list.push_back(coordinate);
                //

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
