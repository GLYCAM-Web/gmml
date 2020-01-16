#include <math.h>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <set>
#include <queue>
#include <stack>

#include "../../../../includes/MolecularModeling/assembly.hpp"
#include "../../../../includes/MolecularModeling/residue.hpp"
#include "../../../../includes/MolecularModeling/atom.hpp"
#include "../../../../includes/MolecularModeling/atomnode.hpp"
#include "../../../../includes/InputSet/CondensedSequenceSpace/condensedsequence.hpp"
#include "../../../../includes/InputSet/CondensedSequenceSpace/condensedsequenceresidue.hpp"
#include "../../../../includes/InputSet/CondensedSequenceSpace/condensedsequenceglycam06residue.hpp"
#include "../../../../includes/InputSet/TopologyFileSpace/topologyfile.hpp"
#include "../../../../includes/InputSet/TopologyFileSpace/topologyassembly.hpp"
#include "../../../../includes/InputSet/TopologyFileSpace/topologyresidue.hpp"
#include "../../../../includes/InputSet/TopologyFileSpace/topologyatom.hpp"
#include "../../../../includes/InputSet/TopologyFileSpace/topologybond.hpp"
#include "../../../../includes/InputSet/TopologyFileSpace/topologybondtype.hpp"
#include "../../../../includes/InputSet/TopologyFileSpace/topologyangle.hpp"
#include "../../../../includes/InputSet/TopologyFileSpace/topologyangletype.hpp"
#include "../../../../includes/InputSet/TopologyFileSpace/topologydihedral.hpp"
#include "../../../../includes/InputSet/TopologyFileSpace/topologydihedraltype.hpp"
#include "../../../../includes/InputSet/TopologyFileSpace/topologyatompair.hpp"
#include "../../../../includes/InputSet/CoordinateFileSpace/coordinatefile.hpp"
#include "../../../../includes/ParameterSet/PrepFileSpace/prepfile.hpp"
#include "../../../../includes/ParameterSet/PrepFileSpace/prepfileresidue.hpp"
#include "../../../../includes/ParameterSet/PrepFileSpace/prepfileatom.hpp"
#include "../../../../includes/InputSet/PdbFileSpace/pdbfile.hpp"
#include "../../../../includes/InputSet/PdbFileSpace/pdbtitlesection.hpp"
#include "../../../../includes/InputSet/PdbFileSpace/pdbmodelcard.hpp"
#include "../../../../includes/InputSet/PdbFileSpace/pdbmodelsection.hpp"
#include "../../../../includes/InputSet/PdbFileSpace/pdbmodelresidueset.hpp"
#include "../../../../includes/InputSet/PdbFileSpace/pdbatomcard.hpp"
#include "../../../../includes/InputSet/PdbFileSpace/pdbheterogenatomsection.hpp"
#include "../../../../includes/InputSet/PdbFileSpace/pdbatomsection.hpp"
#include "../../../../includes/InputSet/PdbFileSpace/pdbconnectsection.hpp"
#include "../../../../includes/InputSet/PdbFileSpace/pdblinkcard.hpp"
#include "../../../../includes/InputSet/PdbFileSpace/pdblinksection.hpp"
#include "../../../../includes/InputSet/PdbFileSpace/pdblinkcardresidue.hpp"
#include "../../../../includes/InputSet/PdbFileSpace/pdbfileprocessingexception.hpp"
#include "../../../../includes/InputSet/PdbqtFileSpace/pdbqtfile.hpp"
#include "../../../../includes/InputSet/PdbqtFileSpace/pdbqtatom.hpp"
#include "../../../../includes/InputSet/PdbqtFileSpace/pdbqtmodel.hpp"
#include "../../../../includes/InputSet/PdbqtFileSpace/pdbqtmodelcard.hpp"
#include "../../../../includes/InputSet/PdbqtFileSpace/pdbqtatomcard.hpp"
#include "../../../../includes/InputSet/PdbqtFileSpace/pdbqtmodelresidueset.hpp"
#include "../../../../includes/ParameterSet/LibraryFileSpace/libraryfile.hpp"
#include "../../../../includes/ParameterSet/LibraryFileSpace/libraryfileatom.hpp"
#include "../../../../includes/ParameterSet/LibraryFileSpace/libraryfileresidue.hpp"
#include "../../../../includes/ParameterSet/ParameterFileSpace/parameterfile.hpp"
#include "../../../../includes/ParameterSet/ParameterFileSpace/parameterfilebond.hpp"
#include "../../../../includes/ParameterSet/ParameterFileSpace/parameterfileangle.hpp"
#include "../../../../includes/ParameterSet/ParameterFileSpace/parameterfiledihedral.hpp"
#include "../../../../includes/ParameterSet/ParameterFileSpace/parameterfiledihedralterm.hpp"
#include "../../../../includes/ParameterSet/ParameterFileSpace/parameterfileatom.hpp"
#include "../../../../includes/ParameterSet/PrepFileSpace/prepfile.hpp"
#include "../../../../includes/ParameterSet/PrepFileSpace/prepfileresidue.hpp"
#include "../../../../includes/ParameterSet/PrepFileSpace/prepfileatom.hpp"
#include "../../../../includes/utils.hpp"
#include "../../../../includes/common.hpp"
#include "../../../../includes/GeometryTopology/grid.hpp"
#include "../../../../includes/GeometryTopology/cell.hpp"

#include <unistd.h>
#include <errno.h>
#include <string.h>

using MolecularModeling::Assembly;

//////////////////////////////////////////////////////////
//                       FUNCTIONS                      //
//////////////////////////////////////////////////////////
PrepFileSpace::PrepFile* Assembly::BuildPrepFileStructureFromAssembly(std::string parameter_file_path)
{
//    std::cout << "Creating prep file ..." << std::endl;
    gmml::log(__LINE__, __FILE__, gmml::INF, "Creating prep file ...");
    PrepFileSpace::PrepFile* prep_file = new PrepFileSpace::PrepFile();
    ResidueVector assembly_residues = this->GetAllResiduesOfAssembly();
    PrepFileSpace::PrepFile::ResidueMap prep_residues = PrepFileSpace::PrepFile::ResidueMap();
    ParameterFileSpace::ParameterFile* parameter_file = new ParameterFileSpace::ParameterFile(parameter_file_path);
    ParameterFileSpace::ParameterFile::DihedralMap dihedrals = parameter_file->GetDihedrals();
    for(ResidueVector::iterator it = assembly_residues.begin(); it != assembly_residues.end(); it++)
    {
        Residue* assembly_residue = *it;
        std::vector<std::string> inserted_improper_dihedral_types = std::vector<std::string>();
        std::vector<std::vector<std::string> > inserted_improper_dihedrals = std::vector<std::vector<std::string> >();
        PrepFileSpace::PrepFileResidue* prep_residue = new PrepFileSpace::PrepFileResidue();
        PrepFileSpace::PrepFileAtomVector prep_atoms = PrepFileSpace::PrepFileAtomVector();
        prep_residue->SetTitle(assembly_residue->GetName());
        prep_residue->SetName(assembly_residue->GetName());
        prep_residue->SetGeometryType(PrepFileSpace::kGeometryCorrect);
        prep_residue->SetCoordinateType(PrepFileSpace::kINT);
        prep_residue->SetDummyAtomOmission(PrepFileSpace::kOmit);
        prep_residue->SetDummyAtomType("DU");
        prep_residue->SetDummyAtomPosition(PrepFileSpace::kPositionBeg);
        prep_residue->SetOutputFormat(PrepFileSpace::kFormatted);

        AtomVector assembly_atoms = assembly_residue->GetAtoms();
        GeometryTopology::CoordinateVector cartesian_coordinate_list;
        int atom_index = 1;
        PrepFileSpace::PrepFileResidue::Loop loops = PrepFileSpace::PrepFileResidue::Loop();
        std::vector<int> bond_index = std::vector<int>();
        for(int i = 0; i < gmml::DEFAULT_DUMMY_ATOMS; i ++)
        {
            PrepFileSpace::PrepFileAtom* dummy_atom = new PrepFileSpace::PrepFileAtom();
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
                cartesian_coordinate_list.push_back(new GeometryTopology::Coordinate(0.0, 0.0, 0.0)); //-5.0, -5.0, -5.0
            }
            else if(i <= 1)
            {
                dummy_atom->SetBondLength(gmml::maxCutOff);
                dummy_atom->SetAngle(0.0);
                dummy_atom->SetDihedral(0.0);
                cartesian_coordinate_list.push_back(new GeometryTopology::Coordinate(gmml::maxCutOff, 0.0, 0.0)); //1.522, 0.0, 0.0 ; -4.5, -5.0, -5.0
            }
            else
            {
                dummy_atom->SetBondLength(gmml::maxCutOff);
                dummy_atom->SetAngle(90.0);
                dummy_atom->SetDihedral(0.0);
                cartesian_coordinate_list.push_back(new GeometryTopology::Coordinate(gmml::maxCutOff, gmml::maxCutOff, 0.0)); //2.43347, 1.34044, 0.0 ; -4.7, -4.2, -4.0
            }
            dummy_atom->SetCharge(0.0);
            dummy_atom->SetTopologicalType(gmml::kTopTypeM);
            atom_index++;
            prep_atoms.push_back(dummy_atom);
        }
        std::vector<gmml::TopologicalType> residue_topological_types = GetAllTopologicalTypesOfAtomsOfResidue(assembly_atoms, loops, bond_index);
        for(AtomVector::iterator it1 = assembly_atoms.begin(); it1 != assembly_atoms.end(); it1++)
        {
            GeometryTopology::CoordinateVector coordinate_list;

            Atom* assembly_atom = (*it1);
            PrepFileSpace::PrepFileAtom* prep_atom = new PrepFileSpace::PrepFileAtom();
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
            prep_atom->SetTopologicalType(residue_topological_types.at(atom_index - gmml::DEFAULT_DUMMY_ATOMS - 1));
            prep_atom->SetType(assembly_atom->GetAtomType());
            int parent_index = prep_atom->GetBondIndex() - 1;
            int grandparent_index = prep_atom->GetAngleIndex() - 1;
            int great_grandparent_index = prep_atom->GetDihedralIndex() - 1;
            coordinate_list.push_back(cartesian_coordinate_list.at(great_grandparent_index));
            coordinate_list.push_back(cartesian_coordinate_list.at(grandparent_index));
            coordinate_list.push_back(cartesian_coordinate_list.at(parent_index));
            cartesian_coordinate_list.push_back(assembly_atom->GetCoordinates().at(model_index_));

            GeometryTopology::Coordinate* internal_coordinate = new GeometryTopology::Coordinate();
            internal_coordinate = internal_coordinate->GeometryTopology::Coordinate::ConvertCartesianCoordinate2InternalCoordinate(
                assembly_atom->GetCoordinates().at(model_index_), coordinate_list);
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

void Assembly::ExtractPrepImproperDihedralTypesFromAssembly(Atom *assembly_atom, std::vector<std::string>
                                                            & inserted_improper_dihedral_types, ParameterFileSpace::ParameterFile::DihedralMap& dihedrals)
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
            std::vector<std::vector<std::string> > all_improper_dihedrals_atom_type_permutations = CreateAllAtomTypePermutationsforImproperDihedralType(neighbor1->GetAtomType(), neighbor2->GetAtomType(),
                                                                                                                                         neighbor3->GetAtomType(), assembly_atom->GetAtomType());
            for(std::vector<std::vector<std::string> >::iterator it = all_improper_dihedrals_atom_type_permutations.begin(); it != all_improper_dihedrals_atom_type_permutations.end(); it++)
            {
                std::vector<std::string> improper_dihedral_permutation = (*it);
                if(dihedrals[improper_dihedral_permutation] != NULL)
                {
                    std::stringstream ss;
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

void Assembly::ExtractPrepImproperDihedralsFromAssembly(Atom *assembly_atom, std::vector<std::vector<std::string> >& inserted_improper_dihedrals, std::vector<std::string> inserted_improper_dihedral_types,
                                                        ParameterFileSpace::ParameterFile::DihedralMap &dihedrals)
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
            std::vector<std::vector<std::string> > all_improper_dihedrals_atom_type_permutations = CreateAllAtomTypePermutationsforImproperDihedralType(neighbor1->GetAtomType(), neighbor2->GetAtomType(),
                                                                                                                                         neighbor3->GetAtomType(), assembly_atom->GetAtomType());
            for(std::vector<std::vector<std::string> >::iterator it = all_improper_dihedrals_atom_type_permutations.begin(); it != all_improper_dihedrals_atom_type_permutations.end(); it++)
            {
                std::vector<std::string> improper_dihedral_permutation = (*it);
                std::stringstream sss;
                sss << improper_dihedral_permutation.at(0) << "_" << improper_dihedral_permutation.at(1) << "_" << improper_dihedral_permutation.at(2) << "_" << improper_dihedral_permutation.at(3);
                if(find(inserted_improper_dihedral_types.begin(), inserted_improper_dihedral_types.end(), sss.str()) != inserted_improper_dihedral_types.end())
                {
                    std::vector<std::string> dihedral1 = std::vector<std::string>();
                    std::vector<std::string> dihedral2 = std::vector<std::string>();
                    std::vector<std::string> dihedral3 = std::vector<std::string>();
                    std::vector<std::string> reverse_dihedral1 = std::vector<std::string>();
                    std::vector<std::string> reverse_dihedral2 = std::vector<std::string>();
                    std::vector<std::string> reverse_dihedral3 = std::vector<std::string>();

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
                        ParameterFileSpace::ParameterFileDihedral* parameter_file_dihedral = NULL;
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

std::vector<gmml::TopologicalType> Assembly::GetAllTopologicalTypesOfAtomsOfResidue(AtomVector assembly_atoms, PrepFileSpace::PrepFileResidue::Loop& loops,
                                                                         std::vector<int> & bond_index, int dummy_atoms)
{
    std::vector<gmml::TopologicalType> topological_types = std::vector<gmml::TopologicalType>();
    std::vector<bool> visited = std::vector<bool>();
    std::vector<int> stack = std::vector<int>();
    for(AtomVector::iterator it = assembly_atoms.begin(); it != assembly_atoms.end(); it++)
    {
        topological_types.push_back(gmml::kTopTypeM);
        bond_index.push_back(0);
        visited.push_back(false);
    }

    for(AtomVector::iterator it = assembly_atoms.begin(); it != assembly_atoms.end(); it++)
    {
        Atom* assembly_atom = *it;
        int index = distance(assembly_atoms.begin(), it);

        // Finding all neighbors' indices of current atom in the assembly
        AtomNode* node = assembly_atom->GetNode();
        std::vector<int> neighbors_atom_index = std::vector<int>();
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
        for(unsigned int i = 0; i < neighbors_atom_index.size(); i++)
            if(visited.at(neighbors_atom_index.at(i)))
                number_of_visited_neighbors++;

        if(stack.empty())                   // Stack is empty
        {
            topological_types.at(index) = gmml::kTopTypeM;
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
                        topological_types.at(index) = gmml::kTopTypeE;
                        visited.at(index) = true;
                        bond_index.at(index + dummy_atoms) = top_stack_atom_index + dummy_atoms + 1;
                        stack.pop_back();
                        break;
                    case 2:
                        switch(number_of_visited_neighbors)
                        {
                            case 1:
                                if(topological_types.at(top_stack_atom_index) == gmml::kTopTypeM)
                                    if(top_stack_index != 0 && stack.at(top_stack_index - 1) == top_stack_atom_index)
                                        topological_types.at(index) = gmml::kTopTypeS;
                                    else
                                        topological_types.at(index) = gmml::kTopTypeM;
                                else
                                    topological_types.at(index) = gmml::kTopTypeS;
                                stack.pop_back();
                                stack.push_back(index);
                                break;
                            case 2:
                                topological_types.at(index) = gmml::kTopTypeE;
                                stack.pop_back();
                                int loop_back_atom_index = 0;
                                for(unsigned int i = 0; i < neighbors_atom_index.size(); i++)
                                    if(neighbors_atom_index.at(i) != top_stack_atom_index && visited.at(neighbors_atom_index.at(i)) == true)
                                        loop_back_atom_index = neighbors_atom_index.at(i);
                                for(unsigned int i = 0; i < stack.size(); i++)
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
                                if(topological_types.at(top_stack_atom_index) == gmml::kTopTypeM)
                                    if(top_stack_index != 0 && stack.at(top_stack_index - 1) == top_stack_atom_index)
                                        topological_types.at(index) = gmml::kTopTypeB;
                                    else
                                        topological_types.at(index) = gmml::kTopTypeM;
                                else
                                    topological_types.at(index) = gmml::kTopTypeB;
                                stack.pop_back();
                                stack.push_back(index);
                                stack.push_back(index);
                                break;
                            case 2:
                                topological_types.at(index) = gmml::kTopTypeS;
                                stack.pop_back();
                                int loop_back_atom_index = 0;
                                for(unsigned int i = 0; i < neighbors_atom_index.size(); i++)
                                    if(neighbors_atom_index.at(i) != top_stack_atom_index && visited.at(neighbors_atom_index.at(i)) == true)
                                        loop_back_atom_index = neighbors_atom_index.at(i);
                                for(unsigned int i = 0; i < stack.size(); i++)
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
                                if(topological_types.at(top_stack_atom_index) == gmml::kTopTypeM)
                                    if(top_stack_index != 0 && stack.at(top_stack_index - 1) == top_stack_atom_index)
                                        topological_types.at(index) = gmml::kTopType3;
                                    else
                                        topological_types.at(index) = gmml::kTopTypeM;
                                else
                                    topological_types.at(index) = gmml::kTopType3;
                                stack.pop_back();
                                stack.push_back(index);
                                stack.push_back(index);
                                stack.push_back(index);
                                break;
                            case 2:
                                topological_types.at(index) = gmml::kTopTypeB;
                                stack.pop_back();
                                int loop_back_atom_index = 0;
                                for(unsigned int i = 0; i < neighbors_atom_index.size(); i++)
                                    if(neighbors_atom_index.at(i) != top_stack_atom_index && visited.at(neighbors_atom_index.at(i)) == true)
                                        loop_back_atom_index = neighbors_atom_index.at(i);
                                for(unsigned int i = 0; i < stack.size(); i++)
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
                    //                    std::cout << "EMPTY" << std::endl;
                    topological_types.at(index) = gmml::kTopTypeM;
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
                            topological_types.at(index) = gmml::kTopTypeE;
                            visited.at(index) = true;
                            bond_index.at(index + dummy_atoms) = stack_neighbor_atom_index + dummy_atoms + 1;
                            stack.erase(stack.begin() + stack_neighbor_index);
                            break;
                        case 2:
                            switch(number_of_visited_neighbors)
                            {
                                case 1:
                                    if(topological_types.at(stack_neighbor_atom_index) == gmml::kTopTypeM)
                                        if(stack_neighbor_index != 0 && stack.at(stack_neighbor_index - 1) == stack_neighbor_atom_index)
                                            topological_types.at(index) = gmml::kTopTypeS;
                                        else
                                            topological_types.at(index) = gmml::kTopTypeM;
                                    else
                                        topological_types.at(index) = gmml::kTopTypeS;
                                    break;
                                    stack.erase(stack.begin() + stack_neighbor_index);
                                    stack.push_back(index);
                                case 2:
                                    topological_types.at(index) = gmml::kTopTypeE;
                                    stack.erase(stack.begin() + stack_neighbor_index);
                                    int loop_back_atom_index = 0;
                                    for(unsigned int i = 0; i < neighbors_atom_index.size(); i++)
                                        if(neighbors_atom_index.at(i) != top_stack_atom_index && visited.at(neighbors_atom_index.at(i)) == true)
                                            loop_back_atom_index = neighbors_atom_index.at(i);
                                    for(unsigned int i = 0; i < stack.size(); i++)
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
                                    if(topological_types.at(stack_neighbor_atom_index) == gmml::kTopTypeM)
                                        if(stack_neighbor_index != 0 && stack.at(stack_neighbor_index - 1) == stack_neighbor_atom_index)
                                            topological_types.at(index) = gmml::kTopTypeB;
                                        else
                                            topological_types.at(index) = gmml::kTopTypeM;
                                    else
                                        topological_types.at(index) = gmml::kTopTypeB;
                                    stack.erase(stack.begin() + stack_neighbor_index);
                                    stack.push_back(index);
                                    stack.push_back(index);
                                    break;
                                case 2:
                                    topological_types.at(index) = gmml::kTopTypeS;
                                    stack.erase(stack.begin() + stack_neighbor_index);
                                    int loop_back_atom_index = 0;
                                    for(unsigned int i = 0; i < neighbors_atom_index.size(); i++)
                                        if(neighbors_atom_index.at(i) != top_stack_atom_index && visited.at(neighbors_atom_index.at(i)) == true)
                                            loop_back_atom_index = neighbors_atom_index.at(i);
                                    for(unsigned int i = 0; i < stack.size(); i++)
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
                                    if(topological_types.at(stack_neighbor_atom_index) == gmml::kTopTypeM)
                                        if(stack_neighbor_index != 0 && stack.at(stack_neighbor_index - 1) == stack_neighbor_atom_index)
                                            topological_types.at(index) = gmml::kTopType3;
                                        else
                                            topological_types.at(index) = gmml::kTopTypeM;
                                    else
                                        topological_types.at(index) = gmml::kTopType3;
                                    stack.erase(stack.begin() + stack_neighbor_index);
                                    stack.push_back(index);
                                    stack.push_back(index);
                                    stack.push_back(index);
                                    break;
                                case 2:
                                    topological_types.at(index) = gmml::kTopTypeB;
                                    stack.erase(stack.begin() + stack_neighbor_index);
                                    int loop_back_atom_index = 0;
                                    for(unsigned int i = 0; i < neighbors_atom_index.size(); i++)
                                        if(neighbors_atom_index.at(i) != top_stack_atom_index && visited.at(neighbors_atom_index.at(i)) == true)
                                            loop_back_atom_index = neighbors_atom_index.at(i);
                                    for(unsigned int i = 0; i < stack.size(); i++)
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
