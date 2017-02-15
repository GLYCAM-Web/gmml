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
#include "../../../../includes/InputSet/CondensedSequenceSpace/condensedsequenceamberprepresidue.hpp"
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
#include "../../../../includes/InputSet/PdbFileSpace/pdbtitlecard.hpp"
#include "../../../../includes/InputSet/PdbFileSpace/pdbmodelcard.hpp"
#include "../../../../includes/InputSet/PdbFileSpace/pdbmodel.hpp"
#include "../../../../includes/InputSet/PdbFileSpace/pdbmodelresidueset.hpp"
#include "../../../../includes/InputSet/PdbFileSpace/pdbatomcard.hpp"
#include "../../../../includes/InputSet/PdbFileSpace/pdbheterogenatomcard.hpp"
#include "../../../../includes/InputSet/PdbFileSpace/pdbatom.hpp"
#include "../../../../includes/InputSet/PdbFileSpace/pdbconnectcard.hpp"
#include "../../../../includes/InputSet/PdbFileSpace/pdblinkcard.hpp"
#include "../../../../includes/InputSet/PdbFileSpace/pdblink.hpp"
#include "../../../../includes/InputSet/PdbFileSpace/pdblinkresidue.hpp"
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
//                       FUNCTIONS                      //
//////////////////////////////////////////////////////////
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

