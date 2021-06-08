#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>

#include <fstream>
#include <set>
#include <queue>
#include <stack>
#include <cstring>
#include <algorithm>

#include "../../../includes/MolecularModeling/assembly.hpp"
#include "../../../includes/MolecularModeling/residue.hpp"
#include "../../../includes/MolecularModeling/atom.hpp"
#include "../../../includes/MolecularModeling/atomnode.hpp"
#include "../../../includes/InputSet/CondensedSequenceSpace/condensedsequence.hpp"
#include "../../../includes/InputSet/CondensedSequenceSpace/condensedsequenceresidue.hpp"
#include "../../../includes/InputSet/CondensedSequenceSpace/condensedsequenceglycam06residue.hpp"
#include "../../../includes/ParameterSet/PrepFileSpace/prepfile.hpp"
#include "../../../includes/ParameterSet/PrepFileSpace/prepfileresidue.hpp"
#include "../../../includes/ParameterSet/PrepFileSpace/prepfileatom.hpp"
#include "../../../includes/ParameterSet/PrepFileSpace/prepfileprocessingexception.hpp"
#include "../../../includes/utils.hpp"
#include "../../../includes/common.hpp"
#include "../../../includes/GeometryTopology/geometrytopology.hpp"
#include "../../../includes/GeometryTopology/coordinate.hpp"
#include "../../../includes/GeometryTopology/grid.hpp"
#include "../../../includes/GeometryTopology/cell.hpp"
#include "../../../includes/MolecularMetadata/GLYCAM/bondlengthbytypepair.hpp"
#include "../../../includes/MolecularMetadata/GLYCAM/amberatomtypeinfo.hpp"

using MolecularModeling::Residue;

//////////////////////////////////////////////////////////
//                       FUNCTIONS                      //
//////////////////////////////////////////////////////////
void Residue::BuildResidueFromPrepFileResidue(PrepFileSpace::PrepFileResidue *prep_residue)
{
//    std::cout << "Building residue from prep residue ..." << std::endl;
//    std::cout << "prep res name: " << prep_residue->GetName() << std::endl;
    gmml::log(__LINE__, __FILE__, gmml::INF, "Building residue from prep residue ...");

    int serial_number = 0;
    std::vector<GeometryTopology::Coordinate*> cartesian_coordinate_list = std::vector<GeometryTopology::Coordinate*>();
    int head_atom_index = (int) INFINITY;
    int tail_atom_index = (int) -INFINITY;
    Atom* head_atom = new Atom();
    Atom* tail_atom = new Atom();
    PrepFileSpace::PrepFileAtomVector prep_atoms = prep_residue->GetAtoms();
    PrepFileSpace::PrepFileAtomVector parent_atoms = prep_residue->GetAtomsParentVector();
    std::map<PrepFileSpace::PrepFileAtom*,Atom*> prep_assembly_atom_map = std::map<PrepFileSpace::PrepFileAtom*,Atom*>();	//associates a prep atom with corresponding assembly atom
    std::map<Atom*,PrepFileSpace::PrepFileAtom*> assembly_prep_parent_atom_map = std::map<Atom*,PrepFileSpace::PrepFileAtom*>();	//associates an assembly atom with corresonding prep atom parent.
    //eventually: assembly atom -> prep parent -> assembly parent
    this->SetName(prep_residue->GetName());
    std::stringstream residue_id;
    //Set id for testing purpose
    //residue_id << prep_residue->GetName() <<"_" << " " << "_" << " " << "_" << "?_" << "?_" << " " << std::endl;
    for(PrepFileSpace::PrepFileAtomVector::iterator it1 = prep_atoms.begin(); it1 != prep_atoms.end(); it1++)
    {
        serial_number++;
        PrepFileSpace::PrepFileAtom* prep_atom = (*it1);
        Atom* assembly_atom = new Atom();
        if(prep_atom->GetType().find("DU") == std::string::npos)
        {
            std::string atom_type = prep_atom->GetType();
            assembly_atom->SetResidue(this);
            this->AddAtom(assembly_atom);
            std::string atom_name = prep_atom->GetName();
            std::string id;
            std::stringstream atom_id;
            assembly_atom->SetName(atom_name);
            assembly_atom->SetNaming("glycam06");

            //Attention: SetAtomType()function is overloaded as MolecularModeling::Atom::SetAtomType() and MolecularModeling::MolecularDynamicAtom::SetAtomType(). You don't really know which one to use.
            //Likewise, GetAtomType() is also overloaded.
            //In my situation, I called MolecularModeling::MolecularModelingAtom::SetAtomType(), but later called MolecularModeling::Atom::GetAtomType(). The result is empty.
            //We need to talk about this later
            assembly_atom->MolecularDynamicAtom::SetAtomType(atom_type);
            assembly_atom->MolecularDynamicAtom::SetCharge(prep_atom->GetCharge());
            assembly_atom->MolecularDynamicAtom::SetMass(gmml::dNotSet);
            assembly_atom->MolecularDynamicAtom::SetRadius(gmml::dNotSet);

// Oliver new code after rework of metadata:
            gmml::MolecularMetadata::GLYCAM::AmberAtomTypeInfoContainer AtomTypeMetaData;
            std::string element = AtomTypeMetaData.GetElementForAtomType(atom_type);
// Oliver rework of Metadata, Oct 2018. Delete this old commented out section if all is well.
//            std::string element = "";
//            int size_of_lookup_map = sizeof(gmml::MolecularMetadata::GLYCAM::Glycam06j1AtomTypes) / sizeof(gmml::MolecularMetadata::GLYCAM::Glycam06j1AtomTypes[0]);
//            for (int i = 0; i < size_of_lookup_map; i++){
//                gmml::MolecularMetadata::GLYCAM::AmberAtomTypeInfo entry = gmml::MolecularMetadata::GLYCAM::Glycam06j1AtomTypes[i];
//                if (atom_type.compare(entry.type_) == 0){
//                    element = entry.element_;
//                }
//            }
            assembly_atom->SetElementSymbol(element);
            int index = std::distance(prep_atoms.begin(), it1);
            PrepFileSpace::PrepFileAtom* parent_prep_atom = parent_atoms.at(index);
            assembly_prep_parent_atom_map[assembly_atom] = parent_prep_atom;
            prep_assembly_atom_map[prep_atom] = assembly_atom;
        }

        if(prep_residue->GetCoordinateType() == PrepFileSpace::kINT)
        {
//serial_number++;
//PrepFileSpace::PrepFileAtom* prep_atom = (*it1);
//Atom* assembly_atom = new Atom();
            std::vector<GeometryTopology::Coordinate*> coordinate_list = std::vector<GeometryTopology::Coordinate*>();
            unsigned int index = distance(prep_atoms.begin(), it1);
            if(index == 0)
            {
//! \todo Add these cout statements to the debugging mechanism once the DebugLevel class (or whatever) is implemented. 
// std::cout << "Starting to get Cartesian coordinates from internal coords for serial number:  " << std::endl;
// std::cout << "    Index is 0 " << std::endl;
// std::cout << "    Serial number:  " << serial_number << " ;  Prep Name:  " << prep_atom->GetName() << " ;  Assembly Name:  " << assembly_atom->GetName() << std::endl;
            }
            if(index == 1)
            {
                int parent_index = parent_atoms.at(index)->GetIndex() - 1;
                GeometryTopology::Coordinate* parent_coordinate = cartesian_coordinate_list.at(parent_index);
                coordinate_list.push_back(parent_coordinate);
// std::cout << "Starting to get Cartesian coordinates from internal coords for serial number:  " << std::endl;
// std::cout << "    Index is 1 " << std::endl;
// std::cout << "    Serial number:  " << serial_number << " ;  Prep Name:  " << prep_atom->GetName() << " ;  Assembly Name:  " << assembly_atom->GetName() << std::endl;
            }
            if(index == 2)
            {
                int parent_index = parent_atoms.at(index)->GetIndex() - 1;
                int grandparent_index = parent_atoms.at(parent_index)->GetIndex() - 1;
                GeometryTopology::Coordinate* grandparent_coordinate = cartesian_coordinate_list.at(grandparent_index);
                GeometryTopology::Coordinate* parent_coordinate = cartesian_coordinate_list.at(parent_index);
                coordinate_list.push_back(grandparent_coordinate);
                coordinate_list.push_back(parent_coordinate);
// std::cout << "Starting to get Cartesian coordinates from internal coords for serial number:  " << std::endl;
// std::cout << "    Index is 2 " << std::endl;
// std::cout << "    Serial number:  " << serial_number << " ;  Prep Name:  " << prep_atom->GetName() << " ;  Assembly Name:  " << assembly_atom->GetName() << std::endl;
            }
            if(index > 2)
            {
                int parent_index = parent_atoms.at(index)->GetIndex() - 1;
                int grandparent_index = parent_atoms.at(parent_index)->GetIndex() - 1;
                int great_grandparent_index = parent_atoms.at(grandparent_index)->GetIndex() - 1;
                GeometryTopology::Coordinate* great_grandparent_coordinate = cartesian_coordinate_list.at(great_grandparent_index);
                GeometryTopology::Coordinate* grandparent_coordinate = cartesian_coordinate_list.at(grandparent_index);
                GeometryTopology::Coordinate* parent_coordinate = cartesian_coordinate_list.at(parent_index);
                coordinate_list.push_back(great_grandparent_coordinate);
                coordinate_list.push_back(grandparent_coordinate);
                coordinate_list.push_back(parent_coordinate);
// std::cout << "Starting to get Cartesian coordinates from internal coords for serial number:  " << std::endl;
// std::cout << "  Index is greater than 2  ;  in prep residue:  " << prep_residue->GetName() << std::endl;
// std::cout << "    Serial number:  " << serial_number << " ;  Prep Name:  " << prep_atom->GetName() << " ;  Assembly Name:  " << assembly_atom->GetName() << std::endl;
// std::cout << "    Parent Name & index  : " << parent_atoms.at(index)->GetName() << "  " << parent_index << std::endl;
// std::cout << "       X:   " << parent_coordinate->GetX() << std::endl;
// std::cout << "       y:   " << parent_coordinate->GetY() << std::endl;
// std::cout << "       Z:   " << parent_coordinate->GetZ() << std::endl;
// std::cout << "    G Parent Name  & index : " << parent_atoms.at(parent_index)->GetName() << "  " << grandparent_index << std::endl;
// std::cout << "       X:   " << grandparent_coordinate->GetX() << std::endl;
// std::cout << "       y:   " << grandparent_coordinate->GetY() << std::endl;
// std::cout << "       Z:   " << grandparent_coordinate->GetZ() << std::endl;
// std::cout << "    G G Parent Name  & index : " << parent_atoms.at(grandparent_index)->GetName()  << "  " << great_grandparent_index << std::endl;
// std::cout << "       X:   " << great_grandparent_coordinate->GetX() << std::endl;
// std::cout << "       y:   " << great_grandparent_coordinate->GetY() << std::endl;
// std::cout << "       Z:   " << great_grandparent_coordinate->GetZ() << std::endl;
            }
// std::cout << "    My Prep Entry Info follows: " << std::endl;
// std::cout << "       Bond:      " << prep_atom->GetBondLength()  << std::endl;
// std::cout << "       Angle:     " << prep_atom->GetAngle()  << std::endl;
// std::cout << "       Dihedral:  " << prep_atom->GetDihedral()  << std::endl;
// std::cout << "    These coordinates are BEFORE the function : " << std::endl;
// std::cout << "the size is:  "  << coordinate_list.size() << std::endl;
// if ( coordinate_list.size() > 2 ) {
// std::cout << "       [0]   X:   " << coordinate_list[0]->GetX() << std::endl;
// std::cout << "             y:   " << coordinate_list[0]->GetY() << std::endl;
// std::cout << "             Z:   " << coordinate_list[0]->GetZ() << std::endl;
// std::cout << "       [1]   X:   " << coordinate_list[1]->GetX() << std::endl;
// std::cout << "             y:   " << coordinate_list[1]->GetY() << std::endl;
// std::cout << "             Z:   " << coordinate_list[1]->GetZ() << std::endl;
// std::cout << "       [2]   X:   " << coordinate_list[2]->GetX() << std::endl;
// std::cout << "             y:   " << coordinate_list[2]->GetY() << std::endl;
// std::cout << "             Z:   " << coordinate_list[2]->GetZ() << std::endl;
// }
            GeometryTopology::Coordinate* coordinate = new GeometryTopology::Coordinate();
            GeometryTopology::Coordinate rawcoordinate;
            double bond_length = prep_atom->GetBondLength();
            double angle_value = prep_atom->GetAngle();
            double dihedral_value = prep_atom->GetDihedral();

            if(index > 200000 )  // normally > 2
            {
                rawcoordinate =  GeometryTopology::get_cartesian_point_from_internal_coords(
                            coordinate_list[0], coordinate_list[1], coordinate_list[2], dihedral_value, angle_value, bond_length);
                // std::cout << "   The NEW RAW coords are:  " << std::endl;
                // std::cout << "         X  :  " << rawcoordinate.GetX() << std::endl;
                // std::cout << "         Y  :  " << rawcoordinate.GetY() << std::endl;
                // std::cout << "         Z  :  " << rawcoordinate.GetZ() << std::endl;
                coordinate->SetX(rawcoordinate.GetX());
                coordinate->SetY(rawcoordinate.GetY());
                coordinate->SetZ(rawcoordinate.GetZ());
            }
            else
            {
                coordinate = coordinate->ConvertInternalCoordinate2CartesianCoordinate(
                            coordinate_list, bond_length, angle_value, dihedral_value);
            }
            cartesian_coordinate_list.push_back(coordinate);
           // GeometryTopology::Coordinate* coordinate = gmml::ConvertInternalCoordinate2CartesianCoordinate(
            //            coordinate_list,
             //           prep_atom->GetBondLength(),
              //          prep_atom->GetAngle(),
               //         prep_atom->GetDihedral());
      //      cartesian_coordinate_list.push_back(coordinate);
            assembly_atom->AddCoordinate(coordinate);
// std::cout << "    My NEW coordinates are : " << std::endl;
// std::cout << "       X:   " << coordinate->GetX() << std::endl;
// std::cout << "       y:   " << coordinate->GetY() << std::endl;
// std::cout << "       Z:   " << coordinate->GetZ() << std::endl;
        }
        else if(prep_residue->GetCoordinateType() == PrepFileSpace::kXYZ)
        {
            assembly_atom->AddCoordinate(new GeometryTopology::Coordinate(prep_atom->GetBondLength(), prep_atom->GetAngle(), prep_atom->GetDihedral()));
        }
        if(prep_atom->GetTopologicalType() == gmml::TopologicalType::kTopTypeM && prep_atom->GetType().find(prep_residue->GetDummyAtomType()) == std::string::npos)
        {
            if(head_atom_index > prep_atom->GetIndex())
            {
                head_atom_index = prep_atom->GetIndex();
                head_atom = assembly_atom; // might be best for this to point to the object rather than copy the pointer
                head_atom->SetName(assembly_atom->GetName());
            }
            if(tail_atom_index < prep_atom->GetIndex())
            {
                tail_atom_index = prep_atom->GetIndex();
                tail_atom = assembly_atom; // might be best for this to point to the object rather than copy the pointer
                tail_atom->SetName(assembly_atom->GetName());
            }
        }
        /*        if(assembly_atom->GetAtomType().find("DU") == std::string::npos)
    {
            this->AddAtom(assembly_atom);
    }*/
    }

    this->AddHeadAtom(head_atom);
    this->AddTailAtom(tail_atom);

    //Convert prep residue loop indices into nodes
    PrepFileSpace::PrepFileResidue::Loop all_loops = prep_residue->GetLoops();
    for (PrepFileSpace::PrepFileResidue::Loop::iterator it = all_loops.begin(); it != all_loops.end(); it++)
    {
        int from = it->first;
        int to = it->second;
        PrepFileSpace::PrepFileAtom* from_prep_atom = NULL;
        PrepFileSpace::PrepFileAtom* to_prep_atom = NULL;
        Atom* from_assembly_atom = NULL;
        Atom* to_assembly_atom = NULL;
        for (unsigned int i = 0; i < prep_atoms.size() ; i++){
            if (prep_atoms[i] -> GetIndex() == from)
            {
                from_prep_atom = prep_atoms[i];
            }
            if (prep_atoms[i] -> GetIndex() == to)
            {
                to_prep_atom = prep_atoms[i];
            }
        }

        if (from_prep_atom != NULL && to_prep_atom != NULL)
        {
            from_assembly_atom = prep_assembly_atom_map[from_prep_atom];
            to_assembly_atom = prep_assembly_atom_map[to_prep_atom];

            if (from_assembly_atom != NULL && to_assembly_atom != NULL)
            {
                AtomNode* from_assembly_atom_node;
                if (from_assembly_atom->GetNode() == NULL)
                {
                    from_assembly_atom_node = new AtomNode();
                    from_assembly_atom_node->SetAtom(from_assembly_atom);
                    from_assembly_atom->SetNode(from_assembly_atom_node);
                }
                else
                    from_assembly_atom_node = from_assembly_atom->GetNode();

                AtomVector from_assembly_atom_neighbors = from_assembly_atom_node->GetNodeNeighbors();
                if (std::find(from_assembly_atom_neighbors.begin(), from_assembly_atom_neighbors.end(), to_assembly_atom) == from_assembly_atom_neighbors.end())
                    from_assembly_atom_node->AddNodeNeighbor (to_assembly_atom);

                AtomNode* to_assembly_atom_node;
                if (to_assembly_atom->GetNode() == NULL)
                {
                    to_assembly_atom_node = new AtomNode();
                    to_assembly_atom_node->SetAtom(to_assembly_atom);
                    to_assembly_atom->SetNode(to_assembly_atom_node);
                }
                else
                    to_assembly_atom_node = to_assembly_atom->GetNode();

                AtomVector to_assembly_atom_neighbors = to_assembly_atom_node->GetNodeNeighbors();
                if (std::find(to_assembly_atom_neighbors.begin(), to_assembly_atom_neighbors.end(), from_assembly_atom) == to_assembly_atom_neighbors.end())
                    to_assembly_atom_node->AddNodeNeighbor (from_assembly_atom);
            }
        }
    }//for

    AtomVector all_atoms_added = this -> GetAtoms();
    //Set atom node based on prep file parent atom:
    for (unsigned int i = 0; i < all_atoms_added.size(); i++)
    {
        Atom* current_assembly_atom = all_atoms_added[i];
        Atom* assembly_parent = NULL;
        PrepFileSpace::PrepFileAtom* current_assembly_atom_prep_parent = assembly_prep_parent_atom_map[current_assembly_atom];
        if (current_assembly_atom_prep_parent != NULL) //When parent atoms are dummy , this equals NULL, since dummy isn't added.
        {
            assembly_parent = prep_assembly_atom_map[current_assembly_atom_prep_parent];
        }
        AtomNode* assembly_atom_node;
        if (current_assembly_atom->GetNode() == NULL)
        {
            assembly_atom_node = new AtomNode();
            assembly_atom_node->SetAtom(current_assembly_atom);
            current_assembly_atom->SetNode(assembly_atom_node);
        }
        else
            assembly_atom_node = current_assembly_atom->GetNode();

        AtomVector existing_assembly_atom_node_neighbors = assembly_atom_node->GetNodeNeighbors();
        if (std::find(existing_assembly_atom_node_neighbors.begin(), existing_assembly_atom_node_neighbors.end(), assembly_parent) == existing_assembly_atom_node_neighbors.end()
                && assembly_parent != NULL)
            assembly_atom_node->AddNodeNeighbor(assembly_parent);

        if (assembly_parent != NULL)
        {
            AtomNode* assembly_parent_node;
            if (assembly_parent->GetNode() == NULL )
            {
                assembly_parent_node = new AtomNode();
                assembly_parent_node->SetAtom(assembly_parent);
                assembly_parent->SetNode(assembly_parent_node);
            }
            else
                assembly_parent_node = assembly_parent->GetNode();

            AtomVector existing_parent_node_neighbors = assembly_parent_node->GetNodeNeighbors();
            if (std::find(existing_parent_node_neighbors.begin(), existing_parent_node_neighbors.end(), current_assembly_atom) == existing_parent_node_neighbors.end())
                assembly_parent_node-> AddNodeNeighbor(current_assembly_atom);
        }
    }//for
}//BuildResidueFromPrepFileResidue

