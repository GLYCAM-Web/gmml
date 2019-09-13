#include "../../../includes/MolecularModeling/atom.hpp"
#include "../../../includes/MolecularModeling/residue.hpp"
#include "../../../includes/MolecularModeling/residuenode.hpp"

#include <cstddef>
#include <iostream>

using MolecularModeling::ResidueNode;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
ResidueNode::ResidueNode() : isVisited_(false)
{

    residuenode_neighbors_= ResidueNodeVector();
    residuenode_connecting_atoms_= ResidueNodeConnectingAtomVector();

}


//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
MolecularModeling::Residue* ResidueNode::GetResidue()
{
    return residue_;
}

MolecularModeling::ResidueNodeVector ResidueNode::GetResidueNodeNeighbors()
{
    return residuenode_neighbors_;
}

MolecularModeling::ResidueVector ResidueNode::GetResidueNeighbors()
{
    MolecularModeling::ResidueVector neighbors;
    MolecularModeling::ResidueNodeVector neighboringNodes = this->GetResidueNodeNeighbors();
    for(auto &node : neighboringNodes)
    {
        neighbors.push_back(node->GetResidue());
    }
    return neighbors;
}

ResidueNode::ResidueNodeConnectingAtomVector ResidueNode::GetResidueNodeConnectingAtoms()
{
    return residuenode_connecting_atoms_;
}

int ResidueNode::GetId()
{
    return id_;
}


bool ResidueNode::GetIsVisited()
{
    return isVisited_;
}


//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void ResidueNode::SetResidue(MolecularModeling::Residue* residue)
{
    residue_ = residue;
}

void ResidueNode::SetResidueNodeNeighbors(ResidueNodeVector residuenode_neighbors)
{
    residuenode_neighbors_.clear();
    for(ResidueNodeVector::iterator it = residuenode_neighbors.begin(); it != residuenode_neighbors.end(); it++)
    {
        residuenode_neighbors_.push_back(*it);
    }
}

void ResidueNode::SetResidueNodeConnectingAtoms(ResidueNodeConnectingAtomVector residuenode_connecting_atoms)
{
    residuenode_connecting_atoms_.clear();
    for(ResidueNodeConnectingAtomVector::iterator it = residuenode_connecting_atoms.begin(); it != residuenode_connecting_atoms.end(); it++)
    {
        residuenode_connecting_atoms_.push_back(*it);
    }
}


void ResidueNode::SetId(int id)
{
    id_ = id;
}

void ResidueNode::SetIsVisited(bool isVisited)
{
    isVisited_=isVisited;
}

//////////////////////////////////////////////////////////
//                       FUNCTIONS                      //
//////////////////////////////////////////////////////////


void ResidueNode::AddResidueNodeConnectingAtom(Atom* connecting_atom)
{
    residuenode_connecting_atoms_.push_back(connecting_atom);
}


void ResidueNode::RemoveResidueNodeConnectingAtom(Atom* connecting_atom)

{
    for(ResidueNodeConnectingAtomVector::iterator it = residuenode_connecting_atoms_.begin(); it != residuenode_connecting_atoms_.end(); it++)
    {
        Atom* current_connecting_atom = (*it);
        if((current_connecting_atom->GetId()).compare(connecting_atom->GetId())==0){
        	residuenode_connecting_atoms_.erase(it);
	}
    }
}


void ResidueNode::AddResidueNodeNeighbor(ResidueNode* residuenode_neighbor)
{
    residuenode_neighbors_.push_back(residuenode_neighbor);
}

void ResidueNode::RemoveNodeNeighbor(ResidueNode* residuenode_neighbor)

{
    for(ResidueNodeVector::iterator it = residuenode_neighbors_.begin(); it != residuenode_neighbors_.end(); it++)
    {
        ResidueNode* current_residuenode = (*it);
        if(current_residuenode->GetId()==(residuenode_neighbor->GetId()))
        residuenode_neighbors_.erase(it);
    }
}

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void ResidueNode::Print(std::ostream &out)
{
    out << "Residue Name:" << residue_->GetName() << std::endl;
    out << "Residuenode ID:" << id_ << std::endl;
    out <<  "Visit Flag:" << isVisited_ << std::endl;
    out <<  "ResidueNode Neighbors:\t id  \t residue name" << std::endl;
    for(ResidueNodeVector::iterator it = residuenode_neighbors_.begin(); it != residuenode_neighbors_.end(); it++)
    {
         ResidueNode* current_residuenode = (*it);
        out<<"\t \t \t"<< current_residuenode->GetId()<<"  ";
        out<< "\t"<<current_residuenode->GetResidue()->GetName()<<std::endl;
    }


    out <<  "Residue Connecting Atoms: \t \t id  \t \t atom name"<<std::endl;
    for(ResidueNodeConnectingAtomVector::iterator it = residuenode_connecting_atoms_.begin(); it != residuenode_connecting_atoms_.end(); it++)
    {
         Atom* connecting_atom = (*it);
         out << "\t \t \t \t" << connecting_atom->GetId() << "  ";
         out << "\t" << connecting_atom->GetName() << std::endl;
    }
}
