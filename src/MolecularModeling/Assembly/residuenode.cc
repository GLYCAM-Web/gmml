#include "../../includes/MolecularModeling/atom.hpp"
#include "../../includes/MolecularModeling/residue.hpp"
#include "../../includes/MolecularModeling/residuenode.hpp"


#include <cstddef>
#include <iostream>

using namespace std;
using namespace MolecularModeling;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
ResidueNode::ResidueNode() : isVisited_(false)
{

    residuenode_neighbors_= ResidueNodeVector();
    cout<<" residuenode_neighbors_ ini"<<endl;
    residuenode_connecting_atoms_= ResidueNodeConnectingAtomVector();
    cout<<" residuenode_connecting_atoms_ ini"<<endl;

}


//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
Residue* ResidueNode::GetResidue()
{
    return residue_;
}

ResidueNode::ResidueNodeVector ResidueNode::GetResidueNodeNeighbors()
{
    return residuenode_neighbors_;
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
void ResidueNode::SetResidue(Residue* residue)
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
    if(residuenode_connecting_atoms_.size()>0){
        cout<<"in if at AddResidueNodeConnectingAtom atom"<<endl;
        for(ResidueNodeConnectingAtomVector::iterator it = residuenode_connecting_atoms_.begin(); it != residuenode_connecting_atoms_.end(); it++)
        {
            Atom* current_connecting_atom = (*it);
             if((current_connecting_atom->GetId()).compare(connecting_atom->GetId())!=0){
                    residuenode_connecting_atoms_.push_back(connecting_atom);
                }
        }
    }else{
         cout<<"in else at AddResidueNodeConnectingAtom atom"<<endl;
        residuenode_connecting_atoms_.push_back(connecting_atom);
        cout<<"after else at AddResidueNodeConnectingAtom atom"<<endl;
    }


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
    for(ResidueNodeVector::iterator it = residuenode_neighbors_.begin(); it != residuenode_neighbors_.end(); it++)
    {
        ResidueNode* current_residuenode = (*it);
        if(current_residuenode->GetId()!=(residuenode_neighbor->GetId())){
             residuenode_neighbors_.push_back(residuenode_neighbor);
        }
    }

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
void ResidueNode::Print(ostream &out)
{
    out << "Residue Name:" << residue_->GetName() << endl;
    out << "Residuenode ID:"<<id_<<endl;
    out <<  "Visit Flag:"<<isVisited_<<endl;
    out <<  "ResidueNode Neighbors: id  residue name"<<endl;
    for(ResidueNodeVector::iterator it = residuenode_neighbors_.begin(); it != residuenode_neighbors_.end(); it++)
    {
         ResidueNode* current_residuenode = (*it);
        out<< current_residuenode->GetId()<<"  ";
        out<< current_residuenode->GetResidue()->GetName()<<endl;
    }


    out <<  "Residue Connecting Atoms: id  atom name"<<endl;
    for(ResidueNodeConnectingAtomVector::iterator it = residuenode_connecting_atoms_.begin(); it != residuenode_connecting_atoms_.end(); it++)
    {
         Atom* connecting_atom = (*it);
         out<< connecting_atom->GetId()<<"  ";
         out<< connecting_atom->GetName()<<endl;
    }
}

