#include "../../includes/MolecularModeling/atomnode.hpp"
#include "../../includes/MolecularModeling/atom.hpp"
#include "../../includes/MolecularModeling/residue.hpp"
#include "../../includes/utils.hpp"

using namespace std;
using namespace MolecularModeling;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
AtomNode::AtomNode() {}
AtomNode::AtomNode(AtomNode *node)
{
    atom_ = new Atom(new Atom(node->GetAtom()));
    node_neighbors_ = AtomVector();
    AtomVector node_neighbors = node->GetNodeNeighbors();
    for(AtomVector::iterator it = node_neighbors.begin(); it != node_neighbors.end(); it++)
        node_neighbors_.push_back(new Atom(*it));
    id_ = node->GetId();
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
Atom* AtomNode::GetAtom()
{
    return atom_;
}

AtomNode::AtomVector AtomNode::GetNodeNeighbors()
{
    return node_neighbors_;
}

int AtomNode::GetId()
{
    return id_;
}

string AtomNode::GetElementLabel()
{
    return element_label_;
}

char AtomNode::GetChiralityLabel()
{
    return chirality_label_;
}

AtomNode::AtomVector AtomNode::GetIntraNodeNeighbors()
{
    return intra_node_neighbors_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void AtomNode::SetAtom(Atom *atom)
{
    atom_ = atom;
}
void AtomNode::SetNodeNeighbors(AtomVector node_neighbors)
{
    node_neighbors_.clear();
    for(AtomVector::iterator it = node_neighbors.begin(); it != node_neighbors.end(); it++)
    {
        node_neighbors_.push_back(*it);
    }
}
void AtomNode::AddNodeNeighbor(Atom *node_neighbor)
{
    node_neighbors_.push_back(node_neighbor);
}

void AtomNode::SetId(int id)
{
    id_ = id;
}
void AtomNode::RemoveNodeNeighbor(Atom *node_neighbor)
{    
    AtomVector new_node_neighbors = AtomVector();
    for(AtomVector::iterator it = node_neighbors_.begin(); it != node_neighbors_.end(); it++)
    {
        Atom* atom = (*it);
        if(atom->GetId().compare(node_neighbor->GetId()) != 0)
                new_node_neighbors.push_back(atom);
    }
    this->SetNodeNeighbors(new_node_neighbors);
}

void AtomNode::SetElementLabel(string element_label)
{
    element_label_ = element_label;
}

void AtomNode::SetChiralityLabel(char chirality_label)
{
    chirality_label_ = chirality_label;
}

void AtomNode::SetIntraNodeNeighbors(AtomVector intra_node_neighbors)
{
    intra_node_neighbors_.clear();
    for(AtomVector::iterator it = intra_node_neighbors.begin(); it != intra_node_neighbors.end(); it++)
        intra_node_neighbors_.push_back(*it);
}

void AtomNode::AddIntraNodeNeighbor(Atom *intra_node_neighbor)
{
    intra_node_neighbors_.push_back(intra_node_neighbor);
}

void AtomNode::RemoveIntraNodeNeighbor(Atom *intra_node_neighbor)
{
    for(AtomVector::iterator it = intra_node_neighbors_.begin(); it != intra_node_neighbors_.end(); it++)
    {
        Atom* atom = (*it);
        if(atom->GetId().compare(intra_node_neighbor->GetId()) == 0)
            intra_node_neighbors_.erase(it);
    }
}

string AtomNode::CreateNeighboringLabel(bool excluding_hydrogen)
{
    AtomVector inter_neighbors = this->GetNodeNeighbors();
    vector<string> inter_neighbors_labels = vector<string>();
    for(AtomVector::iterator it = inter_neighbors.begin(); it != inter_neighbors.end(); it++)
    {
        if((*it)->GetNode()->GetElementLabel() != "H" && excluding_hydrogen)
            inter_neighbors_labels.push_back((*it)->GetNode()->GetElementLabel());
        else if((*it)->GetNode()->GetElementLabel() == "H" && !excluding_hydrogen)
            inter_neighbors_labels.push_back((*it)->GetNode()->GetElementLabel());
    }
    sort(inter_neighbors_labels.begin(), inter_neighbors_labels.end());
    return gmml::ConvertVectorString2String(inter_neighbors_labels);
}

int AtomNode::GetIntraEdgeDegree()
{
    return intra_node_neighbors_.size();
}

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void AtomNode::Print(ostream &out)
{
    out << "Element:" << element_label_ << endl;
    out << atom_->GetId() << ": ";
    for(unsigned int i = 0; i < intra_node_neighbors_.size(); i++)
    {
        out << "\t" << intra_node_neighbors_.at(i)->GetId();
    }
    out << endl;
    int number_of_bonds = node_neighbors_.size();
    switch(number_of_bonds)
    {
        case 0:
//            out << atom_->GetResidue()->GetId() << ":" << atom_->GetName() << endl;
            out << atom_->GetId() << endl;
            break;
        case 1:
//            out << node_neighbors_.at(0)->GetResidue()->GetId() << ":" << node_neighbors_.at(0)->GetName() << " -- "
//                << atom_->GetResidue()->GetId() << ":" << atom_->GetName() << endl;
            out << node_neighbors_.at(0)->GetId() << " -- " << atom_->GetId() << endl;
            break;
        case 2:
//            out << node_neighbors_.at(0)->GetResidue()->GetId() << ":" << node_neighbors_.at(0)->GetName() << " -- "
//                << atom_->GetResidue()->GetId() << ":" << atom_->GetName() << " -- "
//                << node_neighbors_.at(1)->GetResidue()->GetId() << ":" << node_neighbors_.at(1)->GetName() << endl;
            out << node_neighbors_.at(0)->GetId() << " -- " << atom_->GetId() << ":" << " -- "
                << node_neighbors_.at(1)->GetId() << endl;
            break;
        case 3:
//            out << "\t\t" << node_neighbors_.at(0)->GetResidue()->GetId() << ":" << node_neighbors_.at(0)->GetName() << endl
//                << "\t\t" << "  |  " << endl
//                << "\t" << node_neighbors_.at(1)->GetResidue()->GetId() << ":" << node_neighbors_.at(1)->GetName() << " -- "
//                << atom_->GetResidue()->GetId() << ":" << atom_->GetName() << endl
//                << "\t\t" << "  |  " << endl
//                << "\t\t" <<  node_neighbors_.at(2)->GetResidue()->GetId() << ":" << node_neighbors_.at(2)->GetName() << endl;
            out << "\t\t" << node_neighbors_.at(0)->GetId() << endl
                << "\t\t" << "  |  " << endl
                << "\t" << node_neighbors_.at(1)->GetId() << " -- "
                << atom_->GetId() << endl
                << "\t\t" << "  |  " << endl
                << "\t\t" <<  node_neighbors_.at(2)->GetId() << endl;
            break;
        case 4:
//            out << "\t\t\t" << node_neighbors_.at(0)->GetResidue()->GetId() << ":" << node_neighbors_.at(0)->GetName() << endl
//                << "\t\t\t" << "  |  " << endl
//                << "\t" << node_neighbors_.at(1)->GetResidue()->GetId() << ":" << node_neighbors_.at(1)->GetName() << " -- "
//                << atom_->GetResidue()->GetId() << ":" << atom_->GetName() << " -- "
//                << node_neighbors_.at(2)->GetResidue()->GetId() << ":" << node_neighbors_.at(2)->GetName() << endl
//                << "\t\t\t" << "  |  " << endl
//                << "\t\t\t" << node_neighbors_.at(3)->GetResidue()->GetId() << ":" << node_neighbors_.at(3)->GetName() << endl;
            out << "\t\t\t" << node_neighbors_.at(0)->GetId() << endl
                << "\t\t\t" << "  |  " << endl
                << "\t" << node_neighbors_.at(1)->GetId() << " -- "
                << atom_->GetId() << " -- "
                << node_neighbors_.at(2)->GetId() << endl
                << "\t\t\t" << "  |  " << endl
                << "\t\t\t" << node_neighbors_.at(3)->GetId() << endl;
            break;
        case 5:
            break;
    }
}

