#include "../../includes/MolecularModeling/atomnode.hpp"
#include "../../includes/MolecularModeling/atom.hpp"
#include "../../includes/MolecularModeling/residue.hpp"
#include "../../includes/utils.hpp"

using MolecularModeling::AtomNode;
//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
AtomNode::AtomNode(): isVisited_(false){}

AtomNode::AtomNode(AtomNode *node)
{
    atom_ = new Atom(new Atom(node->GetAtom()));
    node_neighbors_ = MolecularModeling::AtomVector();
    MolecularModeling::AtomVector node_neighbors = node->GetNodeNeighbors();
    for(MolecularModeling::AtomVector::iterator it = node_neighbors.begin(); it != node_neighbors.end(); it++)
        node_neighbors_.push_back(new Atom(*it));
    id_ = node->GetId();
    isVisited_= false;
}
AtomNode::AtomNode(AtomNode& node)
{
    MolecularModeling::Atom* tempAtom = new MolecularModeling::Atom(*node.GetAtom());
    this->atom_=tempAtom;

   MolecularModeling::AtomVector node_neighbors =node.GetNodeNeighbors();
   for(MolecularModeling::AtomVector::iterator it = node_neighbors.begin(); it != node_neighbors.end(); it++)
         this->node_neighbors_.push_back(new Atom(*it));

    this->id_=node.GetId();
    this->element_label_=node.GetElementLabel();
    this->chirality_label_=node.GetChiralityLabel();

   MolecularModeling::AtomVector intra_node_neighbors=node.GetIntraNodeNeighbors();
   for(MolecularModeling::AtomVector::iterator it = intra_node_neighbors.begin(); it != intra_node_neighbors.end(); it++)
       this->intra_node_neighbors_.push_back(*it);

   isVisited_= false;
}


//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
MolecularModeling::Atom* AtomNode::GetAtom()
{
    return atom_;
}

MolecularModeling::AtomVector AtomNode::GetNodeNeighbors()
{
    return node_neighbors_;
}

int AtomNode::GetId()
{
    return id_;
}

std::string AtomNode::GetElementLabel()
{
    return element_label_;
}

char AtomNode::GetChiralityLabel()
{
    return chirality_label_;
}

MolecularModeling::AtomVector AtomNode::GetIntraNodeNeighbors()
{
    return intra_node_neighbors_;
}

//Added by Ayush on 04/11/2018 for direction based bonded Atoms in Assembly.
bool AtomNode::GetIsVisited()
{
    return isVisited_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void AtomNode::SetAtom(Atom *atom)
{
    atom_ = atom;
}
void AtomNode::SetNodeNeighbors(MolecularModeling::AtomVector node_neighbors)
{
    node_neighbors_.clear();
    for(MolecularModeling::AtomVector::iterator it = node_neighbors.begin(); it != node_neighbors.end(); it++)
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
    MolecularModeling::AtomVector new_node_neighbors = MolecularModeling::AtomVector();
    for(MolecularModeling::AtomVector::iterator it = node_neighbors_.begin(); it != node_neighbors_.end(); it++)
    {
        MolecularModeling::Atom* atom = (*it);
        if(atom->GetId().compare(node_neighbor->GetId()) != 0)
                new_node_neighbors.push_back(atom);
    }
    this->SetNodeNeighbors(new_node_neighbors);
}

void AtomNode::SetElementLabel(std::string element_label)
{
    element_label_ = element_label;
}

void AtomNode::SetChiralityLabel(char chirality_label)
{
    chirality_label_ = chirality_label;
}

void AtomNode::SetIntraNodeNeighbors(MolecularModeling::AtomVector intra_node_neighbors)
{
    intra_node_neighbors_.clear();
    for(MolecularModeling::AtomVector::iterator it = intra_node_neighbors.begin(); it != intra_node_neighbors.end(); it++)
        intra_node_neighbors_.push_back(*it);
}

void AtomNode::AddIntraNodeNeighbor(Atom *intra_node_neighbor)
{
    intra_node_neighbors_.push_back(intra_node_neighbor);
}

void AtomNode::RemoveIntraNodeNeighbor(Atom *intra_node_neighbor)
{
    for(MolecularModeling::AtomVector::iterator it = intra_node_neighbors_.begin(); it != intra_node_neighbors_.end(); it++)
    {
        Atom* atom = (*it);
        if(atom->GetId().compare(intra_node_neighbor->GetId()) == 0)
            intra_node_neighbors_.erase(it);
    }
}

std::string AtomNode::CreateNeighboringLabel(bool excluding_hydrogen)
{
    MolecularModeling::AtomVector inter_neighbors = this->GetNodeNeighbors();
    std::vector<std::string> inter_neighbors_labels = std::vector<std::string>();
    for(MolecularModeling::AtomVector::iterator it = inter_neighbors.begin(); it != inter_neighbors.end(); it++)
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

//Added by Ayush on 04/11/2018 for direction based bonded Atoms in Assembly.
void AtomNode::SetIsVisited(bool isVisited)
{
    isVisited_= isVisited;
}
//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void AtomNode::Print(std::ostream &out)
{
    out << "Element:" << element_label_ << std::endl;
    out << "Atomnode ID:"<<atom_->GetId() <<std::endl;

    for(unsigned int i = 0; i < intra_node_neighbors_.size(); i++)
    {
        out << "\t" << intra_node_neighbors_.at(i)->GetId();
    }
    out << std::endl;
    int number_of_bonds = node_neighbors_.size();
    switch(number_of_bonds)
    {
        case 0:
//            out << atom_->GetResidue()->GetId() << ":" << atom_->GetName() << std::endl;
            out << atom_->GetId() << std::endl;
            break;
        case 1:
//            out << node_neighbors_.at(0)->GetResidue()->GetId() << ":" << node_neighbors_.at(0)->GetName() << " -- "
//                << atom_->GetResidue()->GetId() << ":" << atom_->GetName() << std::endl;
            out << node_neighbors_.at(0)->GetId() << " -- " << atom_->GetId() << std::endl;
            break;
        case 2:
//            out << node_neighbors_.at(0)->GetResidue()->GetId() << ":" << node_neighbors_.at(0)->GetName() << " -- "
//                << atom_->GetResidue()->GetId() << ":" << atom_->GetName() << " -- "
//                << node_neighbors_.at(1)->GetResidue()->GetId() << ":" << node_neighbors_.at(1)->GetName() << std::endl;
            out << node_neighbors_.at(0)->GetId() << " -- " << atom_->GetId() << ":" << " -- "
                << node_neighbors_.at(1)->GetId() << std::endl;
            break;
        case 3:
//            out << "\t\t" << node_neighbors_.at(0)->GetResidue()->GetId() << ":" << node_neighbors_.at(0)->GetName() << std::endl
//                << "\t\t" << "  |  " << std::endl
//                << "\t" << node_neighbors_.at(1)->GetResidue()->GetId() << ":" << node_neighbors_.at(1)->GetName() << " -- "
//                << atom_->GetResidue()->GetId() << ":" << atom_->GetName() << std::endl
//                << "\t\t" << "  |  " << std::endl
//                << "\t\t" <<  node_neighbors_.at(2)->GetResidue()->GetId() << ":" << node_neighbors_.at(2)->GetName() << std::endl;
            out << "\t\t" << node_neighbors_.at(0)->GetId() << std::endl
                << "\t\t" << "  |  " << std::endl
                << "\t" << node_neighbors_.at(1)->GetId() << " -- "
                << atom_->GetId() << std::endl
                << "\t\t" << "  |  " << std::endl
                << "\t\t" <<  node_neighbors_.at(2)->GetId() << std::endl;
            break;
        case 4:
//            out << "\t\t\t" << node_neighbors_.at(0)->GetResidue()->GetId() << ":" << node_neighbors_.at(0)->GetName() << std::endl
//                << "\t\t\t" << "  |  " << std::endl
//                << "\t" << node_neighbors_.at(1)->GetResidue()->GetId() << ":" << node_neighbors_.at(1)->GetName() << " -- "
//                << atom_->GetResidue()->GetId() << ":" << atom_->GetName() << " -- "
//                << node_neighbors_.at(2)->GetResidue()->GetId() << ":" << node_neighbors_.at(2)->GetName() << std::endl
//                << "\t\t\t" << "  |  " << std::endl
//                << "\t\t\t" << node_neighbors_.at(3)->GetResidue()->GetId() << ":" << node_neighbors_.at(3)->GetName() << std::endl;
            out << "\t\t\t" << node_neighbors_.at(0)->GetId() << std::endl
                << "\t\t\t" << "  |  " << std::endl
                << "\t" << node_neighbors_.at(1)->GetId() << " -- "
                << atom_->GetId() << " -- "
                << node_neighbors_.at(2)->GetId() << std::endl
                << "\t\t\t" << "  |  " << std::endl
                << "\t\t\t" << node_neighbors_.at(3)->GetId() << std::endl;
            break;
        case 5:
            break;
    }

}
