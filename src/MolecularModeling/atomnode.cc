#include "../../includes/MolecularModeling/atomnode.hpp"
#include "../../includes/MolecularModeling/atom.hpp"
#include "../../includes/MolecularModeling/residue.hpp"

using namespace std;
using namespace MolecularModeling;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
AtomNode::AtomNode() {}

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

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void AtomNode::Print(ostream &out)
{
    int number_of_bonds = node_neighbors_.size();
    switch(number_of_bonds)
    {
        case 0:
            out << atom_->GetResidue()->GetId() << ":" << atom_->GetName() << endl;
            break;
        case 1:
            out << node_neighbors_.at(0)->GetResidue()->GetId() << ":" << node_neighbors_.at(0)->GetName() << " -- "
                << atom_->GetResidue()->GetId() << ":" << atom_->GetName() << endl;
            break;
        case 2:
            out << node_neighbors_.at(0)->GetResidue()->GetId() << ":" << node_neighbors_.at(0)->GetName() << " -- "
                << atom_->GetResidue()->GetId() << ":" << atom_->GetName() << " -- "
                << node_neighbors_.at(1)->GetResidue()->GetId() << ":" << node_neighbors_.at(1)->GetName() << endl;
            break;
        case 3:
            out << "\t\t" << node_neighbors_.at(0)->GetResidue()->GetId() << ":" << node_neighbors_.at(0)->GetName() << endl
                << "\t\t" << "  |  " << endl
                << "\t" << node_neighbors_.at(1)->GetResidue()->GetId() << ":" << node_neighbors_.at(1)->GetName() << " -- "
                << atom_->GetResidue()->GetId() << ":" << atom_->GetName() << endl
                << "\t\t" << "  |  " << endl
                << "\t\t" <<  node_neighbors_.at(2)->GetResidue()->GetId() << ":" << node_neighbors_.at(2)->GetName() << endl;
            break;
        case 4:
            out << "\t\t\t" << node_neighbors_.at(0)->GetResidue()->GetId() << ":" << node_neighbors_.at(0)->GetName() << endl
                << "\t\t\t" << "  |  " << endl
                << "\t" << node_neighbors_.at(1)->GetResidue()->GetId() << ":" << node_neighbors_.at(1)->GetName() << " -- "
                << atom_->GetResidue()->GetId() << ":" << atom_->GetName() << " -- "
                << node_neighbors_.at(2)->GetResidue()->GetId() << ":" << node_neighbors_.at(2)->GetName() << endl
                << "\t\t\t" << "  |  " << endl
                << "\t\t\t" << node_neighbors_.at(3)->GetResidue()->GetId() << ":" << node_neighbors_.at(3)->GetName() << endl;
            break;
        case 5:
            break;
    }
}

