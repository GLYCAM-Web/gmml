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
    node_neighbors_ = node_neighbors;
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
            out << atom_->GetResidue()->GetName() << ":" << atom_->GetName() << endl;
            break;
        case 1:
            out << node_neighbors_.at(0)->GetResidue()->GetName() << ":" << node_neighbors_.at(0)->GetName() << "--"
                << atom_->GetResidue()->GetName() << ":" << atom_->GetName() << endl;
            break;
        case 2:
            out << node_neighbors_.at(0)->GetResidue()->GetName() << ":" << node_neighbors_.at(0)->GetName() << "--"
                << atom_->GetResidue()->GetName() << ":" << atom_->GetName() << "--"
                << node_neighbors_.at(1)->GetResidue()->GetName() << ":" << node_neighbors_.at(1)->GetName() << endl;
            break;
        case 3:
            out << "\t" << setw(20) << node_neighbors_.at(0)->GetResidue()->GetName() << ":" << node_neighbors_.at(0)->GetName() << endl
                << "t" << setw(20) << "|" << endl
                << "\t" << setw(20) << atom_->GetResidue()->GetName() << ":" << atom_->GetName() << endl
                << "\t" << setw(10) << right << "/" << setw(10) << left << "\\" << endl
                << "\t" << setw(10) << right << node_neighbors_.at(1)->GetResidue()->GetName() << ":" << node_neighbors_.at(1)->GetName()
                << setw(10) << left << node_neighbors_.at(2)->GetResidue()->GetName() << ":" << node_neighbors_.at(2)->GetName() << endl;
            break;
        case 4:
            out << "\t" << setw(40) << node_neighbors_.at(0)->GetResidue()->GetName() << ":" << node_neighbors_.at(0)->GetName() << endl
                << "t" << setw(40) << "|" << endl
                << "\t" << setw(40) << node_neighbors_.at(1)->GetResidue()->GetName() << ":" << node_neighbors_.at(1)->GetName() << "--"
                << atom_->GetResidue()->GetName() << ":" << atom_->GetName() << "--"
                << node_neighbors_.at(2)->GetResidue()->GetName() << ":" << node_neighbors_.at(2)->GetName() << endl
                << "t" << setw(40) << "|" << endl
                << "\t" << setw(40) << node_neighbors_.at(3)->GetResidue()->GetName() << ":" << node_neighbors_.at(3)->GetName() << endl;
            break;
        case 5:
            break;
    }
}

