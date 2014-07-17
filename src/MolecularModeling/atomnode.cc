#include "../../includes/MolecularModeling/atomnode.hpp"
#include "../../includes/MolecularModeling/atom.hpp"

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
}

