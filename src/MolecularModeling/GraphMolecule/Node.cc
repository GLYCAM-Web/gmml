#include "../../../includes/MolecularModeling/GraphMolecule/Node.hpp"
#include "../../../includes/MolecularModeling/atom.hpp"

using namespace std;
using namespace MolecularModeling;
//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
Node::Node(){}
//Node::Node(Node* node)
//{
//    atom_=node->GetAtom();
//    node_id_ = node->GetAtom()->GetIndex()
//    is_visited_=false;
//    adjNodeList_= NodeVector();
//}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
Atom* Node::GetNode()
{
    return atom_;
}

unsigned long long Node::GetNodeId()
{
    return node_id_;
}


bool Node::GetIsVisited()
{
    return is_visited_;
}


Node::NodeVector Node::GetadjNodeList()
{
    return adjNodeList_;
}


//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void Node::SetNode(Atom *atom)
{
    atom_ = atom;
}

void Node::SetNodeId(long node_id)
{
    node_id_= node_id;
}


void Node::SetIsVisited(bool is_visited)
{
    is_visited_ = is_visited;
}

void Node::SetadjNodeList(NodeVector adjNodeList)
{
    adjNodeList_.clear();
    for(NodeVector::iterator it = adjNodeList.begin(); it != adjNodeList.end(); it++)
    {
      //  Node* adjNode = *it;
        adjNodeList_.push_back(*it);
    }
}

//////////////////////////////////////////////////////////
//                       FUNCTIONS                      //
//////////////////////////////////////////////////////////

void Node::AddadjNode(Node* node)
{
    adjNodeList_.push_back(node);
}


//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void Node::Print(ostream &out)
{
     for(NodeVector::iterator it = adjNodeList_.begin(); it != adjNodeList_.end(); it++)
     {
            Node* adjNode = *it;
            out<<adjNode->GetNodeId()<<",";
     }

}

