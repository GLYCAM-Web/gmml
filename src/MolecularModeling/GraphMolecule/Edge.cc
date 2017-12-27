#include "../../../includes/MolecularModeling/GraphMolecule/Edge.hpp"


using namespace std;
using namespace MolecularModeling;
//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////

Edge::Edge(Node* orgNode,Node* dstNode)
{
    orgNode_= orgNode;
    dstNode_= dstNode;
}


//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
Node* Edge::GetOriginatingNode()
{
    return orgNode_;
}

Node* Edge::GetDestinationNode()
{
    return dstNode_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void Edge::SetOriginatingNode(Node* orgNode)
{
    orgNode_ = orgNode;
}
void Edge::SetDestinationNode(Node* dstNode)
{
    dstNode_ = dstNode;
}
//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void Edge::Print(ostream &out)
{
        out << "------------------------ Displaying the nodes of an edge: --------------------------" << endl;
        out<<"Orginating Node : "<< orgNode_<<endl;
        out<<"Orginating Node : "<< dstNode_<<endl;
}

