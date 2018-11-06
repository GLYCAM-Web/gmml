#include "../../includes/MolecularModeling/Graph/edge.hpp"
#include "../../includes/MolecularModeling/Graph/node.hpp"
#include "../../includes/MolecularModeling/Graph/graph.hpp"


//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////

GraphDS::Graph::Graph(){}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////

GraphDS::Graph::NodeVector GraphDS::Graph::GetGraphNodeList()
{\
    return this->graphNodeList_;
}

GraphDS::Graph::EdgeVector GraphDS::Graph::GetGraphEdgeList()
{\
    return this->graphEdgeList_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////

void GraphDS::Graph::SetGraphNodeList(NodeVector nodeList)
{
    graphNodeList_.clear();
    for(GraphDS::Graph::NodeVector::iterator it = nodeList.begin(); it != nodeList.end(); it++)
    {
        graphNodeList_.push_back(*it);
    }
}

void GraphDS::Graph::SetGraphEdgeList(EdgeVector edgeList)
{
    graphEdgeList_.clear();
    for(GraphDS::Graph::EdgeVector::iterator it = edgeList.begin(); it != edgeList.end(); it++)
    {
        graphEdgeList_.push_back(*it);
    }
}


//////////////////////////////////////////////////////////
//                       FUNCTIONS                      //
//////////////////////////////////////////////////////////

void GraphDS::Graph::AddNewNode(Node *newNode)
{
    graphNodeList_.push_back(newNode);
}

void GraphDS::Graph::AddEdge(Node* firstNode, Node* secondNode)
{
    Edge *newEdge = new Edge(firstNode, secondNode);
    graphEdgeList_.push_back(newEdge);
}
//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
