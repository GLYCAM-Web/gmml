#include "../../../includes/MolecularModeling/Graph/edge.hpp"
#include "../../../includes/MolecularModeling/Graph/node.hpp"
#include "../../../includes/MolecularModeling/Graph/graph.hpp"


//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////

template < class N >
GraphDS::Graph<N>::Graph(){}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////

template < class N >
typename GraphDS::Graph<N>::NodeVector GraphDS::Graph<N>::GetGraphNodeList()
{\
    return graphNodeList_;
}

template < class N >
typename GraphDS::Graph<N>::EdgeVector GraphDS::Graph<N>::GetGraphEdgeList()
{\
    return graphEdgeList_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////

template < class N >
void GraphDS::Graph<N>::SetGraphNodeList(NodeVector nodeList)
{
    graphNodeList_.clear();
    for(typename GraphDS::Graph<N>::NodeVector::iterator it = nodeList.begin(); it != nodeList.end(); it++)
    {
        graphNodeList_.push_back(*it);
    }
}

template < class N >
void GraphDS::Graph<N>::SetGraphEdgeList(EdgeVector edgeList)
{
    graphEdgeList_.clear();
    for(typename GraphDS::Graph<N>::EdgeVector::iterator it = edgeList.begin(); it != edgeList.end(); it++)
    {
        graphEdgeList_.push_back(*it);
    }
}


//////////////////////////////////////////////////////////
//                       FUNCTIONS                      //
//////////////////////////////////////////////////////////

template < class N >
void GraphDS::Graph<N>::AddNewNode(Node<N> *newNode)
{
    graphNodeList_.push_back(newNode);
}

template < class N >
void GraphDS::Graph<N>::AddEdge(Node<N>* firstNode, Node<N>* secondNode)
{
    Edge<N>* newEdge = new Edge<N>(firstNode, secondNode);
    graphEdgeList_.push_back(newEdge);
}
//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
