#include "../../../includes/MolecularModeling/Graph/edge.hpp"
#include "../../../includes/MolecularModeling/Graph/node.hpp"


//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////

template < class N >
GraphDS::Edge<N>::Edge(){}

template < class N >
GraphDS::Edge<N>::Edge(GraphDS::Node<N>* srcNode,GraphDS::Node<N>* dstNode)
{
    this->srcNode_=srcNode;
    this->dstNode_=dstNode;
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////

template < class N >
GraphDS::Node<N>* GraphDS::Edge<N>::GetSourceNode()
{
    return srcNode_;
}

template < class N >
GraphDS::Node<N>* GraphDS::Edge<N>::GetDestinationNode()
{
    return dstNode_;
}

template < class N >
typename GraphDS::Edge<N>::LabelVector GraphDS::Edge<N>::GetEdgeLabels()
{
    return labels_;
}

template < class N >
double GraphDS::Edge<N>::GetEdgeWeight()
{
    return weight_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////

template < class N >
void GraphDS::Edge<N>::SetSourceNode(GraphDS::Node<N>* srcNode)
{
    this->srcNode_= srcNode;
}

template < class N >
void GraphDS::Edge<N>::SetDestinationNode(GraphDS::Node<N>* dstNode)
{
    this->dstNode_= dstNode;
}

template < class N >
void GraphDS::Edge<N>::SetEdgeLabels(GraphDS::Edge<N>::LabelVector labels)
{
    this->labels_.clear();
    for(GraphDS::Edge<N>::LabelVector::iterator it = labels.begin(); it != labels.end(); it++)
    {
        labels_.push_back(*it);
    }
}

template < class N >
void GraphDS::Edge<N>::SetEdgeWeight(double weight)
{
    this->weight_= weight;
}
//////////////////////////////////////////////////////////
//                       FUNCTIONS                      //
//////////////////////////////////////////////////////////

template < class N >
void GraphDS::Edge<N>::AddEdgeLabel(std::string label)
{
        labels_.push_back(label);
}
//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////


