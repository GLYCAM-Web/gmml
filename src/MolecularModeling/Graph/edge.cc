#include "../../../includes/MolecularModeling/Graph/edge.hpp"
#include "../../../includes/MolecularModeling/Graph/node.hpp"


//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////


GraphDS::Edge::Edge(){}


GraphDS::Edge::Edge(GraphDS::Node* srcNode,GraphDS::Node* dstNode)
{
    this->srcNode_=srcNode;
    this->dstNode_=dstNode;
    srcNode->AddEdge(this);
    dstNode->AddEdge(this);
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////


GraphDS::Node* GraphDS::Edge::GetSourceNode()
{
    return srcNode_;
}


GraphDS::Node* GraphDS::Edge::GetDestinationNode()
{
    return dstNode_;
}


typename GraphDS::Edge::LabelVector GraphDS::Edge::GetEdgeLabels()
{
    return labels_;
}


double GraphDS::Edge::GetEdgeWeight()
{
    return weight_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////


void GraphDS::Edge::SetSourceNode(GraphDS::Node* srcNode)
{
    this->srcNode_= srcNode;
}


void GraphDS::Edge::SetDestinationNode(GraphDS::Node* dstNode)
{
    this->dstNode_= dstNode;
}


void GraphDS::Edge::SetEdgeLabels(GraphDS::Edge::LabelVector labels)
{
    this->labels_.clear();
    for(GraphDS::Edge::LabelVector::iterator it = labels.begin(); it != labels.end(); it++)
    {
        labels_.push_back(*it);
    }
}


void GraphDS::Edge::SetEdgeWeight(double weight)
{
    this->weight_= weight;
}
//////////////////////////////////////////////////////////
//                       FUNCTIONS                      //
//////////////////////////////////////////////////////////


void GraphDS::Edge::AddEdgeLabel(std::string label)
{
        labels_.push_back(label);
}
//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////

void GraphDS::Edge::Print(std::ostream &out)
{

//    std::cout<<"Printing Edge details"<<std::endl;
}
