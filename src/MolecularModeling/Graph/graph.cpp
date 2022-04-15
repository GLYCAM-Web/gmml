#include "../../../includes/MolecularModeling/Graph/edge.hpp"
#include "../../../includes/MolecularModeling/Graph/node.hpp"
#include "../../../includes/MolecularModeling/Graph/graph.hpp"


//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////


GraphDS::Graph::Graph()
{
    this->graph_id_ = this->GenerateGraphID();
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////


GraphDS::Graph::NodeVector GraphDS::Graph::GetGraphNodeList()
{\
    return graphNodeList_;
}


GraphDS::Graph::EdgeVector GraphDS::Graph::GetGraphEdgeList()
{\
    return graphEdgeList_;
}


std::string GraphDS::Graph::GetGraphName()
{
    return graph_name_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////


void GraphDS::Graph::SetGraphNodeList(NodeVector nodeList)
{
    graphNodeList_.clear();
    for( GraphDS::Graph::NodeVector::iterator it = nodeList.begin(); it != nodeList.end(); it++)
    {
        graphNodeList_.push_back(*it);
    }
}


void GraphDS::Graph::SetGraphEdgeList(EdgeVector edgeList)
{
    graphEdgeList_.clear();
    for( GraphDS::Graph::EdgeVector::iterator it = edgeList.begin(); it != edgeList.end(); it++)
    {
        graphEdgeList_.push_back(*it);
    }
}

void GraphDS::Graph::SetGraphName(std::string graph_name)
{
       graph_name_=graph_name;
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
    Edge* newEdge = new Edge(firstNode, secondNode);
    graphEdgeList_.push_back(newEdge);
    firstNode->AddEdge(newEdge);
    secondNode->AddEdge(newEdge);
}

void GraphDS::Graph::AddEdge(Edge* newEdge)
{
  graphEdgeList_.push_back(newEdge);
}

void GraphDS::Graph::RemoveNode(GraphDS::Node *delNode)
{
    for(GraphDS::Graph::NodeVector::iterator it = graphNodeList_.begin(); it != graphNodeList_.end(); it++)
    {

            GraphDS::Node* current_node = (*it);
            if(current_node->GetNodeId().compare(delNode->GetNodeId()) != 0)
                 graphNodeList_.erase(it);
    }
}

GraphDS::Node* GraphDS::Graph::FindNodeById(std::string node_id)
{
    for(GraphDS::Graph::NodeVector::iterator it = graphNodeList_.begin(); it != graphNodeList_.end(); it++)
    {
            GraphDS::Node* current_node=NULL;

            current_node = (*it);
            if(current_node->GetNodeId().compare(node_id) != 0)
                return current_node;
            else
                return current_node;
    }
}

std::string GraphDS::Graph::GenerateGraphID()
{
    static unsigned int GraphIndex = 0; // static keyword means it is created only once and persists beyond scope of code block.
    GraphIndex++; // makes copy of GraphIndex, increments the real GraphIndex
    return "Graph_"+std::to_string(GraphIndex); // returns the value in the copy
}
//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////


void GraphDS::Graph::Print(std::ostream &out)
{
    out<<"Printing Graph details"<<std::endl;

    out<<"Format: Node -->Connected Nodes"<<std::endl;
    for(GraphDS::Graph::NodeVector::iterator it = graphNodeList_.begin(); it != graphNodeList_.end(); it++)
    {
            GraphDS::Node* current_node=(*it);
            out<< current_node->GetNodeId();

            for( GraphDS::Graph::EdgeVector::iterator it1 = graphEdgeList_.begin(); it1 != graphEdgeList_.end(); it1++)
            {
                GraphDS::Edge *current_edge = (*it1);
                if(current_edge->GetSourceNode() == current_node)
                {
                  out<< "--> " << current_edge->GetDestinationNode()->GetNodeId();
                }
                
            }
            out<<std::endl;

    }
}
