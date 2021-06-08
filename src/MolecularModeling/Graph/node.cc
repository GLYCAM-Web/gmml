#include "../../../includes/MolecularModeling/Graph/node.hpp"

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////


GraphDS::Node::Node() {

    // this->is_visited_=false;
    // this->node_id_=this->GenerateNodeID();
}

//////////////////////////////////////////////////////////
//                       DECONSTRUCTOR                    //
//////////////////////////////////////////////////////////


// GraphDS::Node::~Node(){}
//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////


void* GraphDS::Node::GetNodeValue()
{
    return node_value_;
}


std::string GraphDS::Node::GetNodeType()
{
    return node_type_;
}


std::string GraphDS::Node::GetNodeId()
{
    return node_id_;
}


bool GraphDS::Node::GetIsVisited()
{
    return is_visited_;
}


 GraphDS::Node::EdgeVector GraphDS::Node::GetEdgeList()
{
    return edgeList_;
}


 GraphDS::Node::TagsVector GraphDS::Node::GetTagList()
{
    return tags_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////


void GraphDS::Node::SetNodeValue(void* node_value)
{
  // std::cout << node_value << "\n";
    this->node_value_= node_value;
}

void GraphDS::Node::SetNodeType(std::string node_type)
{
    this->node_type_= node_type;
}


void GraphDS::Node::SetNodeId(std::string node_id)
{
    this->node_id_= node_id;
}


void GraphDS::Node::SetIsVisited(bool is_visited)
{
    this->is_visited_=is_visited;
}


void GraphDS::Node::SetEdgeList(EdgeVector edgeList)
{
    this->edgeList_.clear();
    for(typename GraphDS::Node::EdgeVector::iterator it = edgeList.begin(); it != edgeList.end(); it++) {
        this->edgeList_.push_back(*it);
    }
}


void GraphDS::Node::SetTagList(TagsVector tags)
{
    this->tags_.clear();
    for(typename GraphDS::Node::TagsVector::iterator it = tags.begin(); it != tags.end(); it++) {
        this->tags_.push_back(*it);
    }
}

//////////////////////////////////////////////////////////
//                       FUNCTIONS                      //
//////////////////////////////////////////////////////////


void GraphDS::Node::AddEdge(Edge *edge)
{
    edgeList_.push_back(edge);
}


void GraphDS::Node::AddTag(std::string tag)
{
    tags_.push_back(tag);
}

std::string GraphDS::Node::GenerateNodeID()
{
    static unsigned long long NodeIndex = 0; // static keyword means it is created only once and persists beyond scope of code block.
    NodeIndex++; // makes copy of NodeIndex, increments the real NodeIndex
    return "Node_"+std::to_string(NodeIndex); // returns the value in the copy
}
//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////


void GraphDS::Node::Print(std::ostream &out)
{
//    std::cout<<"Printing Node details"<<std::endl;
//    std::cout<<"Node ID:"<<node_id_<<std::endl;
}
