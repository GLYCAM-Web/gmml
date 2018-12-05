#include "../../../includes/MolecularModeling/Graph/node.hpp"

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////

template<class N>
GraphDS::Node<N>::Node(): is_visited_(false){}

//////////////////////////////////////////////////////////
//                       DECONSTRUCTOR                    //
//////////////////////////////////////////////////////////

template<class N>
GraphDS::Node<N>::~Node(){}
//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////

template<class N>
N* GraphDS::Node<N>::GetNode()
{
    return node_;
}

template<class N>
std::string GraphDS::Node<N>::GetNodeId()
{
    return node_id_;
}

template<class N>
bool GraphDS::Node<N>::GetIsVisited()
{
    return is_visited_;
}

template<class N>
 typename GraphDS::Node<N>::EdgeVector GraphDS::Node<N>::GetEdgeList()
{
    return edgeList_;
}

template<class N>
typename GraphDS::Node<N>::TagsVector GraphDS::Node<N>::GetTagList()
{
    return tags_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////

template<class N>
void GraphDS::Node<N>::SetNode(N* node)
{
    this->node_= node;
}

template<class N>
void GraphDS::Node<N>::SetNodeId(long node_id)
{
    this->node_id_= node_id;
}

template<class N>
void GraphDS::Node<N>::SetIsVisited(bool is_visited)
{
    this->is_visited_=is_visited;
}

template<class N>
void GraphDS::Node<N>::SetEdgeList(EdgeVector edgeList)
{
    this->edgeList_.clear();
    for(typename GraphDS::Node<N>::EdgeVector::iterator it = edgeList.begin(); it != edgeList.end(); it++) {
        this->edgeList_.push_back(*it);
    }
}

template<class N>
void GraphDS::Node<N>::SetTagList(TagsVector tags)
{
    this->tags_.clear();
    for(typename GraphDS::Node<N>::TagsVector::iterator it = tags.begin(); it != tags.end(); it++) {
        this->tags_.push_back(*it);
    }
}

//////////////////////////////////////////////////////////
//                       FUNCTIONS                      //
//////////////////////////////////////////////////////////

template<class N>
void GraphDS::Node<N>::AddEdge(Edge<N> *edge)
{
    edgeList_.push_back(edge);
}

template<class N>
void GraphDS::Node<N>::AddTag(std::string tag)
{
    tags_.push_back(tag);
}
//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////

template<class N>
void GraphDS::Node<N>::Print(std::ostream &out)
{}
