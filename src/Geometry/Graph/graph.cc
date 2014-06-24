
#include "../../../includes/Geometry/Graph/graph.hpp"
#include "../../../includes/Geometry/Graph/edge.hpp"

using namespace std;
using namespace Geometry;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
template<class T>
Graph<T>::Graph() {}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
template<class T>
T Graph<T>::GetVertex()
{
    return vertex_;
}
template<class T>
typename Graph<T>::EdgeVector Graph<T>::GetEdges()
{
    return edges_;
}
template<class T>
gmml::GraphType Graph<T>::GetType()
{
    return type_;
}
template<class T>
int Graph<T>::GetVertexId()
{
    return vertex_id_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////

template<class T>
void Graph<T>::SetVertex(T vertex)
{
    vertex_ = vertex;
}
template<class T>
void Graph<T>::SetEdges(EdgeVector edges)
{
    edges_ = edges;
}
template<class T>
void Graph<T>::SetType(gmml::GraphType type)
{
    type_ = type;
}
template<class T>
void Graph<T>::SetVertexId(int vertex_id)
{
    vertex_id_ = vertex_id;
}

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
template<class T>
void Graph<T>::Print(ostream &out)
{
}







