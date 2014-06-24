
#include "../../../includes/Geometry/Graph/edge.hpp"

using namespace std;
using namespace Geometry;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
template<class T>
Edge<T>::Edge() {}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
template<class T>
Graph<T>* Edge<T>::GetDst()
{
    return dst_;
}
template<class T>
EdgeAttribute* Edge<T>::GetProperties()
{
    return properties_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
template<class T>
void Edge<T>::SetDst(Graph<T> *dst)
{
    dst_ = dst;
}
template<class T>
void Edge<T>::SetProperties(EdgeAttribute *properties)
{
    properties_ = properties;
}

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
template<class T>
void Edge<T>::Print(ostream &out)
{
}








