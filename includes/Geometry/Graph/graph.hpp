#ifndef GRAPH_HPP
#define GRAPH_HPP

#include <string>
#include <iostream>
#include <vector>

#include "../../common.hpp"

namespace Geometry
{
    template <class T>
    class Edge;
    template <class T>
    class Graph
    {
        public:
            //////////////////////////////////////////////////////////
            //                    TYPE DEFINITION                   //
            //////////////////////////////////////////////////////////
            typedef std::vector<Edge<T >* > EdgeVector;

            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            Graph();

        private:
            //////////////////////////////////////////////////////////
            //                        ATTRIBUTES                    //
            //////////////////////////////////////////////////////////
            T vertex_;
            EdgeVector edges_;
            gmml::GraphType type_;
            int vertex_id_;

    };
}

#endif // GRAPH_HPP
