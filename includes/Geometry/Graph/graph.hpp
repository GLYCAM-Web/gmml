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
            typedef std::vector<Edge<T>* > EdgeVector;

            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            Graph();

            //////////////////////////////////////////////////////////
            //                         ACCESSOR                     //
            //////////////////////////////////////////////////////////
            /*! \fn
              * An accessor function in order to access to the vertex
              * @return vertex_ attribute of the current object of this class
              */
            T GetVertex();
            /*! \fn
              * An accessor function in order to access to the edges
              * @return edges_ attribute of the current object of this class
              */
            EdgeVector GetEdges();
            /*! \fn
              * An accessor function in order to access to the type
              * @return type_ attribute of the current object of this class
              */
            gmml::GraphType GetType();
            /*! \fn
              * An accessor function in order to access to the vertex
              * @return vertex_ attribute of the current object of this class
              */
            int GetVertexId();

            //////////////////////////////////////////////////////////
            //                          MUTATOR                     //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A mutator function in order to set the vertex of the current object
              * Set the vertex_ attribute of the current graph
              * @param vertex The vertex attribute of the current object
              */
            void SetVertex(T vertex);
            /*! \fn
              * A mutator function in order to set the edges of the current object
              * Set the edges_ attribute of the current graph
              * @param edges The edges attribute of the current object
              */
            void SetEdges(EdgeVector edges);
            /*! \fn
              * A mutator function in order to set the edge of the current object
              * Set the edges_ attribute of the current graph
              * @param edge The edge of the current object
              */
            void AddEdge(Edge<T>* edge);
            /*! \fn
              * A mutator function in order to set the type of the current object
              * Set the type_ attribute of the current graph
              * @param type The type attribute of the current object
              */
            void SetType(gmml::GraphType type);
            /*! \fn
              * A mutator function in order to set the vertex_id of the current object
              * Set the vertex_id_ attribute of the current graph
              * @param vertex_id The vertex_id attribute of the current object
              */
            void SetVertexId(int vertex_id);

            //////////////////////////////////////////////////////////
            //                        FUNCTIONS                     //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                      DISPLAY FUNCTION                //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the graph contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cout);

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
