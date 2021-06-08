#ifndef GRAPH_HPP
#define GRAPH_HPP

#include <string>
#include <iostream>
#include <fstream>
#include <vector>

#include "node.hpp"
#include "edge.hpp"

namespace GraphDS
{
class Graph
    {
        public:
            //////////////////////////////////////////////////////////
            //                    TYPE DEFINITION                   //
            //////////////////////////////////////////////////////////
            /*! \typedef
             * List of nodes in a Graph
             */
            typedef std::vector<GraphDS::Node*> NodeVector;
            typedef std::vector<GraphDS::Edge*> EdgeVector;

            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
             Graph();

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////

             /*! \fn
               * An accessor function in order to access the nodes of the graph
               * @return  graphNodeList_ attribute of the current object of Graph class
               */
             NodeVector GetGraphNodeList();

             /*! \fn
               * An accessor function in order to access the edges of the graph
               * @return  graphEdgeList_ attribute of the current object of Graph class
               */
             EdgeVector GetGraphEdgeList();

             /*! \fn
               * An accessor function in order to access the name of the graph
               * @return  graph_name_ attribute of the current object of Graph class
               */
             std::string GetGraphName();

            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////

             /*! \fn
               * A mutator function in order to set nodes to the current Graph
               * @param graphNodeList_ List of nodes of the current Graph class
               */
             void SetGraphNodeList(NodeVector nodeList);

             /*! \fn
               * A mutator function in order to set edges to the current Graph
               * @param graphEdgeList_ List of outgoing/incoming edges of the current Graph class
               */
             void SetGraphEdgeList(EdgeVector edgeList);

             /*! \fn
               * A mutator function in order to set the Graph name
               * @param graph_name string graph name from the user
               */
             void SetGraphName(std::string graph_name);

            //////////////////////////////////////////////////////////
            //                       FUNCTIONS                      //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to add a new node to the Graph
              * Set the node_ attributes of the current graph object
              * @param node_ The node attribute of the current graph object
              */
           void AddNewNode(GraphDS::Node* newNode);

           /*! \fn
             * A function to remove an existing new node from the Graph
             * Remove the node_ attributes of the current graph object
             * @param node_ The node attribute of the current graph object
             */
          void RemoveNode(GraphDS::Node* delNode);

           /*! \fn
             * A function to add an edge between two nodes of the current graph object of this class
             * @param node_ The node attribute of the current graph object
             */
          void AddEdge(GraphDS::Node* firstNode, GraphDS::Node* secondNode);
          
          /*! \fn
             * A function to add an edge to the current graph object of this class
             * @param newEdge The edge to add to the current graph object
             */
          void AddEdge(Edge* newEdge);
          
           /*! \fn
             * A function to find a node in the current graph object using the node id
             * @return node_ attribute of the current graph object of this class
             */

           GraphDS::Node* FindNodeById(std::string node_id);

           /*! \fn
           * A function to generate the ID for a Graph.
           */
           std::string GenerateGraphID();

            /////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the nodes in a graph.
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cerr);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            NodeVector graphNodeList_;         /*!< List of all nodes in the graph >*/
            EdgeVector graphEdgeList_;         /*!< List of all edges in the graph >*/
            std::string graph_id_;             /*!< An identifier for a graph which is a unqiue ID >*/
            std::string graph_name_;           /*!< Name of the Graph assigned >*/
    };
}
#endif // GRAPH_HPP
