#ifndef MOLECULARGRAPH_HPP
#define MOLECULARGRAPH_HPP

#include <string>
#include <iostream>
#include <fstream>
#include <vector>

namespace MolecularModeling
{   class Node;
    class MolecularGraph
    {
        public:
            //////////////////////////////////////////////////////////
            //                    TYPE DEFINITION                   //
            //////////////////////////////////////////////////////////
            /*! \typedef
             * List of nodes in a Graph
             */
            typedef std::vector<Node*> NodeVector;

            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            MolecularGraph();

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////


            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////


            //////////////////////////////////////////////////////////
            //                       FUNCTIONS                      //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to add a new node to the Graph
              * Set the node_ attributes of the current graph object
              * @param node_ The node attribute of the current graph object
              */
           void AddNewNode(Node* newNode);

           /*! \fn
             * A function to remove an existing new node from the Graph
             * Remove the node_ attributes of the current graph object
             * @param node_ The node attribute of the current graph object
             */
          void RemoveNode(Node* delNode);

           /*! \fn
             * A function to add an edge between two nodes of the current graph object of this class
             * @param node_ The node attribute of the current graph object
             */
          void AddEdge(Node* firstNode, Node* secondNode);

           /*! \fn
             * A function to find a node in the current graph object using the node id
             * @return node_ attribute of the current graph object of this class
             */

           Node* FindNodeById(long node_id);

            /////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the nodes in a graph.
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cout);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            Node* node_;         /*!< Pointer to a node of type Node >*/
            NodeVector graphNodeList_;         /*!< List of all nodes in the graph >*/
    };
}

#endif // MOLECULARGRAPH_HPP
