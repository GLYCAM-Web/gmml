#ifndef EDGE_HPP
#define EDGE_HPP

#include <string>
#include <iostream>
#include <fstream>
#include <vector>

namespace MolecularModeling
{   class Node;
    class Edge
    {
        public:
            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            Edge();

            /*! \fn
              * Parameterized constructor
              */
            Edge(Node* orgNode,Node* dstNode);

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
            /*! \fn
              * An accessor function in order to access originating node of an edge.
              * @return orgNode_ attribute of the current object of this edge class
              */
            Node* GetOriginatingNode();

            /*! \fn
              * An accessor function in order to access destination node of an edge.
              * @return orgNode_ attribute of the current object of this edge class
              */
            Node* GetDestinationNode();

            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A mutator function in order to set the origination node of the current edge object
              * Set the orgNode_ attribute of the current edge
              * @param orgNode_ The originating node attribute of the current edge object
              */
            void SetOriginatingNode(Node* orgNode);

            /*! \fn
              * A mutator function in order to set the destination node of the current edge object
              * Set the dstNode_ attribute of the current edge
              * @param dstNode_ The destination node attribute of the current edge object
              */
            void SetDestinationNode(Node* dstNode);

            /////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the edges in a graph.
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cout);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
           Node* orgNode_;                     /*!< Pointer to a originating node of type Node class >*/
           Node* dstNode_;                     /*!< Pointer to a destination node of type Node class >*/
    };
}

#endif // EDGE_HPP
