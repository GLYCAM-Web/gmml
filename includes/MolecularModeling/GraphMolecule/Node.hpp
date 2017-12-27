#ifndef NODE_HPP
#define NODE_HPP

#include <string>
#include <iostream>
#include <fstream>
#include <vector>

namespace MolecularModeling
{   class Atom;
    class Node
    {
        public:
            //////////////////////////////////////////////////////////
            //                    TYPE DEFINITION                   //
            //////////////////////////////////////////////////////////
            /*! \typedef
             * List of nodes adjacent
             */
            typedef std::vector<Node*> NodeVector;

            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            Node();

            /*! \fn
              * Parameterized constructor
              */
         //   Node(Node* node);
           // Node(Atom& atom);
           // Node(Residue* residue);
           // Node(Residue& residue);
           // Node(Molecule* molecule);
           // Node(Molecule& molecule);



            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
            /*! \fn
              * An accessor function in order to access a node.
              * @return atom_ attribute of the current object of this node class
              */
            Atom* GetNode();

            /*! \fn
              * An accessor function in order to access to the graph node id
              * @return id_ attribute of the current object of this class
              */
            unsigned long long GetNodeId();

            /*! \fn
              * An accessor function in order to access whether a node is visited, returns true if a node is visited.
              * @return is_visited_ attribute of the current object of node class
              */
            bool GetIsVisited();

            /*! \fn
              * An accessor function in order to access the nodes adjacent to a node.Returns the list of nodes.
              * @return  adjNodeList_ attribute of the current object of node class
              */
            NodeVector GetadjNodeList();

            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A mutator function in order to set the atom to the current node object
              * Set the atom_ attribute of the current node of type atom
              * @param atom_ The atom attribute of the current object
              */
            void SetNode(Atom* atom);

            /*! \fn
              * A mutator function in order to set the id of the current object
              * Set the id_ attribute of the current graph node
              * @param id_ The id attribute of the current object
              */
            void SetNodeId(long node_id);

            /*! \fn
              * A mutator function in order to set true if a node is visited of the current object
              * Set the is_visited_ attribute of the current node object
              * @param is_visited_ The is_visited_ attribute of the current node object
              */
            void SetIsVisited(bool is_visited);

            /*! \fn
              * A mutator function in order to set an adjacent node to the current node object
              * @param adjNodeList_ The new node to be adjacent to the current node object.This will create an edge.
              */
            void SetadjNodeList(NodeVector adjNodeList);

            //////////////////////////////////////////////////////////
            //                       FUNCTIONS                      //
            //////////////////////////////////////////////////////////
            void AddadjNode(Node* node);

            /////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the node with its adjacent nodes.
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cout);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
           Atom* atom_;                     /*!< Pointer to a node of type Atom >*/
           unsigned long long node_id_;     /*!< An identifier for a graph node which a unqiue index for each atom in an assembly >*/
           bool is_visited_;                /*!< Status of the current node visited/unvisited >*/
           NodeVector adjNodeList_;         /*!< List of nodes adjacent to the current node >*/
         //Edge* edge_;                     /*!< Pointer to a edge for a node of type Edge >*/


    };
}

#endif // NODE_HPP
