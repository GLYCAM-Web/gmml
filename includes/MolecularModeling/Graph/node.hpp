#ifndef NODE_HPP
#define NODE_HPP

#include <string>
#include <iostream>
#include <vector>
#include "../../../includes/MolecularModeling/atom.hpp"
#include "../../../includes/MolecularModeling/residue.hpp"
#include "../../../includes/MolecularModeling/molecule.hpp"

namespace GraphDS
{
    class Edge;
    class Node
    {
        public:
            //////////////////////////////////////////////////////////
            //                    TYPE DEFINITION                   //
            //////////////////////////////////////////////////////////
            /*! \typedef
             * List of edges connecting to node
             */
            typedef std::vector<Edge*> EdgeVector;

            /*! \typedef
             * List of tags assigned to a node
             */
            typedef std::vector<std::string> TagsVector;
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
            Node(Node* atom);
            Node(MolecularModeling::Atom& atom);
            Node(MolecularModeling::Residue* residue);
            Node(MolecularModeling::Residue& residue);
            Node(MolecularModeling::Molecule* molecule);
            Node(MolecularModeling::Molecule& molecule);


            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
            /*! \fn
              * An accessor function in order to access a node.
              * @return node_atom_ attribute of the current object of this node class
              */
            MolecularModeling::Atom* GetAtomNode();

            /*! \fn
              * An accessor function in order to access a node.
              * @return node_residue_ attribute of the current object of this node class
              */
            MolecularModeling::Residue* GetResidueNode();

            /*! \fn
              * An accessor function in order to access a node.
              * @return node_molecule_ attribute of the current object of this node class
              */
            MolecularModeling::Molecule* GetMoleculeNode();

            /*! \fn
              * An accessor function in order to access to the graph node id
              * @return id_ attribute of the current object of this class
              */
             std::string GetNodeId();

            /*! \fn
              * An accessor function in order to access whether a node is visited, returns true if a node is visited.
              * @return is_visited_ attribute of the current object of node class
              */
            bool GetIsVisited();

            /*! \fn
              * An accessor function in order to access the edges of the node
              * @return  edgeList_ attribute of the current object of node class
              */
            EdgeVector GetEdgeList();


            /*! \fn
              * An accessor function in order to access the Tags of the current node
              * @return  tags_ attribute of the current object of node class
              */
            TagsVector GetTagList();

            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A mutator function in order to set the atom to the current node object
              * Set the node_atom_ attribute of the current node of type atom
              * @param node_atom The atom attribute of the current object
              */
            void SetAtomNode(MolecularModeling::Atom* node_atom);

            /*! \fn
              * A mutator function in order to set the residue to the current node object
              * Set the node_residue_ attribute of the current node of type residue
              * @param node_residue The residue attribute of the current object
              */
            void SetResidueNode(MolecularModeling::Residue* node_residue);

            /*! \fn
              * A mutator function in order to set the molecule to the current node object
              * Set the node_molecule_ attribute of the current node of type molecule
              * @param node_molecule The molecule attribute of the current object
              */
            void SetMoleculeNode(MolecularModeling::Molecule* node_molecule);

            /*! \fn
              * A mutator function in order to set the id of the current node object
              * Set the node_id_ attribute of the current graph node
              * @param node_id The id attribute of the current object
              */
            void SetNodeId(long node_id);

            /*! \fn
              * A mutator function in order to set true if a node is visited of the current object
              * Set the is_visited_ attribute of the current node object
              * @param is_visited The is_visited_ attribute of the current node object
              */
            void SetIsVisited(bool is_visited);

            /*! \fn
              * A mutator function in order to set edges to the current node object
              * @param edgeList List of outgoing/incoming edges to the current node object
              */
            void SetEdgeList(EdgeVector edgeList);

            /*! \fn
              * A mutator function in order to set the tags to the current node object
              * @param tags The list of tags assigned to the current node object
              */
            void SetTagList(TagsVector tags);

            //////////////////////////////////////////////////////////
            //                       FUNCTIONS                      //
            //////////////////////////////////////////////////////////

            /*! \fn
              * A mutator function to add an edge incoming/outgoing to the current node object
              * @param edge The edge object of class Edge
              */
            void AddEdge(Edge* edge);

            /*! \fn
              * A mutator function in order to add a tag to the current node object
              * @param tag A tag of type string
              */
            void AddTag(std::string tag);
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
           MolecularModeling::Atom* node_atom_;                        /*!< Pointer to a node of type Atom >*/
           MolecularModeling::Residue* node_residue_;                  /*!< Pointer to a node of type Residue >*/
           MolecularModeling::Molecule* node_molecule_;                /*!< Pointer to a node of type Molecule >*/
           std::string node_id_;                    /*!< An identifier for a graph node which a unqiue index for each atom in an assembly >*/
           bool is_visited_;                        /*!< Status of the current node visited/unvisited >*/
           EdgeVector edgeList_;                    /*!< List of nodes adjacent to the current node >*/
           TagsVector tags_;                        /*!< List of tags assigned to the node >*/
          };
}


#endif // NODE_HPP
