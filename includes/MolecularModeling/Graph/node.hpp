#ifndef NODE_HPP
#define NODE_HPP

#include <string>
#include <iostream>
#include <vector>

namespace GraphDS
{   template<class N> class Edge;
    template<class N> class Node
    {
        public:
            //////////////////////////////////////////////////////////
            //                    TYPE DEFINITION                   //
            //////////////////////////////////////////////////////////
            /*! \typedef
             * List of edges connecting to node
             */
            typedef typename std::vector<Edge<N>*> EdgeVector;

            /*! \typedef
             * List of tags assigned to a node
             */
            typedef typename std::vector<std::string> TagsVector;
            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            Node();

            //////////////////////////////////////////////////////////
            //                       DECONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default deconstructor
              */
            ~Node();

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
            /*! \fn
              * An accessor function in order to access a node.
              * @return node_atom_ attribute of the current object of this node class
              */
            N* GetNode();

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
              * A mutator function in order to set the node of any type to the current node object
              * Set the node_ attribute of the current node
              * @param node The node type attribute of the current object
              */
            void SetNode(N* node);

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
            void AddEdge(Edge<N>* edge);

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

           N* node_;                                 /*!< A template node of type atom/residue/molecule etc. >*/
           std::string node_id_;                    /*!< An identifier for a graph node which a unqiue index >*/
           bool is_visited_;                        /*!< Status of the current node visited/unvisited >*/
           EdgeVector edgeList_;                    /*!< List of nodes adjacent to the current node >*/
           TagsVector tags_;                        /*!< List of tags assigned to the node >*/
          };
}

#include "../../../src/MolecularModeling/Graph/node.cc"
#endif // NODE_HPP
