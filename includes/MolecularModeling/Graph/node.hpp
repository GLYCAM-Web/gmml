#ifndef NODE_HPP
#define NODE_HPP

#include <string>
#include <iostream>
#include <vector>

namespace GraphDS
{   class Edge;
    class Node
    {
        public:
            //////////////////////////////////////////////////////////
            //                    TYPE DEFINITION                   //
            //////////////////////////////////////////////////////////
            /*! \typedef
             * List of edges connecting to node
             */
            typedef std::vector<GraphDS::Edge*> EdgeVector;

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

            //////////////////////////////////////////////////////////
            //                       DECONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default deconstructor
              */
            // ~Node();

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
            /*! \fn
              * An accessor function in order to access a node.
              * @return node_value_ attribute of the current object of this node class
              */
            void* GetNodeValue();

            /*! \fn
              * An accessor function in order to acces node type.
              * @return node_type_ attribute of the current object of this node class
              */
              std::string GetNodeType();

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
              * Set the node_value_ attribute of the current node
              * @param node_value The node type attribute of the current object
              */
            void SetNodeValue(void* node_value);

            /*! \fn
              * A mutator function in order to set the node type to the current node object
              * Set the node_type_ attribute of the current node
              * @param node_type The node type attribute of the current object
              */
            void SetNodeType(std::string node_type);

            /*! \fn
              * A mutator function in order to set the id of the current node object
              * Set the node_id_ attribute of the current graph node
              * @param node_id The id attribute of the current object
              */
            void SetNodeId(std::string node_id);

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
            void AddEdge(GraphDS::Edge* edge);

            /*! \fn
              * A mutator function in order to add a tag to the current node object
              * @param tag A tag of type string
              */
            void AddTag(std::string tag);

            /*! \fn
            * A function to generate the ID for a Node.
            */
            std::string GenerateNodeID();
            /////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the node with its adjacent nodes.
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cerr);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////

           void* node_value_;                       /*!< Value of Node of type atom/residue/molecule etc. >*/
           std::string node_type_;                  /*!< Type(Atom/Residue/Molecule etc.) of Node >*/
           std::string node_id_;                    /*!< An identifier for a graph node which a unqiue index >*/
           bool is_visited_ = false;                        /*!< Status of the current node visited/unvisited >*/
           EdgeVector edgeList_;                    /*!< List of nodes adjacent to the current node >*/
           TagsVector tags_;                        /*!< List of tags assigned to the node >*/
          };
}

#endif // NODE_HPP
