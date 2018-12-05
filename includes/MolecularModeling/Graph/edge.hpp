#ifndef EDGE_HPP
#define EDGE_HPP

#include <string>
#include <iostream>
#include <vector>

#include "node.hpp"

namespace GraphDS
{
    template < class N >
    class Edge
    {
        public:

            //////////////////////////////////////////////////////////
            //                    TYPE DEFINITION                   //
            //////////////////////////////////////////////////////////

            /*! \typedef
             * List of tags assigned to a node
             */
            typedef typename std::vector<std::string> LabelVector;

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
            Edge(GraphDS::Node<N>* srcNode,GraphDS::Node<N>* dstNode);

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
            /*! \fn
              * An accessor function in order to access source node of an edge.
              * @return srcNode_ attribute of the current object of the edge class
              */
            GraphDS::Node<N>* GetSourceNode();

            /*! \fn
              * An accessor function in order to access destination node of an edge.
              * @return orgNode_ attribute of the current object of the edge class
              */
            GraphDS::Node<N>* GetDestinationNode();

            /*! \fn
              * An accessor function in order to access labels of the edge.
              * @return labels_ attribute of the current object of the edge class
              */
            LabelVector GetEdgeLabels();

            /*! \fn
              * An accessor function in order to access weight the edge.
              * @return weight_ attribute of the current object of the edge class
              */
            double GetEdgeWeight();


            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A mutator function in order to set the source node of the current edge object
              * Set the srcNode_ attribute of the current edge
              * @param srcNode The originating node attribute of the current edge object
              */
            void SetSourceNode(GraphDS::Node<N>* srcNode);

            /*! \fn
              * A mutator function in order to set the destination node of the current edge object
              * Set the dstNode_ attribute of the current edge
              * @param dstNode The destination node attribute of the current edge object
              */
            void SetDestinationNode(GraphDS::Node<N>* dstNode);

            /*! \fn
              * A mutator function in order to set the label of the current edge object
              * Set the labels_ attribute of the current edge
              * @param label The LabelVector attribute of the current edge object
              */
            void SetEdgeLabels(LabelVector labels);

            /*! \fn
              * A mutator function in order to set the distance/lenght between two node of the current edge object
              * Set the distance_ attribute of the current edge
              * @param distance The destination node attribute of the current edge object
              */
            void SetEdgeWeight(double weight);

            //////////////////////////////////////////////////////////
            //                       FUNCTIONS                      //
            //////////////////////////////////////////////////////////

            /*! \fn
              * A function to add an edge label to the LabelVector
              * Add a new label to the labels_ attribute
              * @param label A label of type string
              */
             void AddEdgeLabel(std::string label);
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
           GraphDS::Node<N>* srcNode_;         /*!< Pointer to a originating node of type Node class >*/
           GraphDS::Node<N>* dstNode_;         /*!< Pointer to a destination node of type Node class >*/
           LabelVector labels_;                /*!< Label assigned to the edge between two node >*/
           double weight_;                     /*!<  Weight between between two node/ Length of Edge >*/
    };
}
#include "../../../src/MolecularModeling/Graph/edge.cc"
#endif // EDGE_HPP
