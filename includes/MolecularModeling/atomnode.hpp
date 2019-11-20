#ifndef ATOMNODE_HPP
#define ATOMNODE_HPP


#include <string>
#include <iostream>
#include <vector>
#include <cstddef>
#include <iostream>
#include "atom.hpp" // For AtomVector.
#include "../generictypedefs.hpp"

namespace MolecularModeling
{
    class Atom;
    class AtomNode; // Forward declare for the typedef:
    typedef std::vector<AtomNode*>AtomNodeVector;
    class AtomNode
    {
        public:
            //////////////////////////////////////////////////////////
            //                    TYPE DEFINITION                   //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            AtomNode();
            AtomNode(AtomNode* node);
            AtomNode(AtomNode& node);

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
/** \addtogroup Molecular_Data_Structure
               * @{
               */
            /*! \fn
              * An accessor function in order to access to the atom
              * @return atom_ attribute of the current object of this class
              */
            Atom* GetAtom();
            /*! \fn
              * An accessor function in order to access to the node neighbors
              * @return node_neighbors_ attribute of the current object of this class
              */
            AtomVector GetNodeNeighbors();
            /*! \fn
              * An accessor function in order to access to the atom id
              * @return id_ attribute of the current object of this class
              */
            int GetId();

            /// Pattern matching-related accessor functions
            std::string GetElementLabel();
            char GetChiralityLabel();
            AtomVector GetIntraNodeNeighbors();

            /*! \fn                                                  //Added by Ayush on 04/11/2018 for Direction based bonded Atoms in Assembly.
              * An accessor function in order to access whether a Atomnode is visited during Graph Traversal, returns true if a Atomode is visited.
              * @return isVisited_ attribute of the current object of this class
              */
            bool GetIsVisited();
            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////

            /*! \fn
              * A mutator function in order to set the atom of the current object
              * Set the atom_ attribute of the current atom node
              * @param atom_ The atom attribute of the current object
              */
            void SetAtom(Atom* atom);
            /*! \fn
              * A mutator function in order to set the node neighbors of the current object
              * Set the node_neighbors_ attribute of the current atom node
              * @param node_neighbors The node_neighbors attribute of the current object
              */
            void SetNodeNeighbors(AtomVector node_neighbors);
            /*! \fn
              * A function in order to add the node neighbor to the current object
              * Set the node_ attribute of the current atom node
              * @param node neighbor The node neighbor of the current object
              */
            void AddNodeNeighbor(Atom* node_neighbor);
            /*! \fn
              * A mutator function in order to set the id of the current object
              * Set the id_ attribute of the current atom node
              * @param id_ The id attribute of the current object
              */
            void SetId(int id);

            /// Pattern matching-related mutator function
            void RemoveNodeNeighbor(Atom* node_neighbor);
            void SetElementLabel(std::string element_label);
            void SetChiralityLabel(char chirality_label);
            void SetIntraNodeNeighbors(AtomVector intra_node_neighbors);
            void AddIntraNodeNeighbor(Atom* intra_node_neighbor);
            void RemoveIntraNodeNeighbor(Atom* intra_node_neighbor);
            std::string CreateNeighboringLabel(bool excluding_hydrogen = true);
            int GetIntraEdgeDegree();

            /*! \fn                                                            //Added by Ayush on 04/11/2018 for Direction based bonded Atoms in Assembly.
              * A mutator function in order to set true if a Atomnode is visited during Graph Traversal
              * Set the isVisited_ attribute of the current Atomnode object
              * @param isVisited The is_visited_ attribute of the current object
              */
            void SetIsVisited(bool isVisited);
            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the atom node contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cerr);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            Atom* atom_;                        /*!< Pointer back to an atom that this node of graph indicates >*/
            AtomVector node_neighbors_;         /*!< List of all neighbors of this atom in the graph >*/
            int id_;                            /*!< An integer number that indicates the id of a node in the graph >*/
            bool isVisited_;                    /*!< A boolean value to check if the node is visted while traversing atomnodes in the graph >*/ //Added by Ayush on 04/11/2018 for Direction based bonded Atoms in Assembly.

            /// Pattern matching-related attribute
            std::string element_label_;
            char chirality_label_;
            AtomVector intra_node_neighbors_;

    };
}

#endif // ATOMNODE_HPP
