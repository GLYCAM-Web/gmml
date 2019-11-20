#ifndef RESIDUENODE_HPP
#define RESIDUENODE_HPP


#include <string>
#include <iostream>
#include <vector>

namespace MolecularModeling
{
    class Residue;
    typedef std::vector<Residue*> ResidueVector;
    class Atom;
    class ResidueNode; // Forward declare for the typedef:
    typedef std::vector<ResidueNode*> ResidueNodeVector;
    class ResidueNode
    {
        public:
            //////////////////////////////////////////////////////////
            //                    TYPE DEFINITION                   //
            //////////////////////////////////////////////////////////
            typedef std::vector<Atom*> ResidueNodeConnectingAtomVector;

            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            ResidueNode();
            //ResidueNode(ResidueNode* node);
            // ResidueNode(ResidueNode& node);

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
            /*! \fn
              * An accessor function in order to access to the residue
              * @return residue_ attribute of the current object of this class
              */
            Residue* GetResidue();

            /*! \fn
              * An accessor function in order to access to the residue node neighbors
              * @return residuenode_neighbors_ attribute of the current object of this class
              */
            ResidueNodeVector GetResidueNodeNeighbors();

            /*! \fn
              * An accessor function in order to access to the residue neighbors
              * @return A ResidueVector containing residue pointers to each Residue that is connected to this residuenode
              */
            ResidueVector GetResidueNeighbors();

            /*! \fn
              * An accessor function in order to access to the connecting atoms of residuenode
              * @return residuenode_connecting_atoms attribute of the current object of this class
              */
            ResidueNodeConnectingAtomVector GetResidueNodeConnectingAtoms();

            /*! \fn
              * An accessor function in order to access to the residuenode id
              * @return id_ attribute of the current object of this class
              */
            int GetId();

            /*! \fn
              * An accessor function in order to access whether a residuenode is visited, returns true if a residuenode is visited.
              * @return isVisited_ attribute of the current object of this class
              */
            bool GetIsVisited();

            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A mutator function in order to set the residue of the current object
              * Set the residue_ attribute of the current residue node
              * @param residue The residue attribute of the current object
              */
            void SetResidue(Residue* residue);

            /*! \fn
              * A mutator function in order to set the residue node neighbors of the current object
              * Set the residuenode_neighbors_ attribute of the current residue node
              * @param residuenode_neighbors The node_neighbors attribute of the current object
              */
            void SetResidueNodeNeighbors(ResidueNodeVector residuenode_neighbors);

            /*! \fn
              * A mutator function in order to set the connecting atoms to another residuenode of the current object
              * Set the residuenode_connecting_atoms_ attribute of the current residue node
              * @param residuenode_connecting_atoms The connecting atoms of residuenode attribute of the current object
              */
            void SetResidueNodeConnectingAtoms(ResidueNodeConnectingAtomVector residuenode_connecting_atoms);

            /*! \fn
              * A mutator function in order to set the id of the current object
              * Set the id_ attribute of the current residuenode
              * @param id The id attribute of the current object
              */
            void SetId(int id);

            /*! \fn
              * A mutator function in order to set true if a residue node is visited of the current object
              * Set the isVisited_ attribute of the current residuenode object
              * @param isVisited The is_visited_ attribute of the current  object
              */
            void SetIsVisited(bool isVisited);

            //////////////////////////////////////////////////////////
            //                       FUNCTIONS                      //
            //////////////////////////////////////////////////////////

            /*! \fn
              * A function in order to add the connecting atom of residue node to the current object
              * @param connecting_atom The atom which connects a residuenode to another residuenode of the current object
              */
            void AddResidueNodeConnectingAtom(Atom* connecting_atom);

            /*! \fn
              * A function in order to remove the connecting atom of residue node to the current object
              * @param connecting_atom The atom which connects a residuenode to another residuenode of the current object
              */
            void RemoveResidueNodeConnectingAtom(Atom* connecting_atom);

            /*! \fn
              * A function in order to add the residuenode neighbor to the current object
              * @param residuenode_neighbor The residuenode which neighbor of the current object
              */
            void AddResidueNodeNeighbor(ResidueNode* residuenode_neighbor);

            /*! \fn
              * A function in order to remove the residuenode, a neighbor to the current object
              *  @param id The residuenode_neighbor of the residuenode which is neighbor of the current object
              */
            void RemoveNodeNeighbor(ResidueNode* residuenode_neighbor);

            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the residuenode contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cerr);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            Residue* residue_;                                              /*!< Pointer back to a residue that this node of graph indicates >*/
            ResidueNodeVector residuenode_neighbors_;                       /*!< List of all neighbors of this residuenode in the graph >*/
            ResidueNodeConnectingAtomVector residuenode_connecting_atoms_;  /*!< List of all atoms connecting to atoms in another residuenode in the graph>*/
            int id_;                                                        /*!< An integer number that indicates the id of a node in the graph >*/
            bool isVisited_;                                                /*!< A boolean value to check if the node is visted while traversing residuenodes in the graph >*/

    };
}

#endif // RESIDUENODE_HPP
