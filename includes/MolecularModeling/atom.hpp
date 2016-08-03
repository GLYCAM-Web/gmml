#ifndef ATOM_HPP
#define ATOM_HPP

#include <string>
#include <iostream>
#include <vector>

#include "../GeometryTopology/coordinate.hpp"
#include "moleculardynamicatom.hpp"
#include "quantommechanicatom.hpp"
#include "dockingatom.hpp"

namespace MolecularModeling
{
    class Residue;
    class AtomNode;
    class Atom : public MolecularDynamicAtom, public QuantomMechanicAtom, public DockingAtom
    {
        public:
            //////////////////////////////////////////////////////////
            //                    TYPE DEFINITION                   //
            //////////////////////////////////////////////////////////
            /*! \typedef
              * List of coordinates
              */
            typedef std::vector<GeometryTopology::Coordinate*> CoordinateVector;
            typedef std::vector<Atom*> AtomVector;

            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            Atom();
            Atom(Residue* residue, std::string name, CoordinateVector coordinates);
            Atom(Atom* atom);

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
            /*! \fn
              * An accessor function in order to access to the residue
              * @return residue_ attribute of the current object of this class
              */
            Residue* GetResidue();
            /*! \fn
              * An accessor function in order to access to the name
              * @return name_ attribute of the current object of this class
              */
            std::string GetName();
            /*! \fn
              * An accessor function in order to access to the coordinates
              * @return coordinates_ attribute of the current object of this class
              */
            CoordinateVector GetCoordinates();
            /*! \fn
              * An accessor function in order to access to the chemical_type
              * @return chemical_type_ attribute of the current object of this class
              */
            std::string GetChemicalType();
            /*! \fn
              * An accessor function in order to access to the description
              * @return description_ attribute of the current object of this class
              */
            std::string GetDescription();
            /*! \fn
              * An accessor function in order to access to the element symbol
              * @return element_symbol_ attribute of the current object of this class
              */
            std::string GetElementSymbol();
            /*! \fn
              * An accessor function in order to access to the node
              * @return node_ attribute of the current object of this class
              */
            AtomNode* GetNode();
            /*! \fn
              * An accessor function in order to access to the id
              * @return id_ attribute of the current object of this class
              */
            std::string GetId();
            /*! \fn
              * An accessor function in order to access to the is_ring_ attribute of the current object
              * @return is_ring_ attribute of the current object of this class
              */
            bool GetIsRing();
            /*! \fn
              * An accessor function in order to access to the index
              * @return index_ attribute of the current object of this class
              */
            unsigned long long GetIndex();

            //////////////////////////////////////////////////////////
            //                       FUNCTIONS                      //
            //////////////////////////////////////////////////////////
            void FindConnectedAtoms(AtomVector &visitedAtoms);
            double GetDistanceToAtom(Atom *otherAtom);
            unsigned long long generateAtomIndex();

            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A mutator function in order to set the residue of the current object
              * Set the residue_ attribute of the current atom
              * @param residue The residue attribute of the current object
              */
            void SetResidue(Residue* residue);
            /*! \fn
              * A mutator function in order to set the name of the current object
              * Set the name_ attribute of the current atom
              * @param name The name attribute of the current object
              */
            void SetName(std::string name);
            /*! \fn
              * A mutator function in order to set the coordinates of the current object
              * Set the coordinates_ attribute of the current atom
              * @param coordinates The coordinates attribute of the current object
              */
            void SetCoordinates(CoordinateVector coordinates);
            /*! \fn
              * A function in order to add the coordinate to the current object
              * Set the coordinates_ attribute of the current atom
              * @param coordinate The coordinate of the current object
              */
            void AddCoordinate(GeometryTopology::Coordinate* coordinate);
            /*! \fn
              * A mutator function in order to set the chemical type of the current object
              * Set the chemical_type_ attribute of the current atom
              * @param chemical_type The chemical type attribute of the current object
              */
            void SetChemicalType(std::string chemical_type);
            /*! \fn
              * A mutator function in order to set the description of the current object
              * Set the description_ attribute of the current atom
              * @param description The description attribute of the current object
              */
            void SetDescription(std::string description);
            /*! \fn
              * A mutator function in order to set the element symbol of the current object
              * Set the element_symbol_ attribute of the current atom
              * @param element_symbol The element symbol attribute of the current object
              */
            void SetElementSymbol(std::string element_symbol);
            /*! \fn
              * A mutator function in order to set the node of the current object
              * Set the node_ attribute of the current atom
              * @param node The node attribute of the current object
              */
            void SetNode(AtomNode* node);
            /*! \fn
              * A mutator function in order to set the id of the current object
              * Set the id_ attribute of the current atom
              * @param id The identification of the current object
              */
            void SetId(std::string id);
            /*! \fn
              * A mutator function in order to set the is_ring_ attribute of the current object
              * Set the is_ring_ attribute of the current atom
              * @param is_ring The boolean value representing if the current atom object is in aring or not
              */
            void SetIsRing(bool is_ring);

            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the atom contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cout);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            Residue* residue_;                      /*!< A pointer back to the residue that this atom belongs to >*/
            std::string name_;                      /*!< Name of the atom >*/
            CoordinateVector coordinates_;          /*!< Position of the atom >*/
            std::string chemical_type_;             /*!< A descriptor to determines the chemical type of the atom >*/
            std::string description_;               /*!< Short description of the atom >*/
            std::string element_symbol_;            /*!< Element symbol of the atom >*/
            AtomNode* node_;                        /*!< A Pointer to a node of the graph structure that indicates this atom >*/
            std::string id_;                        /*!< An identifier for an atom which is generated based on the type of the input file from which the structure has to be built
                                                      Mostly it is like "residue_name:atom_name" >*/
            bool is_ring_;                          /*!< A boolean value which represents if an atom is involved in a sugar ring or not. This attribute is set during the Sugar ID process >*/
            unsigned long long index_;              /*!< A unqiue index for each atom in an assembly >*/
    };
}

#endif // ATOM_HPP
