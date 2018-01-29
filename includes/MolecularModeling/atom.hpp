#ifndef ATOM_HPP
#define ATOM_HPP

#include <string>
#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>

#include "../GeometryTopology/coordinate.hpp"
#include "moleculardynamicatom.hpp"
#include "quantommechanicatom.hpp"
#include "dockingatom.hpp"
#include "atomnode.hpp"
#include "residue.hpp"

using namespace MolecularModeling;

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
      * List of GeometryTopology::Coordinates
      */
      typedef std::vector<GeometryTopology::Coordinate*> CoordinateVector;
      /*! \typedef
      * List of Atoms
      */
      typedef std::vector<Atom*> AtomVector;
      /*! \typedef
      * List of AtomTypes
      */
      typedef std::vector<std::string> AtomTypeVector;
      //////////////////////////////////////////////////////////
      //                       CONSTRUCTORS                   //
      //////////////////////////////////////////////////////////
      /*! \fn
      * Default constructor
      */
      Atom();
      /*! \fn
      * Constructors
      */
      Atom(Residue*, std::string, CoordinateVector);
      Atom(Residue*, std::string, GeometryTopology::Coordinate);
      /*! \fn
      * Copy constructor(*)
      */
      Atom(const Atom*);
      /*! \fn
      * Copy constructor(&)
      */
      Atom(const Atom&);
      //////////////////////////////////////////////////////////
      //                       DESTRUCTOR                     //
      //////////////////////////////////////////////////////////
      ~Atom();
      //////////////////////////////////////////////////////////
      //                       ACCESSORS                      //
      //////////////////////////////////////////////////////////
      /** \addtogroup Molecular_Data_Structure
      * @{
      */
      /*! \fn
      * An accessor function in order to access to the residue
      * @return residue_ attribute of the current object of this class
      */
      Residue* GetResidue() const;
      /*! \fn
      * An accessor function in order to access to the name
      * @return name_ attribute of the current object of this class
      */
      std::string GetName() const;
      /*! \fn
      * An accessor function in order to access to the coordinates
      * @return coordinates_ attribute of the current object of this class
      */
      CoordinateVector GetCoordinates() const;
      /*! \fn
      * An accessor function in order to access to the chemical_type
      * @return chemical_type_ attribute of the current object of this class
      */
      std::string GetChemicalType() const;
      /*! \fn
      * An accessor function in order to access to the description
      * @return description_ attribute of the current object of this class
      */
      std::string GetDescription() const;
      /*! \fn
      * An accessor function in order to access to the element symbol
      * @return element_symbol_ attribute of the current object of this class
      */
      std::string GetElementSymbol() const;
      /*! \fn
      * An accessor function in order to access to the node
      * @return node_ attribute of the current object of this class
      */
      AtomNode* GetNode() const;
      /*! \fn
      * An accessor function in order to access to the id
      * @return id_ attribute of the current object of this class
      */
      std::string GetId() const;
      /*! \fn
      * An accessor function in order to access to the is_ring_ attribute of the current object
      * @return is_ring_ attribute of the current object of this class
      */
      bool GetIsRing() const;
      /*! \fn
      * An accessor function in order to access to the index
      * @return index_ attribute of the current object of this class
      */
      unsigned long long GetIndex() const;
      /*! \fn                                     //Added by ayush on 12/11/17 for molecules in assembly
      * An accessor function in order to access to the atom types
      * @return atom_types_ attribute of the current object of this class
      */
      AtomTypeVector GetAtomTypes() const;
      /** @}*/
      //////////////////////////////////////////////////////////
      //                       MUTATOR                        //
      //////////////////////////////////////////////////////////
      /** \addtogroup Manipulators
      * @{
      */
      /*! \fn
      * A mutator function in order to set the residue of the current object
      * Set the residue_ attribute of the current atom
      * @param residue The residue attribute of the current object
      */
      void SetResidue(Residue*);
      /*! \fn
      * A mutator function in order to set the name of the current object
      * Set the name_ attribute of the current atom
      * @param name The name attribute of the current object
      */
      void SetName(std::string);
      /*! \fn
      * A mutator function in order to set the coordinates of the current object
      * Set the coordinates_ attribute of the current atom
      * @param coordinates The coordinates attribute of the current object
      */
      void SetCoordinates(CoordinateVector);
      /*! \fn
      * A function in order to add the coordinate to the current object
      * Set the coordinates_ attribute of the current atom
      * @param coordinate The coordinate of the current object
      */
      void AddCoordinate(GeometryTopology::Coordinate*);
      /*! \fn
      * A mutator function in order to set the chemical type of the current object
      * Set the chemical_type_ attribute of the current atom
      * @param chemical_type The chemical type attribute of the current object
      */
      void SetChemicalType(std::string);
      /*! \fn
      * A mutator function in order to set the description of the current object
      * Set the description_ attribute of the current atom
      * @param description The description attribute of the current object
      */
      void SetDescription(std::string);
      /*! \fn
      * A mutator function in order to set the element symbol of the current object
      * Set the element_symbol_ attribute of the current atom
      * @param element_symbol The element symbol attribute of the current object
      */
      void SetElementSymbol(std::string);
      /*! \fn
      * A mutator function in order to set the node of the current object
      * Set the node_ attribute of the current atom
      * @param node The node attribute of the current object
      */
      void SetNode(AtomNode*);
      /*! \fn
      * A mutator function in order to set the id of the current object
      * Set the id_ attribute of the current atom
      * @param id The identification of the current object
      */
      void SetId(std::string);
      /*! \fn
      * A mutator function in order to set the is_ring_ attribute of the current object
      * Set the is_ring_ attribute of the current atom
      * @param is_ring The boolean value representing if the current atom object is in aring or not
      */
      void SetIsRing(bool);
      //Added by ayush on 11/12/17 for molecules in assembly
      /*! \fn
      * A mutator function in order to set the atom types of the current object
      * Set the atom_types_ attribute of the current atom
      * @param atom_types The std::string stypes attribute of the current object
      */
      void SetAtomTypes(AtomTypeVector);
      //Added by ayush on 13/11/17 for molecules in assembly
      /*! \fn
      * A function to add a type to the AtomTypeVector of an atom
      * @param type The std::string type attribute of the atom
      */
      void AddAtomType(std::string);
      /*! \fn
      * A function to generate the index for an Atom.
      */
      unsigned long long generateAtomIndex();
      /** @}*/
      //////////////////////////////////////////////////////////
      //                       FUNCTIONS                      //
      //////////////////////////////////////////////////////////
      /*! \fn
      * A function to find connected Atoms to this Atom.
      * @param visitedAtoms The AtomVector used to find connected atoms.
      */
      void FindConnectedAtoms(AtomVector&);
      /*! \fn
      * A function to aget the distance between this Atom and another.
      * @param otherAtoms The Atom used to get the distance from this Atom.
      */
      double GetDistanceToAtom(Atom*);
      /*! \fn
      * A function to get the distance between this Atom's Coordinates
      * and another set of Coordinates.
      * @param coordinate The GeometryTopology::Coordinate used to get the distance from this Atom's Coordinates.
      */
      double GetDistanceToCoordinate(GeometryTopology::Coordinate*);
      //////////////////////////////////////////////////////////
      //                       DISPLAY FUNCTION               //
      //////////////////////////////////////////////////////////
      /*! \fn
      * A function to print out the atom contents in a structural format
      * Print out the information in a defined structure
      * @param out An output stream, the print result will be written in the given output stream
      */
      void Print(std::ostream& = std::cout); // @TODO DT - See the TODO below for operator<<
			//////////////////////////////////////////////////////////
      //                   OVERLOADED OPERATORS               //
      //////////////////////////////////////////////////////////
			// @TODO DT - Not sure if these are needed, but I am going to leave them
			// 				here for potential furthur discussion.
			//			The potential benefit would be to allow someone to make deep copies of an
			//				Atom object to an already initialized Atom object.
			// void operator=(const Atom&);
			// void operator=(const Atom*);

    private:
      //////////////////////////////////////////////////////////
      //                       ATTRIBUTES                     //
      //////////////////////////////////////////////////////////
      Residue* residue_;                 /*!< A pointer back to the residue that this atom belongs to >*/
      std::string name_;                 /*!< Name of the atom >*/
      CoordinateVector coordinates_;     /*!< Position of the atom >*/
      std::string chemical_type_;        /*!< A descriptor to determines the chemical type of the atom >*/
      std::string description_;          /*!< Short description of the atom >*/
      std::string element_symbol_;       /*!< Element symbol of the atom >*/
      AtomNode* node_;                   /*!< A Pointer to a node of the graph structure that indicates this atom >*/
      std::string id_;                   /*!< An identifier for an atom which is generated based on the type of the input file from which the structure has to be built
                                        			Mostly it is like "residue_name:atom_name" >*/
      bool is_ring_;                     /*!< A boolean value which represents if an atom is involved in a sugar ring or not. This attribute is set during the Sugar ID process >*/
      unsigned long long index_;         /*!< A unqiue index for each atom in an assembly >*/
      AtomTypeVector atom_types_;        /*!< List the atom type in an assembly >*/      //Added by ayush on 13/11/17 for molecules in assembly

			//////////////////////////////////////////////////////////
      //                   HELPER FUNCTIONS                   //
      //////////////////////////////////////////////////////////
			void Copy(const Atom*);
			void SetAttributes(Residue*, std::string, CoordinateVector,
													std::string, std::string, std::string, AtomNode*,
													std::string, bool, AtomTypeVector);
	};
	// @TODO DT - Get this working. For some reason it causes a Seg Fault.
	//			Ideally this function would be where the printing code is and the
	//			class Print function would call this operator.
	//			This would allow someone to do: std::cout << Atom << std::endl;
  // std::ostream& operator<<( std::ostream& out, const Atom& atom );
}

#endif // ATOM_HPP
