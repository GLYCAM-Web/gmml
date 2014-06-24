#ifndef ATOM_HPP
#define ATOM_HPP

#include <string>
#include <iostream>

#include "Geometry/coordinate.hpp"
#include "../../includes/MolecularModeling/md_atom.hpp"
#include "../../includes/MolecularModeling/qm_atom.hpp"
#include "../../includes/MolecularModeling/docking_atom.hpp"

namespace MolecularModeling
{
    class Residue;
    class Atom : public MD_Atom, public QM_Atom, public Docking_Atom
    {
        public:
            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            Atom();

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
              * An accessor function in order to access to the coordinate
              * @return coordinate_ attribute of the current object of this class
              */
            Geometry::Coordinate GetCoordinate();
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
              * A mutator function in order to set the coordinate of the current object
              * Set the coordinate_ attribute of the current atom
              * @param coordinate The coordinate attribute of the current object
              */
            void SetCoordinate(Geometry::Coordinate coordinate);
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
            Residue* residue_;
            std::string name_;
            Geometry::Coordinate coordinate_;
            std::string chemical_type_;
            std::string description_;
            std::string element_symbol_;

    };
}


#endif // ATOM_HPP
