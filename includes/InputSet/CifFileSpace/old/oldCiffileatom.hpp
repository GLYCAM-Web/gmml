#ifndef CIFFILEATOM_HPP
#define CIFFILEATOM_HPP

#include <string>
#include <vector>

#include "../../GeometryTopology/coordinate.hpp"

namespace CifFileSpace{

    class CifFileAtom
    {
        public:

            //////////////////////////////////////////////////////////
            //                    TYPE DEFINITION                   //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                       Constructor                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            CifFileAtom();

            /*! \fn
              * Constructor with given line about an atom name
              * @param name Atom name
              */
//            CifFileAtom(std::string& name);

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
            /*! \fn
              * An accessor function in order to access to the atom name of the current object
              * @return atom_name_ attribute of the current object of this class
              */
            std::string GetAtomName();
            /*! \fn
              * An accessor function in order to access to the alternate atom name of the current object
              * @return atom_name_ attribute of the current object of this class
              */
            std::string GetAlternateAtomName();
            /*! \fn
              * An accessor function in order to access to the element symbol of the current object
              * @return element_symbol_ attribute of the current object of this class
              */
            std::string GetElementSymbol();
            /*! \fn
              * An accessor function in order to access to the charge of the current object
              * @return charge_ attribute of the current object of this class
              */
            double GetCharge();
            /*! \fn
              * An accessor function in order to access to the coordinate of the current object
              * @return coordinate_ attribute of the current object of this class
              */
            GeometryTopology::Coordinate* GetCoordinate();
            /*! \fn
              * An accessor function in order to access to the coordinate ideal of the current object
              * @return coordinate_ideal_ attribute of the current object of this class
              */
            GeometryTopology::Coordinate* GetCoordinateIdeal();

            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function in order to set the atom name of the current object
              * @param atom_name The atom name
              */
            void SetAtomName(std::string atom_name);
            /*! \fn
              * A function in order to set the alternate atom name of the current object
              * @param alternate_atom_name The alternate atom name
              */
            void SetAlternateAtomName(std::string alternate_atom_name);
            /*! \fn
              * A function in order to set the element symbol of the current object
              * @param element_symbol The element symbol of the atom object
              */
            void SetElementSymbol(std::string element_symbol);
            /*! \fn
              * A function in order to set the charge of the current object
              * @param charge The charge of the atom
              */
            void SetCharge(double charge);
            /*! \fn
              * A mutator function in order to set the coordinate of the current object
              * Set the cooridnate_ attribute of the current cif file
              * @param cooridnate The coordinate of the current object
              */
            void SetCoordinate(GeometryTopology::Coordinate* coordinate);
            /*! \fn
              * A mutator function in order to set the ideal coordinate of the current object
              * Set the cooridnate_ideal_ attribute of the current cif file
              * @param cooridnate_ideal The ideal coordinate of the current object
              */
            void SetCoordinateIdeal(GeometryTopology::Coordinate* coordinate_ideal);

            //////////////////////////////////////////////////////////
            //                        FUNCTIONS                     //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to parse the contents of a given stream of a file
              * Parse the given stream and set the attributes of the current object accordingly
              * @param in_file A stream contains whole contents of a cif file
              */
//            bool Read(std::ifstream& in_file);

            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the cif file atom contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cout);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            std::string atom_name_;                           /*!< The The name of the atom object */
            std::string alternate_atom_name_;                 /*!< The alternate name of the current object */
            std::string element_symbol_;                      /*!< The element symbol of the current object */
            double charge_;                                   /*!< The charge of the current object */
            GeometryTopology::Coordinate* coordinate_;        /*!< The geometry coordinate of the current object */
            GeometryTopology::Coordinate* coordinate_ideal_;  /*!< The ideal geometry coordinate of the current object */
    };
}

#endif // CIFFILEATOM_HPP
