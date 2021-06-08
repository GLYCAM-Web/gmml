#ifndef MOLECULARDYNAMICATOM_HPP
#define MOLECULARDYNAMICATOM_HPP

#include <string>
#include <iostream>

namespace MolecularModeling
{
    class MolecularDynamicAtom
    {
        public:
            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            MolecularDynamicAtom();
            MolecularDynamicAtom(MolecularDynamicAtom& moleculardynamicatom);

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
/** \addtogroup Molecular_Data_Structure
               * @{
               */
            /*! \fn
              * An accessor function in order to access to the atom type
              * @return atom_type_ attribute of the current object of this class
              */
            std::string GetAtomType();
            /*! \fn
              * An accessor function in order to access to the charge
              * @return charge_ attribute of the current object of this class
              */
            double GetCharge();
            /*! \fn
              * An accessor function in order to access to the mass
              * @return mass_ attribute of the current object of this class
              */
            double GetMass();
            /*! \fn
              * An accessor function in order to access to the radius
              * @return radius_ attribute of the current object of this class
              */
            double GetRadius();
/** @}*/
            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
/** \addtogroup Molecular_Data_Structure
               * @{
               */
            /*! \fn
              * A mutator function in order to set the atom type of the current object
              * Set the atom_type_ attribute of the current md atom
              * @param atom_type The atom type attribute of the current object
              */
            void SetAtomType(std::string atom_type);
            /*! \fn
              * A mutator function in order to set the charge of the current object
              * Set the charge_ attribute of the current md atom
              * @param charge The charge attribute of the current object
              */
            void SetCharge(double charge);
            /*! \fn
              * A mutator function in order to set the mass of the current object
              * Set the mass_ attribute of the current md atom
              * @param mass The mass attribute of the current object
              */
            void SetMass(double mass);
            /*! \fn
              * A mutator function in order to set the radius of the current object
              * Set the radius_ attribute of the current md atom
              * @param mass The radius attribute of the current object
              */
            void SetRadius(double radius);
/** @}*/
            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the md_atom contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cerr);

            void PrintHet(std::ostream& out = std::cerr);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            std::string atom_type_;             /*!< Atom type >*/
            double charge_;                     /*!< Charge of the atom >*/
            double mass_;                       /*!< Mass of the atom >*/
            double radius_;                     /*!< Radius of the atom >*/

    };
}

#endif // MOLECULARDYNAMICATOM_HPP
