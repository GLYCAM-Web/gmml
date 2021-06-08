// Created by: Delaram Rahbarinia
// Modified by: Alireza Khatamian, Delaram Rahbarinia

#ifndef PDBMODELRESIDUESET_HPP
#define PDBMODELRESIDUESET_HPP

#include <string>
#include <vector>
#include <sstream>
#include <iostream>

namespace PdbFileSpace
{
    class PdbAtomSection;
    class PdbHeterogenAtomSection;

    class PdbModelResidueSet
    {
        public:
            //////////////////////////////////////////////////////////
            //                       TYPE DEFINITION                //
            //////////////////////////////////////////////////////////
            /*! \typedef
              * List of atom cards in a model residue set
              */
            typedef std::vector<PdbAtomSection*> AtomCardVector;
            /*! \typedef
              * List of heterogen atom cards
              */
            typedef std::vector<PdbHeterogenAtomSection*> HeterogenAtomCardVector;

            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            PdbModelResidueSet();
            /*! \fn
              * Constructor with required parameters
              * @param residue_stream_block
              */
            PdbModelResidueSet(std::stringstream& residue_set_block);

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
/** \addtogroup Molecular_Data_Structure
               * @{
               */
            /*! \fn
              * An accessor function in order to access to the atoms in a model residue set
              * @return atoms_ attribute of the current object of this class
              */
            AtomCardVector GetAtomCards();
            /*! \fn
              * An accessor function in order to access to the heterogen atoms in a model residue set
              * @return heterogen_atoms_ attribute of the current object of this class
              */
            HeterogenAtomCardVector GetHeterogenAtomCards();
/** @}*/
            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
/** \addtogroup Manipulators
               * @{
               */
            /*! \fn
              * A mutator function in order to set the atoms of the current object
              * Set the atoms_ attribute of the current model residue set
              * @param atoms The atoms attribute of the current object
              */
            void SetAtomCards(AtomCardVector atoms);
            /*! \fn
              * A function in order to add the atom to the current object
              * Set the atom_ attribute of the current model residue set
              * @param atom The atom attribute of the current object
              */
            void AddAtom(PdbAtomSection* atom);
            /*! \fn
              * A mutator function in order to set the heterogen atoms of the current object
              * Set the heterogen_atoms_ attribute of the current model residue set
              * @param heterogen_atoms The heterogen atoms of the current object
              */
            void SetHeterogenAtoms(HeterogenAtomCardVector heterogen_atoms);
            /*! \fn
              * A function in order to add the heterogen atom to the current object
              * Set the heterogen_atom_ attribute of the current model residue set
              * @param heterogen_atom The heterogen atom attribute of the current object
              */
            void AddHeterogenAtom(PdbHeterogenAtomSection* heterogen_atom);
/** @}*/
            //////////////////////////////////////////////////////////
            //                        FUNCTIONS                     //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the model residue set contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cerr);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            AtomCardVector atoms_;                          /*!< List of atom cards that are in a model >*/
            HeterogenAtomCardVector heterogen_atoms_;       /*!< List of heterogen atom cards that are in a model >*/

    };
}

#endif // PDBMODELRESIDUESET_HPP
