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
        typedef std::vector<PdbAtomSection*> PdbAtomSectionVector;
        /*! \typedef
         * List of heterogen atom cards
         */
        typedef std::vector<PdbHeterogenAtomSection*> HeterogenAtomSectionVector;

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
        PdbAtomSectionVector GetAtomCards();
        /*! \fn
         * An accessor function in order to access to the heterogen atoms in a model residue set
         * @return heterogen_atoms_ attribute of the current object of this class
         */
        HeterogenAtomSectionVector GetHeterogenAtomCards();
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
        void SetAtomCards(PdbAtomSectionVector atoms);
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
        void SetHeterogenAtoms(HeterogenAtomSectionVector heterogen_atoms);
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
        PdbAtomSectionVector atomSectionVector_; /*!< I have no idea what this is, but it isn't what it said it was. >*/
        HeterogenAtomSectionVector heterogenAtomsSectionVeector_; /*!< It is what it is. >*/
    };
} // namespace PdbFileSpace

#endif // PDBMODELRESIDUESET_HPP
