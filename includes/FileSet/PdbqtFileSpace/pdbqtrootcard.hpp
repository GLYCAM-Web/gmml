#ifndef PDBQTROOTCARD_HPP
#define PDBQTROOTCARD_HPP

#include <string>
#include <iostream>
#include <vector>

namespace PdbqtFileSpace
{
    class PdbqtAtomCard;
    class PdbqtRootCard
    {
        public:
            //////////////////////////////////////////////////////////
            //                       TYPE DEFINITION                //
            //////////////////////////////////////////////////////////
            /*! \typedef
              * A vector of pdbqt root cards
              */
            typedef std::vector<PdbqtAtomCard*> AtomCardVector;

            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            PdbqtRootCard();

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
            /*! \fn
              * An accessor function in order to access to the pdb qt root atoms
              * @return root_atoms_ attribute of the current object of this class
              */
            AtomCardVector GetRootAtoms();

            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A mutator function in order to set root atoms of the current object
              * Set the root_cards_ attribute of the current pdbqt root card
              * @param root_atoms The root_atoms attribute of the current object
              */
            void SetRootAtoms(AtomCardVector root_atoms);
            /*! \fn
              * A function in order to add the root atom to the current object
              * Set the root_atoms_ attribute of the current pdbqt root card
              * @param root_atoms_ The root atom of the current object
              */
            void AddRootAtom(PdbqtAtomCard* root_atom);

            //////////////////////////////////////////////////////////
            //                        FUNCTIONS                     //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the root card contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cout);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            AtomCardVector root_atoms_; /*!< List of root atoms in a pdbqt root card >*/
    };
}

#endif // PDBQTROOTCARD_HPP
