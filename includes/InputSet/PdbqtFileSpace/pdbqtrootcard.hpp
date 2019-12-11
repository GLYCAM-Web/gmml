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

            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            PdbqtRootCard();
            /*! \fn
              * A constructor that get a root block of model card and parse the whole block to fill the related fields
              * @param root_block A whole block of atoms belonging to a root card in a pdbqt file
              */
            PdbqtRootCard(std::ifstream& root_block, std::vector<PdbqtFileSpace::PdbqtAtomCard*>& ACV);

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
            /*! \fn
              * An accessor function in order to access to the record name in a root card
              * @return record_name_ attribute of the current object of this class
              */
            std::string GetRecordName();

            /*! \fn
              * An accessor function in order to access to the pdb qt root atoms
              * @return root_atoms_ attribute of the current object of this class
              */
            PdbqtAtomCard* GetRootAtoms();

            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A mutator function in order to set the record name of the current object
              * Set the record_name_ attribute of the current qt root card
              * @param record_name The record name attribute of the current object
              */
            void SetRecordName(const std::string record_name);
            /*! \fn
              * A mutator function in order to set root atoms of the current object
              * Set the root_cards_ attribute of the current pdbqt root card
              * @param root_atoms The root_atoms attribute of the current object
              */
            void SetRootAtoms(PdbqtAtomCard* root_atoms);

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
            void Print(std::ostream& out = std::cerr);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            std::string record_name_;           /*!< Record name of qt atom card: "ATOM" */
            PdbqtAtomCard* root_atoms_; /*!< List of root atoms in a pdbqt root card >*/
    };
}

#endif // PDBQTROOTCARD_HPP
