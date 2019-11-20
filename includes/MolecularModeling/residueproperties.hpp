#ifndef RESIDUEPROPERTIES_HPP
#define RESIDUEPROPERTIES_HPP

#include <string>
#include <iostream>


namespace MolecularModeling
{
    class ResidueProperties
    {
        public:
            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            ResidueProperties();

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
            /*! \fn
              * An accessor function in order to access whether a residue is a Solvent, returns true if a residue is a solvent.
              * @return is_residue_solvent_ attribute of the current object of this class
              */
            bool GetIsResidueSolvent();


            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A mutator function in order to set true if a residue is a solvent of the current object
              * Set the is_residue_solvent_ attribute of the current ResiduePorperties residue
              * @param is_residue_solvent The is_residue_solvent_ attribute of the current object
              */
            void SetIsResidueSolvent(bool is_residue_solvent);

            /////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the md_atom contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cerr);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
           bool is_residue_solvent_;		 /*!< Is Residue Solvent>*/

    };
}

#endif // RESIDUEPROPERTIES_HPP

