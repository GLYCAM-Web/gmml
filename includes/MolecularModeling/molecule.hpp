#ifndef MOLECULE_HPP
#define MOLECULE_HPP

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include "atom.hpp"
#include "residue.hpp"


namespace MolecularModeling
{
    class Atom;
    class Residue;
    class Molecule; // Forward declare for the vector typedef:
    typedef std::vector<MolecularModeling::Molecule*> MoleculeVector;
    class Molecule
    {
        public:
            //////////////////////////////////////////////////////////
            //                    TYPE DEFINITION                   //
            //////////////////////////////////////////////////////////
          //  typedef std::vector<Atom*> AtomVector;
          //  typedef std::vector<Residue*> ResidueVector;

            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            Molecule();

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////

            /*! \fn
              * An accessor function in order to access index of a molecule in an assembly
              * @return molecule_index_ attribute of the current object of this class
              */
            int GetMoleculeIndex();

            /*! \fn
              * An accessor function in order to access to the atoms of the current molecule object
              * @return molecule_atoms_ attribute of the current object of this class
              */
           AtomVector GetMoleculeAtoms();

            /*! \fn
              * An accessor function in order to access to the residues of the current molecule object
              * @return molecule_residues_ attribute of the current object of this class
              */
            ResidueVector GetMoleculeResidues();


            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A mutator function in order to set the molecule index of the current molecule object
              * Set the molecule_index_ attribute of the current molecule
              * @param molecule_index The index attribute of the current molecule object
              */
            void SetMoleculeIndex(int molecule_index);

            /*! \fn
              * A mutator function in order to set the atoms of the current molecule object
              * Set the molecule_atoms_ attribute of the current molecule
              * @param molecule_atoms The molecule_atoms attribute of the current object
              */
            void SetMoleculeAtoms(AtomVector molecule_atoms);

            /*! \fn
              * A mutator function in order to set the residues of the current molecule object
              * Set the molecule_residues_ attribute of the current molecule
              * @param molecule_residues The molecule_residues attribute of the current object
              */
            void SetMoleculeResidues(ResidueVector molecule_residues);


            //////////////////////////////////////////////////////////
            //                       FUNCTIONS                      //
            //////////////////////////////////////////////////////////


            /////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the current molecule object
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cerr);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            int molecule_index_;                            /*!< An index of the molecule in an assembly >*/
            AtomVector molecule_atoms_;                     /*!< List of all atoms in the current molecule object >*/
            ResidueVector molecule_residues_;               /*!< List of residues involved in the current object of molecule >*/
    };
}

#endif // MOLECULE_HPP
