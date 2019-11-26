#ifndef TOPOLOGYRESIDUE_HPP
#define TOPOLOGYRESIDUE_HPP

#include <string>
#include <iostream>
#include <vector>
#include <map>

namespace TopologyFileSpace
{
    class TopologyAtom;

    class TopologyResidue
    {
        public:
            //////////////////////////////////////////////////////////
            //                    TYPE DEFINITION                   //
            //////////////////////////////////////////////////////////
            typedef std::vector<TopologyAtom*> TopologyAtomVector;

            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            TopologyResidue();
            /*! \fn
              * Constructor with required parameters
              * @param residue_name Residue name in a topology file
              * @param atoms Atoms belonging to the current residue
              * @param starting_atom_index Index of the first atom that belongs to the residue among all atoms in a topology file
              */
            TopologyResidue(std::string residue_name, TopologyAtomVector atoms, int index, int starting_atom_index);

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
            /*! \fn
              * An accessor function in order to access to the residue name
              * @return residue_name_ attribute of the current object of this class
              */
            std::string GetResidueName();
            /*! \fn
              * An accessor function in order to access to the atoms
              * @return atoms_ attribute of the current object of this class
              */
            TopologyAtomVector GetAtoms();
            /*! \fn
              * An accessor function in order to access to an atom of the current object using index
              * @param index Index of the desired atom
              * @return atom Atom of the current object with the given index
              */
            TopologyAtom* GetAtomByIndex(int index);
            /*! \fn
              * An accessor function in order to access to the index
              * @return index_ attribute of the current object of this class
              */
            int GetIndex();
            /*! \fn
              * An accessor function in order to access to the starting atom index
              * @return starting_atom_index_ attribute of the current object of this class
              */
            int GetStartingAtomIndex();
            /*! \fn
              * An accessor function in order to access whether a residue is a Solvent, returns true if a residue is a solvent.
              * @return is_residue_solvent_ attribute of the current object of this class
              */
            bool GetIsResidueSolvent();

            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A mutator function in order to set the residue name of the current object
              * Set the residue_name_ attribute of the current topology residue
              * @param residue_name The residue name attribute of the current object
              */
            void SetResidueName(std::string residue_name);
            /*! \fn
              * A function in order to add an atom to the current object
              * Add a new entry to the atoms_ attribute of the current topology residue
              * @param atom A new atom to be added to the current topology residue
              */
            void AddAtom(TopologyAtom* atom);
            /*! \fn
              * A mutator function in order to set the index of the current object
              * Set the index_ attribute of the current topology residue
              * @param index The index attribute of the current object
              */
            void SetIndex(int index);
            /*! \fn
              * A mutator function in order to set the starting atom residue of the current object
              * Set the starting_atom_residue_ attribute of the current topology residue
              * @param starting_atom_index The starting atom index attribute of the current object
              */
            void SetStartingAtomIndex(int starting_atom_index);
            /*! \fn
              * A mutator function in order to set true if a topology residue is a solvent of the current object
              * Set the is_residue_solvent_ attribute of the current  residue
              * @param is_residue_solvent The is_residue_solvent_ attribute of the current object
              */
            void SetIsResidueSolvent(bool is_residue_solvent);

            //////////////////////////////////////////////////////////
            //                        FUNCTIONS                     //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the topology residue contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cerr);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            std::string residue_name_;              /*!< Residue name >*/
            TopologyAtomVector atoms_;                 /*!< List of atoms those belong to this residue >*/
            int index_;                             /*!< Index of this residue among the others in a topology file >*/
            int starting_atom_index_;               /*!< Index of the first atom belongs to this residue among all other atoms in a topology file >*/
            bool is_residue_solvent_;		 /*!< Is Residue Solvent>*/

    };
}

#endif // TOPOLOGYRESIDUE_HPP
