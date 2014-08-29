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
            typedef std::map<std::string, TopologyAtom*> TopologyAtomMap;

            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            TopologyResidue();
            /*! \fn
              * Constructor with required parameters
              * @param residue_name
              * @param atoms
              * @param starting_atom_index
              */
            TopologyResidue(std::string residue_name, TopologyAtomMap atoms, int index, int starting_atom_index);

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
            TopologyAtomMap GetAtoms();
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
            void Print(std::ostream& out = std::cout);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            std::string residue_name_;
            TopologyAtomMap atoms_;
            int index_;
            int starting_atom_index_;

    };
}

#endif // TOPOLOGYRESIDUE_HPP
