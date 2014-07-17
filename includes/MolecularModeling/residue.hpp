#ifndef RESIDUE_HPP
#define RESIDUE_HPP

#include <string>
#include <iostream>
#include <vector>


namespace MolecularModeling
{
    class Assembly;
    class Atom;
    class Residue
    {
        public:
            //////////////////////////////////////////////////////////
            //                    TYPE DEFINITION                   //
            //////////////////////////////////////////////////////////
            typedef std::vector<Atom*> AtomVector;

            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            Residue();

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
            /*! \fn
              * An accessor function in order to access to the assembly
              * @return assembly_ attribute of the current object of this class
              */
            Assembly* GetAssembly();
            /*! \fn
              * An accessor function in order to access to the name
              * @return name_ attribute of the current object of this class
              */
            std::string GetName();
            /*! \fn
              * An accessor function in order to access to the atoms
              * @return atoms_ attribute of the current object of this class
              */
            AtomVector GetAtoms();
            /*! \fn
              * An accessor function in order to access to the head atoms
              * @return head_atoms_ attribute of the current object of this class
              */
            AtomVector GetHeadAtoms();
            /*! \fn
              * An accessor function in order to access to the tail atoms
              * @return tail_atoms_ attribute of the current object of this class
              */
            AtomVector GetTailAtoms();
            /*! \fn
              * An accessor function in order to access to the chemical type
              * @return chemical_residue_ attribute of the current object of this class
              */
            std::string GetChemicalType();
            /*! \fn
              * An accessor function in order to access to the description
              * @return description_ attribute of the current object of this class
              */
            std::string GetDescription();            

            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A mutator function in order to set the assembly of the current object
              * Set the assembly_ attribute of the current residue
              * @param assembly The assembly attribute of the current object
              */
            void SetAssembly(Assembly* assembly);
            /*! \fn
              * A mutator function in order to set the name of the current object
              * Set the name_ attribute of the current residue
              * @param name The name attribute of the current object
              */
            void SetName(std::string name);
            /*! \fn
              * A mutator function in order to set the atoms of the current object
              * Set the atoms_ attribute of the current residue
              * @param atoms The atoms attribute of the current object
              */
            void SetAtoms(AtomVector atoms);
            /*! \fn
              * A function in order to add the atom to the current object
              * Set the atoms_ attribute of the current residue
              * @param atom The atom of the current object
              */
            void AddAtom(Atom* atom);
            /*! \fn
              * A mutator function in order to set the head atoms of the current object
              * Set the head_atoms_ attribute of the current residue
              * @param head_atoms The head atoms attribute of the current object
              */
            void SetHeadAtoms(AtomVector head_atoms);
            /*! \fn
              * A function in order to add the head atom to the current object
              * Set the head_atom_ attribute of the current residue
              * @param head_atom The head atom of the current object
              */
            void AddHeadAtom(Atom* head_atom);
            /*! \fn
              * A mutator function in order to set the tail atoms of the current object
              * Set the tail_atoms_ attribute of the current residue
              * @param tail_atoms The head atoms attribute of the current object
              */
            void SetTailAtoms(AtomVector tail_atoms);
            /*! \fn
              * A function in order to add the tail atom to the current object
              * Set the tail_atoms_ attribute of the current residue
              * @param tail_atom The tail atom of the current object
              */
            void AddTailAtom(Atom* tail_atom);
            /*! \fn
              * A mutator function in order to set the chemical type of the current object
              * Set the chemical_type_ attribute of the current residue
              * @param chemical_type The chemical type attribute of the current object
              */
            void SetChemicalType(std::string chemical_type);
            /*! \fn
              * A mutator function in order to set the description of the current object
              * Set the description_ attribute of the current residue
              * @param description The description attribute of the current object
              */
            void SetDescription(std::string description);

            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the md_atom contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cout);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            Assembly* assembly_;
            std::string name_;
            AtomVector atoms_;
            AtomVector head_atoms_;
            AtomVector tail_atoms_;
            std::string chemical_type_;
            std::string description_;

    };
 }


#endif // RESIDUE_HPP
