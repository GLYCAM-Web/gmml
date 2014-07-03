#ifndef TOPOLOGYATOM_HPP
#define TOPOLOGYATOM_HPP

#include <string>
#include <iostream>
#include <vector>
#include <map>

namespace TopologyFileSpace
{
    class TopologyAtom
    {
        public:
            //////////////////////////////////////////////////////////
            //                       TYPE DEFINITION                //
            //////////////////////////////////////////////////////////
            typedef std::vector<std::string> ExcludedAtomNames;

            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            TopologyAtom();
            TopologyAtom(int atom_index, std::string atom_name, std::string type, double atom_charge, int atomic_number, double atom_mass, ExcludedAtomNames excluded_atoms,
                         int number_of_excluded_atoms, double radii, double screen, std::string tree_chain_classification);

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
            /*! \fn
              * An accessor function in order to access to the atom name
              * @return atom_name_ attribute of the current object of this class
              */
            std::string GetAtomName();
            /*! \fn
              * An accessor function in order to access to the atom charge
              * @return atom_charge_ attribute of the current object of this class
              */
            double GetAtomCharge();
            /*! \fn
              * An accessor function in order to access to the atomic number
              * @return atomic_number_ attribute of the current object of this class
              */
            int GetAtomicNumber();
            /*! \fn
              * An accessor function in order to access to the atom mass
              * @return atom_mass_ attribute of the current object of this class
              */
            double GetAtomMass();
            /*! \fn
              * An accessor function in order to access to the excluded atoms
              * @return excluded_atoms_ attribute of the current object of this class
              */
            std::vector<std::string> GetExcludedAtoms();
            /*! \fn
              * An accessor function in order to access to the atom radii
              * @return radii_ attribute of the current object of this class
              */
            double GetRadii();
            /*! \fn
              * An accessor function in order to access to the atom screen
              * @return screen_ attribute of the current object of this class
              */
            double GetScreen();
            /*! \fn
              * An accessor function in order to access to the tree chain classification
              * @return tree_chain_classification_ attribute of the current object of this class
              */
            std::string GetTreeChainClassification();
            /*! \fn
              * An accessor function in order to access to the number of excluded atoms
              * @return number_of_excluded_atoms_ attribute of the current object of this class
              */
            int GetNumberOfExcludedAtoms();
            /*! \fn
              * An accessor function in order to access to the type
              * @return type_ attribute of the current object of this class
              */
            std::string GetType();
            /*! \fn
              * An accessor function in order to access to the index
              * @return index_ attribute of the current object of this class
              */
            int GetIndex();

            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A mutator function in order to set the atom name of the current object
              * Set the atom_name_ attribute of the current topology atom
              * @param atom_name The bond type attribute of the current object
              */
            void SetAtomName(std::string atom_name);
            /*! \fn
              * A mutator function in order to set the atom charge of the current object
              * Set the index_ attribute of the current topology atom
              * @param atom_charge The atom charge attribute of the current object
              */
            void SetAtomCharge(double atom_charge);
            /*! \fn
              * A mutator function in order to set the atomic number of the current object
              * Set the atomic_number_ attribute of the current topology atom
              * @param atomic_number The atomic_number attribute of the current object
              */
            void SetAtomicNumber(int atomic_number);
            /*! \fn
              * A mutator function in order to set the atom mass of the current object
              * Set the atom_mass_ attribute of the current topology atom
              * @param atom_mass The atom mass attribute of the current object
              */
            void SetAtomMass(double atom_mass);
            /*! \fn
              * A mutator function in order to set the excluded atoms of the current object
              * Set the excluded_atoms_ attribute of the current topology atom
              * @param excluded_atoms The excluded atoms attribute of the current object
              */
            void SetExcludedAtoms(std::vector<std::string> excluded_atoms);
            /*! \fn
              * A mutator function in order to set the radii of the current object
              * Set the radii_ attribute of the current topology atom
              * @param radii The radii attribute of the current object
              */
            void SetRadii(double radii);
            /*! \fn
              * A mutator function in order to set the screen of the current object
              * Set the screen_ attribute of the current topology atom
              * @param screen The screen attribute of the current object
              */
            void SetScreen(double screen);
            /*! \fn
              * A mutator function in order to set the tree chain classification of the current object
              * Set the tree_chain_classification_ attribute of the current topology atom
              * @param tree_chain_classification The tree chain classification attribute of the current object
              */
            void SetTreeChainClasification(std::string tree_chain_classification);
            /*! \fn
              * A mutator function in order to set the number of excluded atoms of the current object
              * Set the number_of_excluded_atoms_ attribute of the current topology atom
              * @param number_of_excluded_atoms The number of excluded atoms attribute of the current object
              */
            void SetNumberOfExcludedAtoms(int number_of_excluded_atoms);
            /*! \fn
              * A mutator function in order to set the type of the current object
              * Set the type_ attribute of the current topology atom
              * @param type The type attribute of the current object
              */
            void SetType(std::string type);
            /*! \fn
              * A mutator function in order to set the index of the current object
              * Set the index_ attribute of the current topology atom
              * @param index The index attribute of the current object
              */
            void SetIndex(int index);

            //////////////////////////////////////////////////////////
            //                        FUNCTIONS                     //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the atom contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cout);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            std::string atom_name_;
            double atom_charge_;
            int atomic_number_;
            double atom_mass_;
            std::vector<std::string> excluded_atoms_;
            double radii_;
            double screen_;
            std::string tree_chain_classification_;
            int number_of_excluded_atoms_;
            std::string type_;
            int index_;
    };
}

#endif // TOPOLOGYATOM_HPP
