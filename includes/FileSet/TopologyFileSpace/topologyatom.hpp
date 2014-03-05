#ifndef TOPOLOGYATOM_HPP
#define TOPOLOGYATOM_HPP

#include <string>
#include <iostream>
#include <vector>
#include <map>

namespace TopologyFileSpace
{
    class TopologyAtomType;

    class TopologyAtom
    {
        public:
            //////////////////////////////////////////////////////////
            //                       TYPE DEFINITION                //
            //////////////////////////////////////////////////////////
            typedef std::map<std::string, double> TopologyLennardJonesMap;

            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            TopologyAtom();

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
              * An accessor function in order to access to the atom type
              * @return atom_type_ attribute of the current object of this class
              */
            TopologyAtomType* GetAtomType();
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
            char GetTreeChainClassification();
            /*! \fn
              * An accessor function in order to access to the lennard jones A
              * @return lennard_jones_a_ attribute of the current object of this class
              */
            TopologyLennardJonesMap GetLennardJonesA();
            /*! \fn
              * An accessor function in order to access to the lennard jones b
              * @return lennard_jones_a_ attribute of the current object of this class
              */
            TopologyLennardJonesMap GetLennardJonesB();

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
              * A mutator function in order to set the atom type of the current object
              * Set the atom_type_ attribute of the current topology atom
              * @param atom_type The atom type attribute of the current object
              */
            void SetAtomType(TopologyAtomType atom_type);
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
            void SetTreeChainClasification(char tree_chain_classification);

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
            TopologyAtomType* atom_type_;
            std::vector<std::string> excluded_atoms_;
            double radii_;
            double screen_;
            char tree_chain_classification_;
            TopologyLennardJonesMap lennard_jones_a_;
            TopologyLennardJonesMap lennard_jones_b_;
    };
}

#endif // TOPOLOGYATOM_HPP
