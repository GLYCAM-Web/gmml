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
            /*! \typedef
              * Vector of atom names in this cases excluded from interaction with another atom
              */
            typedef std::vector<std::string> ExcludedAtomNames;

            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            TopologyAtom();
            /*! \fn
              * Constructor with required parameters
              * @param atom_index Index of appearing of the atom in a topology file
              * @param atom_name Name of the atom that appears in the topology file
              * @param type Amber type of the atom that appears in the topology file
              * @param atom_charge Charge of the atom that appears in the topology file
              * @param atomic_number Atomic number of the atom that appears in the topology file
              * @param atom_mass Atom mass that appears in the topology file
              * @param excluded_atoms List of atoms that are excluded from interaction with the current atom
              * @param number_of_excluded_atoms The number of excluded atoms
              * @param radii Radius of the atom appears in the topology file
              * @param screen Screen value of the atom appears in the topology file
              * @param tree_chain_classification Tree chain classification of the atom in order to extract bonding information of the atoms in a residue
              * @param residue_name Name of the residue that the atom belongs to
              */
            TopologyAtom(int atom_index, std::string atom_name, std::string type, double atom_charge, int atomic_number, double atom_mass, ExcludedAtomNames excluded_atoms,
                         int number_of_excluded_atoms, double radii, double screen, std::string tree_chain_classification, std::string residue_name);

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
            /*! \fn
              * An accessor function in order to access to the residue name
              * @return residue_name_ attribute of the current object of this class
              */
            std::string GetResidueName();

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
              * A function in order to add the excluded atom to the current object
              * Set the excluded_atoms_ attribute of the current topology atom
              * @param excluded_atom The excluded atom to be added to the current object
              */
            void AddExcludedAtom(std::string excluded_atom);
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
            /*! \fn
              * A mutator function in order to set the residue name of the current object
              * Set the residue_name_ attribute of the current topology atom
              * @param residue_name The residue name attribute of the current object
              */
            void SetResidueName(std::string residue_name);

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
            void Print(std::ostream& out = std::cerr);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            std::string atom_name_;                     /*!< Name of the atom appears in a topology file >*/
            double atom_charge_;                        /*!< Charge of the atom that appears in a topology file >*/
            int atomic_number_;                         /*!< Atomic number of the atom that appears in a topology file >*/
            double atom_mass_;                          /*!< Mass of the atom  that appears in a topology file >*/
            std::vector<std::string> excluded_atoms_;   /*!< List of atom names that are excluded from interaction with the current atom >*/
            double radii_;                              /*!< Radius of the atom that appears in a topology file >*/
            double screen_;                             /*!< Screen value of the atom that appears in a topology file >*/
            std::string tree_chain_classification_;     /*!< Tree chain classification of the atom appears in a topology file in order to extract the bonding information inside a residue >*/
            int number_of_excluded_atoms_;              /*!< Number of atoms that are excluded from interaction with the current atom >*/
            std::string type_;                          /*!< Amber type of the atom that appears in a topology file >*/
            int index_;                                 /*!< Index of appearing o the atom in a topology file >*/
            std::string residue_name_;                  /*!< Name of the residue that the current atom belongs to >*/
    };
}

#endif // TOPOLOGYATOM_HPP
