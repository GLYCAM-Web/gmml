// Created by: Delaram Rahbarinia
// Modified by: Alireza Khatamian, Delaram Rahbarinia

#ifndef PDBATOMCARD_HPP
#define PDBATOMCARD_HPP

#include <string>
#include <sstream>
#include <map>
#include <vector>
#include <iostream>

namespace PdbFileSpace
{
    class PdbAtom;
    class PdbAtomCard
    {
        public:
            //////////////////////////////////////////////////////////
            //                       TYPE DEFINITION                //
            //////////////////////////////////////////////////////////
            /*! \typedef
              * A mapping between atom serial number and its its belonging factors and information in a model in a pdb file
              */
            typedef std::map<int, PdbAtom*> PdbAtomMap;
            /*! \typedef
              * A list of pdb atoms in the order of input pdb file
              */
            typedef std::vector<PdbAtom*> PdbAtomOrderVector;

            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            PdbAtomCard();
            /*! \fn
              * A constructor that get a stream block of atom card and parse the whole block to fill the related fields
              * @param stream_block A whole block of atoms belonging to a model in a pdb file
              * @param index
              */
            PdbAtomCard(std::stringstream& stream_block, std::string index);

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
/** \addtogroup Molecular_Data_Structure
              * @{
              */
            /*! \fn
              * An accessor function in order to access to the record name in a atom card
              * @return record_name_ attribute of the current object of this class
              */
            std::string GetRecordName();
            /*! \fn
              * An accessor function in order to access to the atoms in a atom card
              * @return atoms_ attribute of the current object of this class
              */
            PdbAtomMap GetAtoms();
            /*! \fn
              * An accessor function in order to access to the atoms in the atom card in the order of the original pdb file
              * @return ordered_atoms_ attribute of the current object of this class
              */
            PdbAtomOrderVector GetOrderedAtoms();
/** @}*/
            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
/** \addtogroup Manipulators
              * @{
              */
            /*! \fn
              * A mutator function in order to set the record name of the current object
              * Set the record_name_ attribute of the current atom card
              * @param record_name The record name attribute of the current object
              */
            void SetRecordName(const std::string record_name);
            /*! \fn
              * A mutator function in order to set the atoms of the current object
              * Set the atoms_ attribute of the current atom card
              * @param atoms Atom attribute of the current object
              */
            void SetAtoms(PdbAtomMap atoms);
            /*! \fn
              * A mutator function in order to set the atoms of the current object in the order of the original pdb file
              * Set the ordered_atoms_ attribute of the current atom card
              * @param ordered_atoms attribute of the current object
              */
            void SetOrderedAtoms(PdbAtomOrderVector ordered_atoms);
/** @}*/
            //////////////////////////////////////////////////////////
            //                        FUNCTIONS                     //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the atom card contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cout);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            std::string record_name_;           /*!< Record name of atom card: "ATOM" */
            PdbAtomMap atoms_;                  /*!< Map of all atoms informatin that belong to a specific model in a pdb file by their serial numbers */
            PdbAtomOrderVector ordered_atoms_;  /*!< A list of atoms in order of the input pdb file */

    };
}

#endif // PDBATOMCARD_HPP
