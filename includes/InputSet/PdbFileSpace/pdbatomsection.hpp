// Created by: Delaram Rahbarinia
// Modified by: Alireza Khatamian, Delaram Rahbarinia, Dave Montgomery

#ifndef PDBATOMSECTION_HPP
#define PDBATOMSECTION_HPP

#include <string>
#include <sstream>
#include <map>
#include <vector>
#include <iostream>

namespace PdbFileSpace
{
    class PdbAtomCard;
    class PdbAtomSection
    {
        public:
            //////////////////////////////////////////////////////////
            //                       TYPE DEFINITION                //
            //////////////////////////////////////////////////////////
            /*! \typedef
              * A mapping between atom serial number and its its belonging factors and information in a model in a pdb file
              */
            typedef std::map<int, PdbAtomCard*> PdbAtomMap;
            /*! \typedef
              * A list of pdb atoms in the order of input pdb file
              */
            typedef std::vector<PdbAtomCard*> PdbAtomCardOrderVector;

            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            PdbAtomSection();
            /*! \fn
              * A constructor that get a stream block of atom card and parse the whole block to fill the related fields
              * @param stream_block A whole block of atoms belonging to a model in a pdb file
              * @param index
              */
            PdbAtomSection(std::stringstream& stream_block, std::string index);

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
            PdbAtomMap GetAtomCards();
            /*! \fn
              * An accessor function in order to access to the atoms in the atom card in the order of the original pdb file
              * @return ordered_atoms_ attribute of the current object of this class
              */
            PdbAtomCardOrderVector GetOrderedAtomCards();
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
            void SetAtomCards(PdbAtomMap atom_cards);
            /*! \fn
              * A mutator function in order to set the atoms of the current object in the order of the original pdb file
              * Set the ordered_atoms_ attribute of the current atom card
              * @param ordered_atoms attribute of the current object
              */
            void SetOrderedAtomCards(PdbAtomCardOrderVector ordered_atom_cards);
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
            void Print(std::ostream& out = std::cerr);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            std::string record_name_;           /*!< Record name of atom card: "ATOM" */
            PdbAtomMap atom_cards_;                  /*!< Map of all atom cards information that belong to a specific model in a pdb file by their serial numbers */
            PdbAtomCardOrderVector ordered_atom_cards_;  /*!< A list of atom cards in order of the input pdb file */

    };
}

#endif // PDBATOMSECTION_HPP
