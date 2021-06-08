// Created by: Delaram Rahbarinia
// Modified by: Alireza Khatamian, Delaram Rahbarinia, Dave Montgomery

#ifndef PDBLINKSECTION_HPP
#define PDBLINKSECTION_HPP

#include <string>
#include <vector>
#include <sstream>
#include <iostream>

namespace PdbFileSpace
{
    class PdbLinkCard;

    class PdbLinkSection
    {
        public:
            //////////////////////////////////////////////////////////
            //                       TYPE DEFINITION                //
            //////////////////////////////////////////////////////////
            /*! \typedef
              * List of links
              */
            typedef std::vector< PdbLinkCard* > LinkCardVector;

            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            PdbLinkSection();
            /*! \fn
              * Constructor with required parameters
              * @param stream_block
              */
            PdbLinkSection(std::stringstream& stream_block);

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
/** \addtogroup Molecular_Data_Structure
              * @{
              */
            /*! \fn
              * An accessor function in order to access to the record name in a link card
              * @return record_name_ attribute of the current object of this class
              */
            std::string GetRecordName();
            /*! \fn
              * An accessor function in order to access to the resdiue links in a link card
              * @return resdiue_links_ attribute of the current object of this class
              */
            LinkCardVector GetResidueLinkCards();
/** @}*/
            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
/** \addtogroup Manipulators
              * @{
              */
            /*! \fn
              * A mutator function in order to set the record name of the current object
              * Set the record_name_ attribute of the current link card
              * @param record_name The record name of the current object
              */
            void SetRecordName(std::string record_name);
            /*! \fn
              * A mutator function in order to set the residue links of the current object
              * Set the residue_link_cards_ attribute of the current link card
              * @param residue_links The residue links of the current object
              */
            void SetResidueLinkCards(LinkCardVector residue_link_cards);
            /*! \fn
              * A function in order to add the residue link to the current object
              * Set the residue_link_ attribute of the current link card
              * @param residue_link The residue link of the current object
              */
            void AddResidueLinkCard(PdbLinkCard* residue_link_card);
/** @}*/
            //////////////////////////////////////////////////////////
            //                        FUNCTIONS                     //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the link card contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cerr);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            std::string record_name_;       /*!< Record name of a link card in a pdb file >*/
            LinkCardVector residue_link_cards_;      /*!< List of links involving in a pdb file >*/

    };
}

#endif // PDBLINKSECTION_HPP
