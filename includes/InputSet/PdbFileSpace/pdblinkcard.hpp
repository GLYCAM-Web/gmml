// Created by: Delaram Rahbarinia
// Modified by: Alireza Khatamian, Delaram Rahbarinia

#ifndef PDBLINKCARD_HPP
#define PDBLINKCARD_HPP

#include <string>
#include <vector>
#include <sstream>
#include <iostream>

namespace PdbFileSpace
{
    class PdbLink;

    class PdbLinkCard
    {
        public:
            //////////////////////////////////////////////////////////
            //                       TYPE DEFINITION                //
            //////////////////////////////////////////////////////////
            /*! \typedef
              * List of links
              */
            typedef std::vector< PdbLink* > LinkVector;

            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            PdbLinkCard();
            /*! \fn
              * Constructor with required parameters
              * @param stream_block
              */
            PdbLinkCard(std::stringstream& stream_block);

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
            LinkVector GetResidueLinks();
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
              * Set the residue_links_ attribute of the current link card
              * @param residue_links The residue links of the current object
              */
            void SetResidueLinks(LinkVector residue_links);
            /*! \fn
              * A function in order to add the residue link to the current object
              * Set the residue_link_ attribute of the current link card
              * @param residue_link The residue link of the current object
              */
            void AddResidueLink(PdbLink* residue_link);
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
            void Print(std::ostream& out = std::cout);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            std::string record_name_;       /*!< Record name of a link card in a pdb file >*/
            LinkVector residue_links_;      /*!< List of links involving in a pdb file >*/

    };
}

#endif // PDBLINKCARD_HPP
