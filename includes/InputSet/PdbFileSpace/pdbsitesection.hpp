// Created by: Delaram Rahbarinia
// Modified by: Alireza Khatamian, Delaram Rahbarinia, Dave Montgomery

#ifndef PDBSITESECTION_HPP
#define PDBSITESECTION_HPP

#include <string>
#include <map>
#include <sstream>
#include <iostream>

namespace PdbFileSpace
{
    class PdbSiteCard;

    class PdbSiteSection
    {
        public:
            //////////////////////////////////////////////////////////
            //                       TYPE DEFINITION                //
            //////////////////////////////////////////////////////////
            /*! \typedef
              * A mapping between between site identifier and the site itself
              */
            typedef std::map<std::string, PdbSiteCard*> PdbSiteCardMap;

            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            PdbSiteSection();
            /*! \fn
              * Constructor with required parameters
              * @param stream_block The whole block of site card as a stream
              */
            PdbSiteSection(std::stringstream& stream_block);

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
/** \addtogroup Molecular_Data_Structure
               * @{
               */
            /*! \fn
              * An accessor function in order to access to the record name in a site card
              * @return record_name_ attribute of the current object of this class
              */
            std::string GetRecordName();
            /*! \fn
              * An accessor function in order to access to the residue site_cards in a site card
              * @return residue_site_cards_ attribute of the current object of this class
              */
            PdbSiteCardMap GetResidueSiteCards();
/** @}*/
            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
/** \addtogroup Manipulators
               * @{
               */
            /*! \fn
              * A mutator function in order to set the record name of the current object
              * Set the record_name_ attribute of the current site card
              * @param record_name The record name of the current object
              */
            void SetRecordName(const std::string record_name);
/** @}*/
            //////////////////////////////////////////////////////////
            //                        FUNCTIONS                     //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the pdb site card contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cerr);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            std::string record_name_;       /*!< Name of the record for site card that appears in first column of each line of a pdb file >*/
            PdbSiteCardMap residue_site_cards_;      /*!< Mapping of all site ids in a pdb file and their corresponding site residue object >*/

    };
}


#endif // PDBSITESECTION_HPP
