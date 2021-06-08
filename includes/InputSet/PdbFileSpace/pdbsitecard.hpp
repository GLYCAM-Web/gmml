// Created by: Delaram Rahbarinia
// Modified by: Alireza Khatamian, Delaram Rahbarinia, Dave Montgomery

#ifndef PDBSITECARD_HPP
#define PDBSITECARD_HPP

#include <string>
#include <sstream>
#include <vector>
#include <iostream>

namespace PdbFileSpace
{
    class PdbSiteResidue;

    class PdbSiteCard
    {
        public:
            //////////////////////////////////////////////////////////
            //                       TYPE DEFINITION                //
            //////////////////////////////////////////////////////////
            /*! \typedef
              * Vector of site residues
              */
            typedef std::vector< PdbSiteResidue* > SiteResidueVector;

            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            PdbSiteCard();
            /*! \fn
              * Constructor with required parameters
              * @param site_block
              */
            PdbSiteCard(std::stringstream& site_block);

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
/** \addtogroup Molecular_Data_Structure
               * @{
               */
            /*! \fn
              * An accessor function in order to access to the site id in a pdb site
              * @return site_id_ attribute of the current object of this class
              */
            std::string GetSiteId();
            /*! \fn
              * An accessor function in order to access to the residues in a pdb site
              * @return residues_ attribute of the current object of this class
              */
            SiteResidueVector GetResidues();
            /*! \fn
              * An accessor function in order to access to the number of residues in a pdb site
              * @return number_of_residues_ attribute of the current object of this class
              */
            int GetNumberOfResidues();
/** @}*/
            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
/** \addtogroup Manipulators
               * @{
               */
            /*! \fn
              * A mutator function in order to set the site id of the current object
              * Set the site_id_ attribute of the current site
              * @param site_id The site id of the current object
              */
            void SetSiteId(const std::string site_id);
            /*! \fn
              * A mutator function in order to set the residues of the current object
              * Set the residues_ attribute of the current site
              * @param residues The residues of the current object
              */
            void SetResidues(SiteResidueVector residues);
            /*! \fn
              * A mutator function in order to add the residue to the current object
              * Set the residues_ attribute of the current site
              * @param residue The residue of the current object
              */
            void AddResidue(PdbSiteResidue* residue);
            /*! \fn
              * A mutator function in order to set the number of residues of the current object
              * Set the number_of_residues_ attribute of the current site
              * @param number_of_residues The number of residues of the current object
              */
            void SetNumberOfResidues(int number_of_residues);
/** @}*/
            //////////////////////////////////////////////////////////
            //                        FUNCTIONS                     //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the pdb site contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cerr);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            std::string site_id_;           /*!< Site identifier that appears in a site card for each site >*/
            SiteResidueVector residues_;    /*!< Residues that are in the site card of a pdb file >*/
            int number_of_residues_;        /*!< Number of residues in the current site object >*/
    };
}


#endif // PDBSITECARD_HPP
