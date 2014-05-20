// Created by: Delaram Rahbarinia
// Modified by: Alireza Khatamian, Delaram Rahbarinia

#ifndef PDBSITE_HPP
#define PDBSITE_HPP

#include <string>
#include <sstream>
#include <vector>
#include <iostream>

namespace PdbFileSpace
{
    class PdbSiteResidue;

    class PdbSite
    {
        public:
            //////////////////////////////////////////////////////////
            //                       TYPE DEFINITION                //
            //////////////////////////////////////////////////////////
            typedef std::vector< PdbSiteResidue* > SiteResidueVector;

            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            PdbSite();
            PdbSite(std::stringstream& site_block);

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
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

            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
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

            //////////////////////////////////////////////////////////
            //                        FUNCTIONS                     //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            void Print(std::ostream& out = std::cout);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            std::string site_id_;
            SiteResidueVector residues_;
            int number_of_residues_;
    };
}


#endif // PDBSITE_HPP
