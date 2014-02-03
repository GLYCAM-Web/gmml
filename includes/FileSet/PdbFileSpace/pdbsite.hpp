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
            PdbSite();
            PdbSite(std::stringstream& site_block);

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
            std::string GetSiteId();
            SiteResidueVector GetResidues();
            int GetNumberOfResidues();

            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
            void SetSiteId(const std::string site_id);
            void SetResidues(const SiteResidueVector residues);
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
