#ifndef PDBSITE_HPP
#define PDBSITE_HPP

#include <string>
#include <vector>

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

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
            int GetSiteId();
            SiteResidueVector GetResidues();
            int GetNumberOfResidues();

            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
            void SetSiteId(int site_id);
            void SetResidues(const SiteResidueVector residues);
            void SetNumberOfResidues(int number_of_residues);

            //////////////////////////////////////////////////////////
            //                        FUNCTIONS                     //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            int site_id_;
            SiteResidueVector residues_;
            int number_of_residues_;
    };
}


#endif // PDBSITE_HPP
