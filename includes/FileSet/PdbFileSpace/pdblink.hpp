#ifndef PDBLINK_HPP
#define PDBLINK_HPP

#include <string>
#include <vector>

namespace PdbFileSpace
{
    class PdbLinkResidue;

    class PdbLink
    {
        public:
            //////////////////////////////////////////////////////////
            //                       TYPE DEFINITION                //
            //////////////////////////////////////////////////////////
            typedef std::vector< PdbLinkResidue* > LinkResidueVector;

            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            PdbLink();
            PdbLink(std::string& line);

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
            LinkResidueVector GetResidues();
            double GetLinkLength();

            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
            void SetResidues(const LinkResidueVector residues);
            void SetLinkLength(double link_length);

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
            LinkResidueVector residues_;
            double link_length_;

    };
}

#endif // PDBLINK_HPP
