// Created by: Delaram Rahbarinia
// Modified by: Alireza Khatamian, Delaram Rahbarinia

#ifndef PDBLINK_HPP
#define PDBLINK_HPP

#include <string>
#include <vector>
#include <iostream>

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
            /*! \fn
              * Default constructor
              */
            PdbLink();
            PdbLink(std::string& line);

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
            /*! \fn
              * An accessor function in order to access to the residues in a link class
              * @return residues_ attribute of the current object of this class
              */
            LinkResidueVector GetResidues();
            /*! \fn
              * An accessor function in order to access to the link length in a link class
              * @return link_length_ attribute of the current object of this class
              */
            double GetLinkLength();

            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A mutator function in order to set the residues of the current object
              * Set the residues_ attribute of the current link card
              * @param residues The residues of the current object
              */
            void SetResidues(const LinkResidueVector residues);
            /*! \fn
              * A function in order to add residue to the current object
              * Set the residues_ attribute of the current link card
              * @param residue The residue of the current object
              */
            void AddResidue(PdbLinkResidue* residue);
            /*! \fn
              * A mutator function in order to set the link length of the current object
              * Set the link_length_ attribute of the current link card
              * @param link_length The link length of the current object
              */
            void SetLinkLength(double link_length);

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
            LinkResidueVector residues_;
            double link_length_;

    };
}

#endif // PDBLINK_HPP
