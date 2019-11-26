// Created by: Delaram Rahbarinia
// Modified by: Alireza Khatamian, Delaram Rahbarinia, Dave Montgomery

#ifndef PDBLINKCARD_HPP
#define PDBLINKCARD_HPP

#include <string>
#include <vector>
#include <iostream>

namespace PdbFileSpace
{
    class PdbLinkCardResidue;

    class PdbLinkCard
    {
        public:
            //////////////////////////////////////////////////////////
            //                       TYPE DEFINITION                //
            //////////////////////////////////////////////////////////
            /*! \typedef
              * List of residues in a link
              */
            typedef std::vector< PdbLinkCardResidue* > LinkResidueVector;

            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            PdbLinkCard();
            /*! \fn
              * Constructor with required parameters
              * @param line
              */
            PdbLinkCard(std::string& line);

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
/** \addtogroup Molecular_Data_Structure
              * @{
              */
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
/** @}*/
            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
/** \addtogroup Manipulators
              * @{
              */
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
            void AddResidue(PdbLinkCardResidue* residue);
            /*! \fn
              * A mutator function in order to set the link length of the current object
              * Set the link_length_ attribute of the current link card
              * @param link_length The link length of the current object
              */
            void SetLinkLength(double link_length);
/** @}*/
            //////////////////////////////////////////////////////////
            //                        FUNCTIONS                     //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the link contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cerr);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            LinkResidueVector residues_;        /*!< List of residues in a link >*/
            double link_length_;                /*!< Link length >*/

    };
}

#endif // PDBLINKCARD_HPP
