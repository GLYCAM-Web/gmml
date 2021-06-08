// Created by: Dave Montgomery

#ifndef PDBCISPEPTIDESECTION_HPP
#define PDBCISPEPTIDESECTION_HPP

#include <string>
#include <vector>
#include <iostream>

namespace PdbFileSpace
{
    class PdbCISPeptideCard;

    class PdbCISPeptideSection
    {
        public:
            //////////////////////////////////////////////////////////
            //                       TYPE DEFINITION                //
            //////////////////////////////////////////////////////////
            /*! \typedef
              * List of CIS peptide cards
              */
            typedef std::vector<PdbCISPeptideCard*> CISPeptideCardVector;


            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            PdbCISPeptideSection();
            /*! \fn
              * Constructor with required parameters
              * @param stream_block
              */
            PdbCISPeptideSection(std::stringstream& stream_block);

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
/** \addtogroup Molecular_Data_Structure
              * @{
              */
            /*! \fn
              * An accessor function in order to access to the cis_peptide in a cis_peptide card
              * @return cis_peptide_ attribute of the current object of this class
              */
            CISPeptideCardVector GetCISPeptideCards();
/** @}*/
            //////////////////////////////////////////////////////////
            //                       MUTATOR                 PdbCISPeptideSection       //
            //////////////////////////////////////////////////////////
/** \addtogroup Manipulators
              * @{
              */
            /*! \fn
              * A mutator function in order to set the cis_peptide attribute of the current object
              * Set the cis_peptide_ attribute of the current cis_peptide card
              * @param cis_peptide The cis_peptide attribute of the current object
              */
            void SetCISPeptideCards(CISPeptideCardVector cis_peptide);
            /*! \fn
              * A function in order to add the cis_peptide attribute to the current object
              * Set the cis_peptide_ attribute of the current cis_peptide card
              * @param cis_peptide The cis_peptide attribute of the current object
              */
            void AddCISPeptideCards(PdbCISPeptideCard *cis_peptide);
/** @}*/
            //////////////////////////////////////////////////////////
            //                        FUNCTIONS                     //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the cis_peptide section contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cerr);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            CISPeptideCardVector cis_peptide_;      /*!< List of CIS peptide cards >*/

    };
}

#endif // PDBCISPEPTIDESECTION_HPP
