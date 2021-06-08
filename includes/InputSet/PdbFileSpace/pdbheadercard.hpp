// Created by: Alireza Khatamian
// Modified by: Alireza Khatamian, Delaram Rahbarinia

#ifndef PDBHEADERCARD_HPP
#define PDBHEADERCARD_HPP

#include <string>
#include <sstream>
#include <iostream>

namespace PdbFileSpace
{
    class PdbHeaderCard
    {
        public:
            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            PdbHeaderCard();
            /*! \fn
              * Constructor with required parameters
              * @param record_name Record name of header card appears in the first column of each line in a pdb file
              * @param classification Classification of the pdb file
              * @param deposition_date Date that the file has been deposited
              * @param identifier_code Identifier code of the pdb file
              */
            PdbHeaderCard(const std::string& record_name, const std::string& classification, const std::string& deposition_date, const std::string& identifier_code);
            /*! \fn
              * Constructor with required parameters
              * @param stream_block
              */
            PdbHeaderCard(std::stringstream& stream_block);

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
/** \addtogroup Molecular_Data_Structure
              * @{
              */
            /*! \fn
              * An accessor function in order to access to the record name in a header card
              * @return record_name_ attribute of the current object of this class
              */
            std::string GetRecordName();
            /*! \fn
              * An accessor function in order to access to the classification attribute of a header card
              * @return classification_ attribute of the current object of this class
              */
            std::string GetClassification();
            /*! \fn
              * An accessor function in order to access to the deposition date in a header card
              * @return deposition_date_ attribute of the current object of this class
              */
            std::string GetDepositionDate();
            /*! \fn
              * An accessor function in order to access to the identifier code in a header card
              * @return identifier_code_ attribute of the current object of this class
              */
            std::string GetIdentifierCode();
/** @}*/
            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
/** \addtogroup Manipulators
              * @{
              */
            /*! \fn
              * A mutator function in order to set the record name of the current object
              * Set the record_name_ attribute of the current header card
              * @param record_name The record name of the current object
              */
            void SetRecordName(const std::string record_name);
            /*! \fn
              * A mutator function in order to set the classification attribute of the current object
              * Set the classification_ attribute of the current header card
              * @param classification The classification of the current object
              */
            void SetClassification(const std::string classification);
            /*! \fn
              * A mutator function in order to set the deposition date of the current object
              * Set the deposition_date_ attribute of the current header card
              * @param deposition_date The rdeposition date of the current object
              */
            void SetDepositionDate(const std::string deposition_date);
            /*! \fn
              * A mutator function in order to set the identifier code of the current object
              * Set the identifier_code_ attribute of the current header card
              * @param identifier_code The identifier code of the current object
              */
            void SetIdentificationCode(const std::string identifier_code);
/** @}*/
            //////////////////////////////////////////////////////////
            //                        FUNCTIONS                     //
            //////////////////////////////////////////////////////////            

            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the header card contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cerr);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            std::string record_name_;           /*!< Record name of headr card in a pdb file >*/
            std::string classification_;        /*!< Classification of the pdb file >*/
            std::string deposition_date_;       /*!< Date of deposition >*/
            std::string identifier_code_;       /*!< Identifier code of the pdb file >*/

    };
}
#endif // PDBHEADERCARD_HPP
